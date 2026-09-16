"""Derive every logo asset from pkgdown/logo-source.png.

    python pkgdown/make-logo.py

Writes man/figures/logo.png and the four files in pkgdown/favicon/. The source
is never modified. Requires Pillow; nothing in the package or the site build
depends on this script, it is only here so the assets can be rebuilt.

Two corrections happen on the way, both worth recording because both were
invisible until measured.

CROPPING. The source carries a handful of pixels at alpha 1-15 out at the
canvas edge -- invisible against any background, but enough that
Image.getbbox(), which is alpha > 0, reports the hexagon as touching the edge
when the mid-line is actually transparent for 90px inward. Cropping on that
left 12% dead width in the file; since the README sizes the image by width,
everything drawn inside came out 12% small. Every threshold from 16 to 250
agrees on the hexagon to within two pixels, so ALPHA_FLOOR is not a tuned
number.

GEOMETRY. The source is a generated image, so its hexagon is drawn rather than
constructed, and it is not quite regular:

  * 0.94% too wide for its height -- 1077px where a regular hexagon of that
    height is 1067 -- so the two vertical sides sit ~5px outside where they
    belong while the top and bottom vertices are correct;
  * the lower sloping edges run at 60.25 degrees against 59.92 for the upper
    pair, so the shape is very slightly keystoned rather than mirror-symmetric.

Together: outline deviation of rms 3.4px, worst 6.8px, on a 616px circumradius.
Under 1%, but a constant-width rim running round a regular polygon turns a
shape error into a visible thick-and-thin, and it read as a wonky hexagon.

Rescaling each row horizontally so its span matches the exact hexagon and its
midpoint sits on the axis removes the stretch and the keystone in one pass. The
median correction is 0.94%, so nothing in the artwork visibly moves, and the
deviation drops to rms 0.7px.
"""

import base64
import io
import math
import os

from PIL import Image

ALPHA_FLOOR = 16          # see CROPPING above
WEB_WIDTH = 560           # ~2x the 272px the README draws it at on GitHub

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "pkgdown", "logo-source.png")
LOGO = os.path.join(ROOT, "man", "figures", "logo.png")
FAVICON = os.path.join(ROOT, "pkgdown", "favicon")


def cropped():
    im = Image.open(SRC).convert("RGBA")
    mask = im.split()[3].point(lambda p: 255 if p >= ALPHA_FLOOR else 0)
    return im.crop(mask.getbbox())


def squared(im):
    """Map the drawn hexagon onto an exact one of the same height."""
    W, H = im.size
    alpha = im.split()[3].point(lambda p: 255 if p >= ALPHA_FLOOR else 0).load()
    R = H / 2.0
    IW = int(round(math.sqrt(3) * R))

    def ideal_half(y):
        dy = abs(y + 0.5 - H / 2.0)
        if dy <= R / 2:
            return math.sqrt(3) * R / 2
        return math.sqrt(3) * max(R - dy, 0.0)

    out = Image.new("RGBA", (IW, H), (0, 0, 0, 0))
    for y in range(H):
        xs = [x for x in range(W) if alpha[x, y]]
        if not xs:
            continue
        iw = int(round(2 * ideal_half(y)))
        if iw < 1:
            continue
        strip = im.crop((xs[0], y, xs[-1] + 1, y + 1)).resize((iw, 1), Image.LANCZOS)
        out.paste(strip, (int(round(IW / 2 - iw / 2)), y))
    print("  geometry: %dx%d drawn -> %dx%d exact" % (W, H, IW, H))
    return out


def write_favicons(im):
    os.makedirs(FAVICON, exist_ok=True)
    side = max(im.size)
    sq = Image.new("RGBA", (side, side), (0, 0, 0, 0))
    sq.paste(im, ((side - im.width) // 2, (side - im.height) // 2))

    for size, name in ((96, "favicon-96x96.png"), (180, "apple-touch-icon.png")):
        p = os.path.join(FAVICON, name)
        sq.resize((size, size), Image.LANCZOS).save(p, "PNG", optimize=True)
        print("  %-22s %d bytes" % (name, os.path.getsize(p)))

    p = os.path.join(FAVICON, "favicon.ico")
    sq.resize((64, 64), Image.LANCZOS).save(
        p, "ICO", sizes=[(16, 16), (32, 32), (48, 48), (64, 64)])
    print("  %-22s %d bytes" % ("favicon.ico", os.path.getsize(p)))

    # pkgdown's head template links favicon.svg unconditionally, so a raster
    # wrapped in SVG is here only to stop that being a 404 on every page.
    buf = io.BytesIO()
    sq.resize((96, 96), Image.LANCZOS).save(buf, "PNG", optimize=True)
    svg = ('<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 96 96" '
           'width="96" height="96"><title>fleet</title>'
           '<image width="96" height="96" href="data:image/png;base64,%s"/></svg>\n'
           % base64.b64encode(buf.getvalue()).decode("ascii"))
    p = os.path.join(FAVICON, "favicon.svg")
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(svg)
    print("  %-22s %d bytes" % ("favicon.svg", os.path.getsize(p)))


def main():
    im = cropped()
    print("source %s" % os.path.basename(SRC))
    hexagon = squared(im)
    web = hexagon.resize(
        (WEB_WIDTH, int(round(WEB_WIDTH * hexagon.height / hexagon.width))),
        Image.LANCZOS)
    web.save(LOGO, "PNG", optimize=True)
    print("  %-22s %s  %d bytes"
          % ("man/figures/logo.png", web.size, os.path.getsize(LOGO)))
    write_favicons(web)


if __name__ == "__main__":
    main()
