// Make horizontally scrolling regions reachable without a mouse.
//
// Wide tables and block equations are put in their own scroll box so they do
// not push the page sideways. A box that scrolls but holds nothing focusable
// cannot be scrolled by keyboard at all, so on a narrow screen its right-hand
// side is simply unreachable -- WCAG 2.1.1, and on the evidence page that is
// the verdict column. Giving it tabindex="0" makes it a tab stop the arrow
// keys then scroll; role and label stop it being an unnamed landmark.
//
// Done here rather than in the Rmd because it has to apply to every generated
// table on every page, including the reference index, and because whether a
// region actually overflows is only knowable once it is laid out.
(function () {
  "use strict";
  function label(el) {
    var cap = el.querySelector("caption");
    if (cap && cap.textContent.trim()) return cap.textContent.trim();
    var h = el.closest("section, .section");   // pkgdown emits div.section
    h = h && h.querySelector("h1, h2, h3");
    return h && h.textContent.trim()
      ? "Scrollable content: " + h.textContent.trim()
      : "Scrollable content";
  }

  function mark() {
    var sel = "main table, main .math.display, main math[display='block'], main pre";
    document.querySelectorAll(sel).forEach(function (el) {
      var overflows = el.scrollWidth > el.clientWidth + 1;
      if (!overflows) {
        // a previously marked region can stop overflowing when the window grows
        if (el.dataset.scrollRegion === "1") {
          el.removeAttribute("tabindex");
          el.removeAttribute("role");
          el.removeAttribute("aria-label");
          delete el.dataset.scrollRegion;
        }
        return;
      }
      if (el.dataset.scrollRegion === "1") return;
      // don't steal focus from something already focusable inside
      if (el.querySelector("a[href], button, input, select, textarea")) return;
      el.setAttribute("tabindex", "0");
      el.setAttribute("role", "region");
      el.setAttribute("aria-label", label(el));
      el.dataset.scrollRegion = "1";
    });
  }

  if (document.readyState === "loading") {
    document.addEventListener("DOMContentLoaded", mark);
  } else {
    mark();
  }
  // images load late and can change what overflows; so can a resize
  window.addEventListener("load", mark);
  var t;
  window.addEventListener("resize", function () {
    clearTimeout(t);
    t = setTimeout(mark, 150);
  });
})();
