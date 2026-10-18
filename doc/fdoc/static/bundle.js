/* fdoc bundle mode: show one of a module's automatic function pages,
   selected by the URL fragment "#name" or "#name:symbol". */
(function () {
  "use strict";
  var B = window.FDOC_BUNDLE;
  if (!B) return;
  var side = document.getElementById("sidebar");
  var main = document.getElementById("main");

  function show() {
    var h = decodeURIComponent(location.hash.slice(1));
    var i = h.indexOf(":");
    var name = i >= 0 ? h.slice(0, i) : h, sym = i >= 0 ? h.slice(i + 1) : "";
    var p = B.pages[name];
    if (!p) {
      for (var k in B.pages) { p = B.pages[k]; name = k; break; }
      if (!p) return;
    }
    document.title = p[0];
    side.innerHTML = B.sidebars[p[1]];
    main.innerHTML = p[2];
    /* in-page anchors (#x) refer to the page shown: make them #page:x */
    var inpage = main.querySelectorAll('a[href^="#"]');
    for (var q = 0; q < inpage.length; q++) {
      var a = inpage[q].getAttribute("href").slice(1);
      if (a && a.indexOf(":") < 0) inpage[q].setAttribute("href", "#" + name + ":" + a);
    }
    var here = location.href.split("#")[0];
    var links = side.querySelectorAll("a");
    for (var j = 0; j < links.length; j++) {
      if (links[j].href === here + "#" + name) {
        links[j].setAttribute("aria-current", "page");
        var top = links[j].getBoundingClientRect().top - side.getBoundingClientRect().top + side.scrollTop;
        if (top > side.clientHeight * 0.6 || top < side.scrollTop) side.scrollTop = top - side.clientHeight / 2;
      }
    }
    var target = sym && document.getElementById(sym);
    if (target) target.scrollIntoView(); else window.scrollTo(0, 0);
  }
  window.addEventListener("hashchange", show);
  show();
})();
