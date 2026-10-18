/* fdoc: search, theme toggle, mobile navigation. No dependencies. */
(function () {
  "use strict";
  var root = window.FDOC_ROOT || "";
  var index = window.FDOC_INDEX || [];

  /* theme */
  var themeBtn = document.getElementById("theme-toggle");
  if (themeBtn) themeBtn.addEventListener("click", function () {
    var cur = document.documentElement.dataset.theme ||
      (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
    var next = cur === "dark" ? "light" : "dark";
    document.documentElement.dataset.theme = next;
    try { localStorage.setItem("fdoc-theme", next); } catch (e) {}
  });

  /* mobile sidebar */
  var sbBtn = document.getElementById("sb-toggle");
  if (sbBtn) sbBtn.addEventListener("click", function () {
    var open = document.body.classList.toggle("sb-open");
    sbBtn.setAttribute("aria-expanded", open ? "true" : "false");
  });
  document.addEventListener("click", function (e) {
    if (document.body.classList.contains("sb-open") && !e.target.closest(".sidebar, .sb-toggle")) {
      document.body.classList.remove("sb-open");
      if (sbBtn) sbBtn.setAttribute("aria-expanded", "false");
    }
  });
  /* bring the current entry into view inside the sidebar only */
  var sb = document.getElementById("sidebar");
  var cur = sb && sb.querySelector("[aria-current]");
  if (cur) {
    var top = cur.getBoundingClientRect().top - sb.getBoundingClientRect().top + sb.scrollTop;
    if (top > sb.clientHeight * 0.6) sb.scrollTop = top - sb.clientHeight / 2;
  }

  /* "Cite this page": switch between the text and BibTeX forms, copy.
     Delegated, since bundle mode replaces the page content. */
  document.addEventListener("click", function (e) {
    var box = e.target.closest && e.target.closest(".cite-page");
    if (!box) return;
    var tab = e.target.closest(".cite-tab");
    if (tab) {
      var which = tab.getAttribute("data-cite");
      Array.prototype.forEach.call(box.querySelectorAll(".cite-tab"), function (t) {
        t.setAttribute("aria-selected", t === tab ? "true" : "false");
      });
      Array.prototype.forEach.call(box.querySelectorAll("[data-cite-panel]"), function (p) {
        p.hidden = p.getAttribute("data-cite-panel") !== which;
      });
      return;
    }
    var btn = e.target.closest(".cite-copy");
    if (!btn) return;
    var panel = box.querySelector("[data-cite-panel]:not([hidden])");
    if (!panel) return;
    var text = panel.textContent;
    function done(ok) {
      btn.textContent = ok ? "Copied" : "Select and copy";
      setTimeout(function () { btn.textContent = "Copy"; }, 1600);
    }
    function fallback() {
      var range = document.createRange();
      range.selectNodeContents(panel);
      var sel = window.getSelection();
      sel.removeAllRanges(); sel.addRange(range);
      var ok = false;
      try { ok = document.execCommand("copy"); } catch (err) {}
      done(ok);
    }
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(function () { done(true); }, fallback);
    } else {
      fallback();
    }
  });

  /* search */
  var input = document.getElementById("search");
  var list = document.getElementById("search-results");
  if (!input || !list) return;
  var sel = -1, results = [];

  function score(entry, q, words) {
    var name = entry[0].toLowerCase();
    var isSym = entry[2] !== "algorithm" && entry[2] !== "topic" && entry[2] !== "module" && entry[1].indexOf("#") >= 0 || /^[a-z_0-9]+$/.test(entry[0]);
    if (name === q) return 1000;
    if (name.indexOf(q) === 0) return 800 - name.length;
    var i = name.indexOf(q);
    if (i >= 0) return 600 - i - name.length / 4;
    // all words match the name (e.g. "mulmod preinv")
    if (words.length > 1 && words.every(function (w) { return name.indexOf(w) >= 0; })) return 450 - name.length / 4;
    var hay = (entry[0] + " " + entry[3]).toLowerCase();
    if (words.every(function (w) { return hay.indexOf(w) >= 0; })) return (isSym ? 200 : 250) - name.length / 8;
    return -1;
  }

  function esc(s) {
    return s.replace(/[&<>"]/g, function (c) { return { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]; });
  }

  function render() {
    var q = input.value.trim().toLowerCase();
    if (!q) { list.hidden = true; return; }
    var words = q.split(/\s+/);
    results = [];
    for (var i = 0; i < index.length; i++) {
      var s = score(index[i], q, words);
      if (s >= 0) results.push([s, index[i]]);
    }
    results.sort(function (a, b) { return b[0] - a[0]; });
    results = results.slice(0, 14).map(function (r) { return r[1]; });
    sel = results.length ? 0 : -1;
    if (!results.length) {
      list.innerHTML = '<li class="r-none">No matches for “' + esc(input.value.trim()) + '”</li>';
    } else {
      list.innerHTML = results.map(function (e, j) {
        var page = !/^[A-Za-z_][A-Za-z_0-9]*$/.test(e[0]);
        return '<li role="option" aria-selected="' + (j === sel) + '"><a href="' + root + e[1] + '">' +
          '<span class="r-name' + (page ? " page" : "") + '">' + esc(e[0]) + '</span>' +
          '<span class="r-kind">' + esc(e[2]) + '</span>' +
          (e[3] ? '<span class="r-desc">' + esc(e[3]) + '</span>' : "") + '</a></li>';
      }).join("");
    }
    list.hidden = false;
  }

  function move(d) {
    if (!results.length) return;
    sel = (sel + d + results.length) % results.length;
    Array.prototype.forEach.call(list.children, function (li, j) { li.setAttribute("aria-selected", j === sel); });
    var li = list.children[sel];
    if (li && li.scrollIntoView) li.scrollIntoView({ block: "nearest" });
  }

  input.addEventListener("input", render);
  input.addEventListener("focus", render);
  input.addEventListener("keydown", function (e) {
    if (e.key === "ArrowDown") { e.preventDefault(); move(1); }
    else if (e.key === "ArrowUp") { e.preventDefault(); move(-1); }
    else if (e.key === "Enter" && sel >= 0) { e.preventDefault(); location.href = root + results[sel][1]; }
    else if (e.key === "Escape") { list.hidden = true; input.blur(); }
  });
  document.addEventListener("keydown", function (e) {
    if (e.key === "/" && document.activeElement !== input && !/INPUT|TEXTAREA/.test(document.activeElement.tagName)) {
      e.preventDefault(); input.focus(); input.select();
    }
  });
  document.addEventListener("click", function (e) {
    if (!e.target.closest(".search")) list.hidden = true;
  });
})();
