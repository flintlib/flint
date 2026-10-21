"""
HTML backend.  Math is pre-rendered with KaTeX (one batched node process,
cached between builds); if node is not available the pages fall back to
client-side rendering with the same vendored KaTeX.
"""

import html
import json
import os
import re
import shutil
import subprocess

from .markup import Node, plain_text
from .csyntax import C_KEYWORDS, C_TYPES, tokenize
from .outputs import describe as describe_output
from .model import THEOREM_KINDS, COLLECTIONS, PROSE_KINDS, anchor_url

HERE = os.path.dirname(os.path.abspath(__file__))
KATEX_DIR = os.path.join(HERE, "vendor", "katex")

esc = html.escape


# ---------------------------------------------------------------------------
# KaTeX
# ---------------------------------------------------------------------------

KATEX_JS = r"""
const katex = require(process.argv[2]);
let data = '';
process.stdin.on('data', d => data += d);
process.stdin.on('end', () => {
  const req = JSON.parse(data);
  const out = req.items.map(([tex, display]) => {
    try {
      return [true, katex.renderToString(tex, {displayMode: display, macros: Object.assign({}, req.macros),
                                               throwOnError: true, strict: 'ignore'})];
    } catch (e) {
      return [false, String(e.message || e)];
    }
  });
  process.stdout.write(JSON.stringify(out));
});
"""


class Katex:
    def __init__(self, cache_path, macros, warn):
        self.items = {}
        self.cache_path = cache_path
        self.macros = macros
        self.warn = warn
        self.cache = {}
        self.available = shutil.which("node") is not None
        if os.path.exists(cache_path):
            try:
                with open(cache_path) as f:
                    self.cache = json.load(f)
            except Exception:
                self.cache = {}
        self.where = {}

    def expand(self, tex):
        """Expand the user macros (macros.tex) textually, so that the stored
        TeX, the MathML annotation and any client-side fallback are all
        self-contained."""
        for _ in range(5):
            new = tex
            for name, body in self.macros.items():
                new = re.sub(re.escape(name) + r"(?![A-Za-z])", lambda m: body, new)
            if new == tex:
                break
            tex = new
        return tex

    def register(self, key, where=""):
        """The index of a formula given by its key (as made by placeholder)."""
        if key not in self.items:
            self.items[key] = len(self.items)
            if where:
                self.where.setdefault(key, where)
        return self.items[key]

    def placeholder(self, tex, display, where=""):
        tex = self.expand(tex)
        key = ("D" if display else "I") + tex
        if key not in self.items:
            self.items[key] = len(self.items)
            self.where[key] = where
        return "\x00K%d\x00" % self.items[key]

    def render_all(self):
        macro_sig = json.dumps(self.macros, sort_keys=True)
        if self.cache.get("__macros__") != macro_sig:
            self.cache = {"__macros__": macro_sig}
        todo = [k for k in self.items if k not in self.cache]
        if todo and self.available:
            script = os.path.join(os.path.dirname(self.cache_path), ".katex_batch.js")
            with open(script, "w") as f:
                f.write(KATEX_JS)
            req = {"macros": self.macros, "items": [[k[1:], k[0] == "D"] for k in todo]}
            res = subprocess.run(["node", script, os.path.join(KATEX_DIR, "katex.min.js")],
                                 input=json.dumps(req).encode(), capture_output=True)
            if res.returncode != 0:
                self.warn("KaTeX failed: " + res.stderr.decode()[:500])
                self.available = False
            else:
                for k, (ok, out) in zip(todo, json.loads(res.stdout)):
                    if ok:
                        self.cache[k] = out
                    else:
                        self.warn("math error in %s: %s" % (self.where.get(k, "?"), out))
                        self.cache[k] = ('<span class="math-error" title="%s">%s</span>'
                                         % (esc(out), esc(k[1:])))
                with open(self.cache_path, "w") as f:
                    json.dump(self.cache, f)
        self.out = {}
        for k, i in self.items.items():
            if k in self.cache:
                self.out[i] = self.cache[k]
            else:   # client-side fallback
                tag = "div" if k[0] == "D" else "span"
                self.out[i] = '<%s class="tex%s">%s</%s>' % (tag, " display" if k[0] == "D" else "",
                                                             esc(k[1:]), tag)

    @property
    def prerendered(self):
        return all(k in self.cache for k in self.items)

    def substitute(self, s):
        return re.sub("\x00K(\\d+)\x00", lambda m: self.out[int(m.group(1))], s)


# ---------------------------------------------------------------------------
# Renderer for one page
# ---------------------------------------------------------------------------

class HtmlRenderer:
    def __init__(self, site, page, katex):
        self.site = site
        self.page = page
        self.katex = katex
        self.depth = page.url.count("/")
        self.cites = []

    def rel(self, url):
        if url is None:
            return "#"
        if "://" in url:
            return url
        base, _, anchor = url.partition("#")
        if base == self.page.url and anchor:
            return "#" + anchor
        return "../" * self.depth + url

    def where(self, node=None):
        line = node.get("line") if node is not None else None
        return "%s%s" % (os.path.relpath(self.page.path or self.page.id, self.site.srcroot),
                         ":%d" % line if line else "")

    # -- inline ----------------------------------------------------------------

    def inline(self, nodes):
        out = []
        for nd in nodes:
            k = nd.kind
            if k == "text":
                out.append(esc(nd["text"], quote=False))
            elif k == "code":
                out.append(self.code_span(nd["text"]))
            elif k == "math":
                out.append(self.katex.placeholder(nd["tex"], False, self.where()))
            elif k == "em":
                out.append("<em>%s</em>" % self.inline(nd.children))
            elif k == "strong":
                out.append("<strong>%s</strong>" % self.inline(nd.children))
            elif k == "cite":
                out.append(self.cite(nd["keys"], nd.get("locs")))
            elif k == "link":
                out.append(self.link(nd))
            else:
                out.append(self.inline(nd.children))
        return "".join(out)

    def plain_inline(self, nodes):
        """Inline content without links, for use inside another link (cards)."""
        self.nolinks = True
        try:
            return self.inline(nodes)
        finally:
            self.nolinks = False

    def code_span(self, text):
        name = text[:-2] if text.endswith("()") else text
        sym = self.site.symbols.get(name)
        if sym is not None and not getattr(self, "nolinks", False):
            cls = "sym self" if sym.block.page is self.page else "sym"
            return '<a class="%s" href="%s"><code>%s</code></a>' % (cls, self.rel(sym.url), esc(text))
        return "<code>%s</code>" % esc(text)

    def link(self, nd):
        target = nd["target"]
        if getattr(self, "nolinks", False):
            return self.inline(nd.children) if nd.children else esc(target)
        url, default, kind = self.site.resolve(target, self.page)
        if kind == "unresolved":
            self.site.warn("unresolved link target %r" % target, self.where())
            return '<span class="broken">%s</span>' % (self.inline(nd.children) or esc(target))
        text = self.inline(nd.children) if nd.children else (self.inline(default) if default else esc(target))
        if kind == "external":
            cls = ' class="ext"' if not nd.get("bare") else ' class="ext url"'
            return '<a%s href="%s">%s</a>' % (cls, esc(url), text)
        return '<a href="%s">%s</a>' % (self.rel(url), text)

    def cite(self, keys, locs=None):
        parts = []
        locs = locs or [""] * len(keys)
        for key, loc in zip(keys, locs):
            if key not in self.site.bib:
                self.site.warn("unknown citation [@%s]" % key, self.where())
                parts.append('<span class="broken">%s</span>' % esc(key))
                continue
            if key not in self.cites:
                self.cites.append(key)
            title = plain_text(self.site.bib[key])
            parts.append('<a href="#ref-%s" title="%s">%s</a>%s' % (esc(key), esc(title), esc(key),
                                                                  ", " + esc(loc) if loc else ""))
        return '<span class="cite">[%s]</span>' % "; ".join(parts)

    # -- blocks ------------------------------------------------------------------

    def blocks(self, nodes):
        out = []
        run = []        # consecutive '::: summary' cards, laid out as one grid
        for nd in nodes:
            if nd.kind == "div" and nd.get("kind") == "summary":
                card = self.block(nd)
                if card:
                    run.append(card)
                continue
            if run:
                out.append(self.card_grid(run))
                run = []
            out.append(self.block(nd))
        if run:
            out.append(self.card_grid(run))
        return "\n".join(out)

    @staticmethod
    def card_grid(cards):
        if len(cards) == 1:
            return cards[0]
        return '<div class="cards">%s</div>' % "".join(cards)

    def block(self, nd):
        k = nd.kind
        if k == "para":
            return "<p>%s</p>" % self.inline(nd.children)
        if k == "heading":
            lvl = min(nd["level"] + (0 if self.page.kind in ("module", "home", "refs") else 0), 6)
            hid = nd["id"] or slug(plain_text(nd.children))
            return ('<h%d id="%s">%s<a class="hl" href="#%s" aria-label="Link to this section">#</a></h%d>'
                    % (lvl, hid, self.inline(nd.children), hid, lvl))
        if k == "list":
            tag = "ol" if nd["ordered"] else "ul"
            start = ' start="%d"' % nd["start"] if nd["ordered"] and nd["start"] != 1 else ""
            items = []
            for it in nd["items"]:
                if len(it) == 1 and it[0].kind == "para":
                    items.append("<li>%s</li>" % self.inline(it[0].children))
                else:
                    items.append("<li>%s</li>" % self.blocks(it))
            return "<%s%s>%s</%s>" % (tag, start, "".join(items), tag)
        if k == "displaymath":
            return '<div class="dmath">%s</div>' % self.katex.placeholder(nd["tex"], True, self.where(nd))
        if k == "code":
            return self.code_block(nd)
        if k == "table":
            return self.table(nd)
        if k == "div":
            return self.div(nd)
        return ""

    def table(self, nd):
        al = nd["aligns"]
        def row(cells, tag):
            return "<tr>%s</tr>" % "".join(
                '<%s%s>%s</%s>' % (tag, ' class="r"' if i < len(al) and al[i] == "right" else
                                   ' class="c"' if i < len(al) and al[i] == "center" else "",
                                   self.inline(c), tag) for i, c in enumerate(cells))
        return ('<div class="table-wrap"><table><thead>%s</thead><tbody>%s</tbody></table></div>'
                % (row(nd["head"], "th"), "".join(row(r, "td") for r in nd["rows"])))

    def code_block(self, nd):
        lang = nd["lang"]
        code = self.highlight(nd["text"]) if lang in ("c", "C", "h") else esc(nd["text"])
        if nd.get("example_name"):
            out = ""
            if nd.get("output") is not None:
                note = describe_output(nd.get("output_flags"))
                out = ('<div class="ex-out"><div class="ex-label">Output%s</div><pre>%s</pre></div>'
                       % ('<span class="ex-note">%s</span>' % esc(note) if note else "", esc(nd["output"])))
            req = nd["flags"].get("requires")
            reqs = ('<span class="ex-req">%s</span>' % esc(", ".join(
                        ("%s-bit only" % r) if r.isdigit() else ("needs %s" % r)
                        for r in (x.strip() for x in str(req).split(",")) if r)) if req else "")
            return ('<figure class="example" id="example-%s"><figcaption><span class="ex-title">Example</span>'
                    '<span class="ex-file">%s.c</span>%s<span class="ex-check" title="This program is '
                    'compiled and run by make checkdoc; its output must match exactly.">checked by make checkdoc</span>'
                    '</figcaption><pre class="code"><code>%s</code></pre>%s</figure>'
                    % (esc(nd["example_name"]), esc(nd["example_name"]), reqs, code, out))
        return '<pre class="code"><code>%s</code></pre>' % code

    def highlight(self, src, link_self=True):
        out = []
        for kind, tok in tokenize(src):
            if kind == "comment":
                out.append('<span class="c-com">%s</span>' % esc(tok))
            elif kind == "pre":
                out.append('<span class="c-pre">%s</span>' % esc(tok))
            elif kind == "string":
                out.append('<span class="c-str">%s</span>' % esc(tok))
            elif kind == "number":
                out.append('<span class="c-num">%s</span>' % esc(tok))
            elif kind == "ident":
                out.append(self.ident(tok))
            else:
                out.append(esc(tok))
        return "".join(out)

    def ident(self, tok, cls_fn="c-fn"):
        sym = self.site.symbols.get(tok)
        if sym is not None:
            cls = "c-ty" if sym.kind == "type" else cls_fn
            return '<a class="%s" href="%s">%s</a>' % (cls, self.rel(sym.url), esc(tok))
        if tok in C_KEYWORDS:
            return '<span class="c-kw">%s</span>' % tok
        if tok in C_TYPES or re.match(r"^\w+_(t|ptr|srcptr|struct)$", tok):
            return '<span class="c-ty">%s</span>' % tok
        if re.match(r"^[A-Z][A-Z0-9_]{2,}$", tok):
            return '<span class="c-mac">%s</span>' % tok
        return esc(tok)

    def proto(self, p, name_href=None, anchor=False, src=True):
        """Render one prototype with the function name emphasised."""
        out = []
        if p.ret:
            out.append(self.type_html(p.ret) + " ")
        name = esc(p.name)
        if name_href:
            name = '<a href="%s">%s</a>' % (self.rel(name_href), name)
        out.append('<span class="p-name">%s</span>' % name)
        if p.params is not None:
            params = []
            for prm in p.params:
                params.append('<span class="p-arg">%s</span>' % self.param_html(prm))
            out.append('(<span class="p-args">%s</span>)' % ", ".join(params))
        idattr = ' id="%s"' % esc(p.name) if anchor else ""
        srcurl = self.site.source_url(p.name) if src else None
        srclink = ('<a class="p-src" href="%s" title="View the definition of %s on GitHub">source</a>'
                   % (esc(srcurl), esc(p.name))) if srcurl else ""
        return '<div class="proto"%s><code>%s</code>%s</div>' % (idattr, "".join(out), srclink)

    def type_html(self, t):
        return re.sub(r"[A-Za-z_]\w*", lambda m: self.ident(m.group(0)), esc(t))

    def param_html(self, prm):
        m = re.match(r"^(.*?)([A-Za-z_]\w*)(\s*\[[^\]]*\])?$", prm)
        if m and m.group(1).strip() and m.group(2) not in C_TYPES:
            return '%s<span class="p-var">%s</span>%s' % (self.type_html(m.group(1)), esc(m.group(2)),
                                                          esc(m.group(3) or ""))
        return self.type_html(prm)

    def div(self, nd):
        kind = nd["kind"]
        did = ' id="%s"' % esc(nd["id"]) if nd["id"] else ""
        if kind in THEOREM_KINDS:
            head = '<span class="env-head">%s %d</span>' % (kind.capitalize(), nd["number"])
            if nd["args"]:
                head += ' <span class="env-title">(%s)</span>' % self.inline(_parse(nd["args"]))
            body = self.blocks_with_lead(nd.children, head + ". ")
            return '<div class="env env-%s"%s>%s</div>' % (kind, did, body)
        if kind == "proof":
            title = "Proof" + (" of %s" % self.inline(_parse(nd["args"])) if nd["args"] else "")
            body = self.blocks_with_lead(nd.children, '<span class="proof-head">%s.</span> ' % title)
            return '<div class="proof"%s>%s<span class="qed" aria-label="end of proof"></span></div>' % (did, body)
        if kind == "algorithm":
            title = 'Algorithm %d' % nd["number"]
            if nd["args"]:
                title += ': %s' % self.inline(_parse(nd["args"]))
            return ('<div class="algo"%s><div class="algo-title">%s</div><div class="algo-body">%s</div></div>'
                    % (did, title, self.blocks(nd.children)))
        if kind in ("note", "warning"):
            label = nd["args"] or kind.capitalize()
            return ('<aside class="admon admon-%s"%s><div class="admon-label">%s</div>%s</aside>'
                    % (kind, did, esc(label), self.blocks(nd.children)))
        if kind == "summary":
            return self.summary_card(nd)
        if kind == "functions":
            items = []
            for name in nd["args"].split():
                sym = self.site.symbols.get(name)
                if sym is None:
                    self.site.warn("unknown function %s in ::: functions" % name, self.where(nd))
                    continue
                p = next(p for p in sym.block.protos if p.name == name)
                items.append(self.proto(p, sym.url, src=False))
            return '<div class="fnlist">%s</div>' % "".join(items)
        if kind == "modules":
            return self.module_listing()
        if kind == "details":
            return self.blocks(nd.children)
        if kind in ("function", "type", "macro", "variable"):
            return self.overview_block(nd)
        self.site.warn("unknown block ::: %s" % kind, self.where(nd))
        return self.blocks(nd.children)

    def module_listing(self):
        from .markup import meta_list
        site = self.site
        cats = meta_list(self.page.meta.get("module-categories")) if isinstance(
            self.page.meta.get("module-categories"), list) else []
        by_cat = {}
        for name, mp in site.modules_in_order():
            by_cat.setdefault(mp.meta.get("category", "Other"), []).append(mp)
        order = cats + sorted(c for c in by_cat if c not in cats)
        out = []
        for c in order:
            if c not in by_cat:
                continue
            items = "".join(
                '<a class="modcard k-module" href="%s"><code>%s</code><span>%s</span></a>'
                % (self.rel(mp.url), esc(mp.meta.get("header", mp.module + ".h")),
                   self.plain_inline(_parse(mp.meta.get("description", ""))))
                for mp in by_cat[c])
            out.append('<section class="modcat"><h3>%s</h3><div class="modcards">%s</div></section>'
                       % (esc(c), items))
        return "".join(out)

    def blocks_with_lead(self, nodes, lead):
        if nodes and nodes[0].kind == "para":
            first = "<p>%s%s</p>" % (lead, self.inline(nodes[0].children))
            return first + "\n" + self.blocks(nodes[1:])
        return "<p>%s</p>" % lead + self.blocks(nodes)

    def summary_card(self, nd):
        target = nd["args"].split()[0] if nd["args"] else ""
        tp = self.site.resolve_page(target)
        if tp is None:
            self.site.warn("unknown page in ::: summary %s" % target, self.where(nd))
            return ""
        return self.page_card(tp)

    def page_card(self, tp):
        """A card linking to a page: its kind, title and abstract. The kind
        also selects the card's background tone (class k-<kind>)."""
        summ = tp.summary()
        tone, kind = card_kind(tp)
        return ('<a class="card k-%s" href="%s"><span class="card-kind">%s</span><span class="card-title">%s</span>'
                '<span class="card-text">%s</span></a>'
                % (tone, self.rel(tp.url), kind, self.plain_inline(tp.title), self.plain_inline(summ) if summ else ""))

    def overview_block(self, nd):
        b = nd["block"]
        sigs = "".join(self.proto(p, self.site.symbols[p.name].url, anchor=True, src=False)
                       for p in b.protos)
        more = ""
        if b.details or b.page.kind == "detail":
            more = ('<a class="more" href="%s">%s</a>'
                    % (self.rel(anchor_url(b.page.url, b.names[0] if len(b.page.blocks) > 1 else "")),
                       "Full documentation" if b.page.kind == "detail" else "More"))
        return ('<div class="fblock fblock-%s">%s<div class="fdesc">%s%s</div></div>'
                % (b.kind, sigs, self.blocks(b.brief), more))

    def references(self):
        if not self.cites:
            return ""
        items = "".join('<li id="ref-%s"><span class="ref-key">[%s]</span> %s</li>'
                        % (esc(k), esc(k), self.inline(self.site.bib[k])) for k in self.cites)
        return ('<section class="refs"><h2 id="references">References<a class="hl" href="#references">#</a></h2>'
                '<ul class="reflist">%s</ul></section>' % items)


def _parse(s):
    from .markup import parse_inline
    return parse_inline(s)


def html_unescape(s):
    import html as _html
    return _html.unescape(s)


def wbr(s):
    """Allow long identifiers to wrap after underscores."""
    return re.sub(r"(?<=[A-Za-z0-9])_", "_<wbr>", s)


def slug(s):
    s = re.sub(r"[^\w\s-]", "", s.lower())
    return re.sub(r"[\s_]+", "-", s).strip("-")


# ---------------------------------------------------------------------------
# Page assembly
# ---------------------------------------------------------------------------

class HtmlBuilder:
    def __init__(self, site, outdir):
        self.site = site
        self.outdir = outdir
        os.makedirs(outdir, exist_ok=True)
        self.katex = Katex(os.path.join(outdir, ".katex-cache.json"), site.macros, site.warn)

    def build(self):
        # remove pages left over from earlier builds
        for dirpath, _, files in os.walk(self.outdir):
            for fn in files:
                if fn.endswith(".html"):
                    os.remove(os.path.join(dirpath, fn))
        from .progress import Progress
        from .parallel import pmap
        ids = [page.id for page in self.site.order]
        # render the pages in parallel; each worker numbers the formulas of
        # its pages itself, and the numbers are made global here
        pr = Progress("rendering pages", len(ids))
        results = pmap(_render_page_worker, ids, pr, state=self)
        rendered = {}
        for pid, (text, keys) in zip(ids, results):
            index = [self.katex.register(k, w) for k, w in keys]
            rendered[pid] = re.sub("\x00K(\\d+)\x00", lambda m: "\x00K%d\x00" % index[int(m.group(1))], text)
        pr.done()
        pr = Progress("typesetting math")
        self.katex.render_all()
        pr.done("KaTeX with Node.js" if self.katex.prerendered else "in the browser: Node.js not found")
        self._rendered = rendered
        pr = Progress("writing pages", len(ids))
        written = pmap(_write_page_worker, ids, pr, state=self)
        self._rendered = None
        bundles = {}
        for pid, text in zip(ids, written):
            if text is not None:
                page = self.site.pages[pid]
                bundles.setdefault(page.url.split("#")[0].rsplit("/", 1)[0], []).append((page, text))
        for module, pages in bundles.items():
            self.write_bundle(module, pages)
        self.copystatic()
        self.write_search_index()
        pr.done()
        return len(rendered)

    def write_bundle(self, module, pages):
        """Write <module>/f.html, a shell page, and <module>/f.js with the
        title, sidebar and main content of every bundled page. Sidebars are
        stored once per distinct content; the current entry is marked by the
        shell's script."""
        sidebars, index, data = [], {}, {}
        shell = None
        for page, text in pages:
            m_side = re.search(r'<aside class="sidebar" id="sidebar">(.*?)</aside>', text, re.S)
            m_main = re.search(r'<main class="main" id="main">\n(.*?)\n  </main>', text, re.S)
            m_title = re.search(r"<title>(.*?)</title>", text, re.S)
            side = m_side.group(1).replace(' aria-current="page"', "")
            if side not in index:
                index[side] = len(sidebars)
                sidebars.append(side)
            name = page.url.split("#", 1)[1]
            data[name] = [html_unescape(m_title.group(1)), index[side], m_main.group(1)]
            if shell is None:
                shell = (text[:m_side.start(1)] + text[m_side.end(1):m_main.start(1)]
                         + '<p class="stub">Loading…</p>' + text[m_main.end(1):])
        js = "window.FDOC_BUNDLE=" + json.dumps({"sidebars": sidebars, "pages": data},
                                                 ensure_ascii=False, separators=(",", ":")) + ";"
        d = os.path.join(self.outdir, module)
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "f.js"), "w", encoding="utf-8") as f:
            f.write(js)
        shell = re.sub(r'<script src="((?:\.\./)*)static/fdoc.js"></script>',
                       lambda m: ('<script src="f.js"></script><script src="%sstatic/bundle.js"></script>'
                                  '<script src="%sstatic/fdoc.js"></script>' % (m.group(1), m.group(1))),
                       shell)
        with open(os.path.join(d, "f.html"), "w", encoding="utf-8") as f:
            f.write(shell)

    def katex_js(self, text):
        def repl(m):
            if self.katex.prerendered:
                return ""
            return ('<script src="%sstatic/katex/katex.min.js"></script><script>'
                    'document.querySelectorAll(".tex").forEach(e=>katex.render(e.textContent,e,'
                    '{displayMode:e.classList.contains("display"),throwOnError:false}));</script>' % m.group(1))
        return re.sub("\x00KATEXJS([./]*)\x00", repl, text)

    def copystatic(self):
        st = os.path.join(self.outdir, "static")
        os.makedirs(os.path.join(st, "katex", "fonts"), exist_ok=True)
        for fn in ("fdoc.css", "fdoc.js", "bundle.js"):
            shutil.copy(os.path.join(HERE, "static", fn), os.path.join(st, fn))
        shutil.copy(os.path.join(KATEX_DIR, "katex.min.css"), os.path.join(st, "katex", "katex.min.css"))
        if not self.katex.prerendered:
            shutil.copy(os.path.join(KATEX_DIR, "katex.min.js"), os.path.join(st, "katex", "katex.min.js"))
        for fn in os.listdir(os.path.join(KATEX_DIR, "fonts")):
            shutil.copy(os.path.join(KATEX_DIR, "fonts", fn), os.path.join(st, "katex", "fonts", fn))

    def write_search_index(self):
        entries = []
        for name, sym in sorted(self.site.symbols.items()):
            brief = first_sentence(plain_text(sym.block.brief))
            entries.append([name, sym.url, sym.block.module, brief])
        for page in self.site.order:
            if page.kind in ("algorithm", "topic", "group", "tutorial", "module", "detail", "guide", "manual"):
                label = {"algorithm": "algorithm", "topic": "overview", "group": "overview",
                         "tutorial": "tutorial", "module": "module", "detail": page.module, "manual": "user guide",
                         "guide": "%s guide" % page.module}[page.kind]
                if page.kind == "manual" and page.collection != "guide":
                    label = COLLECTIONS.get(page.collection, label).lower()
                entries.append([page.title_text, page.url, label,
                                plain_text(page.summary() or []) if page.kind != "module" else
                                page.meta.get("description", "")])
        with open(os.path.join(self.outdir, "static", "searchindex.js"), "w", encoding="utf-8") as f:
            f.write("window.FDOC_INDEX=" + json.dumps(entries, ensure_ascii=False) + ";")

    # -- shared chrome ------------------------------------------------------------

    def render_page(self, page, r):
        kind = page.kind
        if kind == "module":
            main = self.module_page(page, r)
        elif kind in ("detail", "auto"):
            main = self.function_page(page, r)
        elif kind == "guide":
            main = self.guide_page(page, r)
        elif kind in PROSE_KINDS:
            main = self.prose_page(page, r)
        elif kind == "refs":
            main = self.refs_page(page, r)
        else:
            main = self.home_page(page, r)
        return self.layout(page, r, main, self.sidebar(page, r))

    def layout(self, page, r, main, sidebar):
        site = self.site
        cfg = site.config
        up = "../" * r.depth
        title = page.title_text
        full_title = "%s — %s" % (title, cfg["project"]) if page.kind != "home" else cfg["project"] + " documentation"
        src_rel = os.path.relpath(page.path, site.flintroot) if page.path else None
        edit = ""
        if src_rel:
            edit = ('<a href="%s/blob/main/%s">Edit this page</a> · ' % (cfg["github"], src_rel.replace(os.sep, "/")))
        nav = []
        for label, url, test in (("User guide", "guide/index.html", lambda p: p.collection == "guide"),
                                 ("Modules", "index.html", lambda p: p.kind in ("home", "module", "detail", "auto", "guide")),
                                 ("Algorithms", "algorithms/index.html", lambda p: p.collection == "algorithms"),
                                 ("Examples &amp; tutorials", "tutorials/index.html", lambda p: p.collection == "tutorials"),
                                 ("History", "history/index.html", lambda p: p.collection == "history"),
                                 ("References", "references.html", lambda p: p.kind == "refs")):
            cur = ' aria-current="page"' if test(page) else ""
            nav.append('<a href="%s%s"%s>%s</a>' % (up, url, cur, label))
        katex_js = "\x00KATEXJS%s\x00" % up
        cite = self.citation(page)
        # a site-wide notice ('notice:' in the front matter of the home
        # page), shown in a box at the top of every page
        home = site.pages.get("index")
        notice = home.meta.get("notice") if home else None
        notice = ('<aside class="notice" role="note">%s</aside>\n' % r.inline(_parse(notice))
                  if notice else "")
        # the PDF is published next to the HTML (flint.pdf in the base URL)
        pdf = cfg.get("baseurl", "") + "flint.pdf"
        return TEMPLATE.format(cite=cite, pdf=esc(pdf), notice=notice,
            title=esc(full_title), up=up, nav="".join(nav), sidebar=sidebar, main=main,
            project=esc(cfg["project"]), version=esc(cfg["version"]), edit=edit,
            commit=esc(site.commit[:10]), katex_js=katex_js,
            description=esc(first_sentence(plain_text(page.summary() or [])) or title))

    def citation(self, page):
        """'Cite this page' box: a ready-to-paste reference and a BibTeX entry."""
        cfg = self.site.config
        url = getattr(page, "canonical_url", page.url)
        if url == "index.html":
            url = ""
        full_url = cfg.get("baseurl", "") + url
        date = cfg.get("date", "")
        year = date[:4]
        version = cfg.get("version", "")
        title = page.title_text if page.kind != "home" else "FLINT documentation"
        mpage = self.site.modules.get(page.module) if page.module else None
        header = mpage.meta.get("header", page.module + ".h") if mpage else None
        if page.kind == "module":
            title = "%s (module overview)" % title
        elif header and page.kind in ("detail", "auto", "guide") and header not in title:
            title = "%s (%s)" % (title, header)
        series = "FLINT: Fast Library for Number Theory, documentation"
        ver = (", version %s" % version) if version else ""
        text = ('The FLINT team, \u201c%s,\u201d %s%s, %s. %s (generated %s).'
                % (title, series, ver, year, full_url, date))
        slug = re.sub(r"[^A-Za-z0-9]+", "-", url[:-5] if url.endswith(".html") else url).strip("-") or "index"
        bib = ("@misc{flint-doc-%s,\n"
               "  author       = {{The FLINT team}},\n"
               "  title        = {%s},\n"
               "  howpublished = {%s%s},\n"
               "  year         = {%s},\n"
               "  url          = {%s},\n"
               "  note         = {Generated %s}\n"
               "}" % (slug, bibtex_title(title), series, ver, year, full_url, date))
        return ('<section class="cite-page" aria-label="How to cite this page">'
                '<div class="cite-head"><span class="cite-label">Cite this page</span>'
                '<span class="cite-tabs" role="tablist">'
                '<button type="button" class="cite-tab" role="tab" aria-selected="true" data-cite="text">Text</button>'
                '<button type="button" class="cite-tab" role="tab" aria-selected="false" data-cite="bibtex">BibTeX</button>'
                '</span><button type="button" class="cite-copy">Copy</button></div>'
                '<p class="cite-body" data-cite-panel="text">%s</p>'
                '<pre class="cite-body" data-cite-panel="bibtex" hidden>%s</pre></section>'
                % (esc(text), esc(bib)))

    def sidebar(self, page, r):
        site = self.site
        out = []
        mod = page.module
        if page.kind in ("module", "detail", "auto", "guide") and mod in site.modules:
            mpage = site.modules[mod]
            out.append('<div class="sb-head"><a href="%s">%s</a></div>'
                       % (r.rel(mpage.url), esc(mpage.meta.get("header", mod + ".h"))))
            cur_section = page.blocks[0].section[0] if page.blocks else None
            guides = site.module_guides(mod)
            if guides:
                out.append('<div class="sb-label">Guides</div><ul class="sb-list sb-guides">')
                for g in guides:
                    cur = ' aria-current="page"' if g is page else ""
                    out.append('<li><a href="%s"%s>%s</a></li>' % (r.rel(g.url), cur, r.inline(g.title)))
                out.append('</ul><div class="sb-label">Reference</div>')
            out.append('<ul class="sb-list">')
            pages = site.module_pages(mod)
            for node in mpage.body:
                if node.kind == "heading" and node["level"] == 2:
                    sid = node["id"]
                    open_ = sid == cur_section
                    out.append('<li class="sb-sec%s"><a href="%s">%s</a>'
                               % (" open" if open_ else "", r.rel(anchor_url(mpage.url, sid)), r.inline(node.children)))
                    if open_:
                        out.append('<ul>')
                        for p in pages:
                            if p.blocks and p.blocks[0].section[0] == sid:
                                cur = ' aria-current="page"' if p is page else ""
                                label = (r.inline(p.title) if p.kind == "detail"
                                         else "<code>%s</code>" % wbr(esc(p.blocks[0].names[0])))
                                out.append('<li><a href="%s"%s>%s</a></li>' % (r.rel(p.url), cur, label))
                        out.append('</ul>')
                    out.append('</li>')
            out.append('</ul>')
        elif page.collection:
            out.append(self.collection_sidebar(page, r, page.collection))
        else:
            out.append('<div class="sb-head"><a href="%s">Modules</a></div><ul class="sb-list">' % r.rel("index.html"))
            for name, mp in sorted(site.modules.items()):
                out.append('<li class="sb-plain"><a href="%s"><code>%s</code></a></li>'
                           % (r.rel(mp.url), esc(mp.meta.get("header", name + ".h"))))
            out.append('</ul>')
            for coll in COLLECTIONS:
                if self.site.collection_root(coll):
                    out.append(self.collection_sidebar(page, r, coll))
        return "".join(out)

    def collection_sidebar(self, page, r, coll):
        site = self.site
        root = site.collection_root(coll)
        out = ['<div class="sb-head"><a href="%s">%s</a></div><ul class="sb-list">'
               % (r.rel(root.url), COLLECTIONS[coll])]
        for g, pages in site.collection_groups(coll):
            items = []
            for p in pages:
                cur = ' aria-current="page"' if p is page else ""
                items.append('<li><a href="%s"%s>%s</a></li>' % (r.rel(p.url), cur, r.inline(p.title)))
            if g is None:
                out.append("".join(items).replace("<li>", '<li class="sb-plain">'))
            else:
                cur = ' aria-current="page"' if g is page else ""
                label = esc(g.meta["short"]) if g.meta.get("short") else r.inline(g.title)
                sub = "<ul>%s</ul>" % "".join(items) if items else ""
                out.append('<li class="sb-sec open"><a href="%s"%s>%s</a>%s</li>'
                           % (r.rel(g.url), cur, label, sub))
        out.append("</ul>")
        return "".join(out)

    def card(self, r, tp):
        return r.page_card(tp)

    def crumbs(self, r, items):
        parts = ['<a href="%s">%s</a>' % (r.rel(u), t) if u else '<span>%s</span>' % t for t, u in items]
        return '<nav class="crumbs" aria-label="Breadcrumb">%s</nav>' % '<span class="sep">/</span>'.join(parts)

    # -- page kinds ------------------------------------------------------------------

    def module_page(self, page, r):
        hdr = page.meta.get("header", page.module + ".h")
        sections = [nd for nd in page.body if nd.kind == "heading" and nd["level"] == 2]
        toc = "".join('<li><a href="#%s">%s</a></li>' % (nd["id"], r.inline(nd.children)) for nd in sections)
        first = next((i for i, nd in enumerate(page.body) if nd.kind == "heading" and nd["level"] == 2),
                     len(page.body))
        intro = r.blocks(page.body[:first])
        body = r.blocks(page.body[first:])
        stats = '%d functions' % sum(len(b.protos) for b in self.site.blocks if b.module == page.module)
        guides = self.site.module_guides(page.module)
        if guides:
            intro += ('<section class="group"><h2 id="guides">Guides<a class="hl" href="#guides">#</a></h2>'
                      '<div class="cards">%s</div></section>' % "".join(self.card(r, g) for g in guides))
        return ('%s<header class="page-head"><div class="eyebrow">Module · %s</div>'
                '<h1><code>%s</code></h1>'
                '<pre class="include"><code><span class="c-pre">#include &lt;flint/%s&gt;</span></code></pre></header>'
                '<div class="prose intro">%s</div>'
                '<nav class="toc" aria-label="Sections"><div class="toc-label">Sections</div><ol>%s</ol></nav>'
                '<div class="prose overview">%s</div>%s'
                % (self.crumbs(r, [("FLINT", "index.html"), (esc(page.module), None)]), stats,
                   esc(hdr), esc(hdr), intro, toc, body,
                   r.references()))

    def function_page(self, page, r):
        site = self.site
        mpage = site.modules[page.module]
        b0 = page.blocks[0]
        trail = [("FLINT", "index.html"),
                 ("<code>%s</code>" % esc(page.module), mpage.url),
                 (r.inline(b0.section[1]), anchor_url(mpage.url, b0.section[0] or ""))]
        subs = {b.subsection[0] for b in page.blocks if getattr(b, "subsection", None)}
        if len(subs) == 1 and all(getattr(b, "subsection", None) for b in page.blocks):
            trail.append((r.inline(b0.subsection[1]), anchor_url(mpage.url, b0.subsection[0])))
        crumbs = self.crumbs(r, trail)
        synopsis = []
        # With several blocks on one page, each block's extended notes stay
        # with its prototypes so that it is clear which function they describe.
        inline_details = len(page.blocks) > 1
        for b in page.blocks:
            sigs = "".join(r.proto(p, None, anchor=True) for p in b.protos)
            more = ('<div class="syn-details">%s</div>' % r.blocks(b.details)
                    if inline_details and b.details else "")
            synopsis.append('<div class="syn-block">%s<div class="syn-brief">%s</div>%s</div>'
                            % (sigs, r.blocks(b.brief), more))
        if page.kind == "detail":
            title = '<h1>%s</h1>' % r.inline(page.title)
        else:
            title = '<h1 class="fn-title"><code>%s</code></h1>' % ", ".join(esc(n) for n in b0.names)
        details = "" if inline_details else "".join(r.blocks(b.details) for b in page.blocks)
        body = r.blocks(page.body)
        extended = details + body
        # neighbours in module order
        pages = site.module_pages(page.module)
        i = pages.index(page)
        prev_ = pages[i - 1] if i > 0 else None
        next_ = pages[i + 1] if i + 1 < len(pages) else None
        def label(p):
            text = r.inline(p.title) if p.kind == "detail" else "<code>%s</code>" % esc(p.blocks[0].names[0])
            if p.blocks[0].section[0] != b0.section[0]:
                text = '<small class="pn-sec">%s</small>%s' % (r.plain_inline(p.blocks[0].section[1]), text)
            return text
        guides = [g for g in site.order if g.kind == "guide"
                  and any(n in site.guide_mentions.get(g.id, ()) for b in page.blocks for n in b.names)]
        guide_links = ""
        if guides:
            guide_links = ('<section class="usedby"><h2 id="guides">Guides<a class="hl" href="#guides">#</a></h2>'
                           '<ul>%s</ul></section>' % "".join(
                               '<li><a href="%s">%s</a>%s</li>' % (r.rel(g.url), r.inline(g.title),
                                   "" if g.module == page.module else
                                   ' <span class="muted">(<code>%s</code>)</span>' % esc(g.module))
                               for g in guides))
        pn = '<nav class="pn">%s%s</nav>' % (
            '<a class="prev" href="%s"><span>Previous</span>%s</a>' % (r.rel(prev_.url), label(prev_)) if prev_ else "<span></span>",
            '<a class="next" href="%s"><span>Next</span>%s</a>' % (r.rel(next_.url), label(next_)) if next_ else "<span></span>")
        hdr = mpage.meta.get("header", page.module + ".h")
        return ('%s<header class="page-head">%s</header>'
                '<section class="synopsis" aria-label="Synopsis"><div class="syn-head"><span>Synopsis</span>'
                '<code class="c-pre">#include &lt;flint/%s&gt;</code></div>%s</section>'
                '<div class="prose">%s</div>%s%s%s'
                % (crumbs, title, esc(hdr), "".join(synopsis), extended, guide_links, r.references(), pn))

    def guide_page(self, page, r):
        site = self.site
        mpage = site.modules[page.module]
        crumbs = self.crumbs(r, [("FLINT", "index.html"),
                                 ("<code>%s</code>" % esc(page.module), mpage.url),
                                 ("Guides", anchor_url(mpage.url, "guides"))])
        hdr = mpage.meta.get("header", page.module + ".h")
        return ('%s<header class="page-head"><div class="eyebrow">Guide · <code>%s</code></div><h1>%s</h1></header>'
                '<div class="prose">%s</div>%s'
                % (crumbs, esc(hdr), r.inline(page.title), r.blocks(page.body), r.references()))

    def prose_page(self, page, r):
        site = self.site
        body = r.blocks(page.body)
        used = ""
        links = site.backlinks.get(page.id, [])
        if links and page.kind in ("algorithm", "topic", "group"):
            items = []
            for p, b in links:
                if b is not None:
                    names = ", ".join('<a class="sym" href="%s"><code>%s</code></a>'
                                      % (r.rel(site.symbols[n].url), esc(n)) for n in b.names)
                    items.append("<li>%s</li>" % names)
                elif p.kind in ("detail", "auto", "guide"):
                    items.append('<li><a href="%s">%s</a></li>' % (r.rel(p.url), r.inline(p.title)))
                elif p.kind in PROSE_KINDS:
                    kind = {"tutorial": "tutorial", "algorithm": "algorithm"}.get(p.kind, "overview")
                    items.append('<li><a href="%s">%s</a> <span class="muted">(%s)</span></li>'
                                 % (r.rel(p.url), r.inline(p.title), kind))
                elif p.kind == "module":
                    items.append('<li><a href="%s">%s</a> <span class="muted">(module)</span></li>'
                                 % (r.rel(p.url), r.inline(p.title)))
            used = ('<section class="usedby"><h2 id="used-by">Referenced from<a class="hl" href="#used-by">#</a></h2>'
                    '<ul>%s</ul></section>' % "".join(items))
        coll = page.collection
        root = site.collection_root(coll)
        trail = [("FLINT", "index.html")]
        if page.kind == "collection":
            trail.append((COLLECTIONS[coll], None))
        else:
            trail.append((COLLECTIONS[coll], root.url))
            if page.group and page.kind != "group":
                g = site.pages.get("%s/%s/index" % (coll, page.group))
                if g:
                    trail.append((esc(g.meta.get("short") or g.title_text), g.url))
        crumbs = self.crumbs(r, trail)
        eyebrow = {"algorithm": "Algorithm", "topic": "Overview", "group": "Overview",
                   "collection": "", "tutorial": "Tutorial", "manual": ""}[page.kind]
        meta = ""
        if page.meta.get("complexity"):
            meta = '<dl class="facts"><dt>Cost</dt><dd>%s</dd></dl>' % r.inline(_parse(page.meta["complexity"]))
        listing = ""
        if page.kind == "collection":
            parts = []
            for g, pages in site.collection_groups(coll):
                cards = "".join(self.card(r, p) for p in pages)
                if g is None:
                    parts.append('<div class="cards">%s</div>' % cards)
                else:
                    parts.append('<section class="group"><h2 id="%s"><a href="%s">%s</a></h2><p class="group-summary">%s</p>'
                                 '<div class="cards">%s</div></section>'
                                 % (esc(g.group), r.rel(g.url), r.inline(g.title),
                                    r.inline(g.summary()) if g.summary() else "", cards))
            listing = "".join(parts)
        elif page.kind == "group":
            pages = next((ps for g, ps in site.collection_groups(coll) if g is page), [])
            if pages:
                listing = ('<section class="group"><h2 id="pages">In this section<a class="hl" href="#pages">#</a></h2>'
                           '<div class="cards">%s</div></section>' % "".join(self.card(r, p) for p in pages))
        eb = '<div class="eyebrow">%s</div>' % eyebrow if eyebrow else ""
        return ('%s<header class="page-head">%s<h1>%s</h1>%s</header>'
                '<div class="prose">%s</div>%s%s%s'
                % (crumbs, eb, r.inline(page.title), meta, body, listing, used,
                   r.references()))

    def refs_page(self, page, r):
        items = "".join('<li id="ref-%s"><span class="ref-key">[%s]</span> %s</li>'
                        % (esc(k), esc(k), r.inline(v)) for k, v in sorted(self.site.bib.items(), key=lambda kv: kv[0].lower()))
        return ('<header class="page-head"><h1>%s</h1></header><ul class="reflist all">%s</ul>'
                % (r.inline(page.title), items))

    def home_page(self, page, r):
        return '<div class="prose home">%s</div>' % r.blocks(page.body)


def bibtex_title(s):
    """A page title as a BibTeX title: braces keep the capitalisation,
    LaTeX special characters are escaped."""
    out = []
    for ch in s:
        if ch in "&%$#_{}":
            out.append("\\" + ch)
        elif ch == "\\":
            out.append("\\textbackslash{}")
        elif ch in "^~":
            out.append("\\%s{}" % ("textasciicircum" if ch == "^" else "textasciitilde"))
        else:
            out.append(ch)
    return "{" + "".join(out) + "}"


class _MathCollector:
    """Stands in for Katex while a page is rendered in a worker process:
    records the formulas in the order of their placeholders."""

    def __init__(self, katex):
        self.katex = katex
        self.items = {}
        self.keys = []

    def placeholder(self, tex, display, where=""):
        tex = self.katex.expand(tex)
        key = ("D" if display else "I") + tex
        if key not in self.items:
            self.items[key] = len(self.keys)
            self.keys.append((key, where))
        return "\x00K%d\x00" % self.items[key]


def _render_page_worker(pid):
    from .parallel import state
    builder = state()
    page = builder.site.pages[pid]
    math = _MathCollector(builder.katex)
    text = builder.render_page(page, HtmlRenderer(builder.site, page, math))
    return text, math.keys


def _write_page_worker(pid):
    """Write one page; for a bundled page, return its text instead."""
    from .parallel import state
    builder = state()
    page = builder.site.pages[pid]
    text = builder.katex_js(builder.katex.substitute(builder._rendered[pid]))
    if "#" in page.url:
        return text
    path = os.path.join(builder.outdir, page.url)
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write(text)
    return None


def card_kind(page):
    """(tone, label) of a page's card: the label names the kind of page,
    the tone (a CSS class suffix) gives each kind its own shade of gray."""
    if page.kind == "manual":
        if page.collection == "guide":
            return "manual", "User guide"
        return page.collection or "manual", COLLECTIONS.get(page.collection, "")
    return {"algorithm": ("algorithm", "Algorithm"), "topic": ("overview", "Overview"),
            "group": ("overview", "Overview"), "tutorial": ("tutorial", "Tutorial"),
            "module": ("module", "Module"), "guide": ("guide", "Guide")}.get(page.kind, ("other", ""))


def first_sentence(s):
    s = " ".join(s.split())
    m = re.match(r"^(.+?[.!?])(\s|$)", s)
    out = m.group(1) if m else s
    return out if len(out) < 200 else out[:197] + "…"


TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">
<meta name="claude:viewport-fit-sides" content="cover">
<title>{title}</title>
<meta name="description" content="{description}">
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:ital,wght@0,400;0,500;0,600;1,400&family=IBM+Plex+Sans:wght@400;500;600&family=Source+Serif+4:ital,opsz,wght@0,8..60,400;0,8..60,600;1,8..60,400&display=swap">
<link rel="stylesheet" href="{up}static/katex/katex.min.css">
<link rel="stylesheet" href="{up}static/fdoc.css">
<script>try{{var t=localStorage.getItem("fdoc-theme");if(t)document.documentElement.dataset.theme=t}}catch(e){{}}</script>
</head>
<body>
<header class="topbar">
  <div class="topbar-in">
    <button class="sb-toggle" type="button" aria-label="Toggle navigation" aria-expanded="false" id="sb-toggle"><span></span></button>
    <a class="brand" href="{up}index.html"><span class="brand-name">{project}</span><span class="brand-ver">{version}</span></a>
    <nav class="topnav">{nav}</nav>
    <div class="search">
      <input id="search" type="search" placeholder="Search functions" autocomplete="off" spellcheck="false" aria-label="Search functions and pages">
      <kbd>/</kbd>
      <ul id="search-results" class="search-results" role="listbox" hidden></ul>
    </div>
    <button class="theme" type="button" id="theme-toggle" aria-label="Toggle dark mode"></button>
  </div>
</header>
<div class="shell">
  <aside class="sidebar" id="sidebar">{sidebar}</aside>
  <main class="main" id="main">
{notice}{main}
    {cite}
    <footer class="foot">{edit}<a href="{pdf}">PDF version</a> · Built with fdoc from <code>{commit}</code></footer>
  </main>
</div>
<script>window.FDOC_ROOT="{up}";</script>
<script src="{up}static/searchindex.js"></script>
<script src="{up}static/fdoc.js"></script>
{katex_js}
</body>
</html>
"""
