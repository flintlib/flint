"""
The documentation model: pages, function blocks, symbols, bibliography and
cross-references.  Everything here is backend independent.
"""

import os
import re
import subprocess

from .outputs import check_flags as check_output_flags
from .markup import (BlockParser, MarkupError, Node, meta_list, parse_inline,
                     plain_text, split_front_matter)
from .csyntax import Prototype, SourceIndex

THEOREM_KINDS = ("theorem", "lemma", "proposition", "corollary", "definition",
                 "conjecture", "remark")
COLLECTIONS = {"guide": "User guide", "algorithms": "Algorithms", "tutorials": "Examples & tutorials",
               "history": "History"}
PROSE_KINDS = ("algorithm", "topic", "group", "collection", "tutorial", "manual")

DIV_KINDS = set(THEOREM_KINDS) | {"proof", "algorithm", "note", "warning", "details",
                                  "function", "type", "macro", "variable", "summary", "functions", "modules",
                                  "figure", "steps"}


class Warn:
    def __init__(self):
        self.items = []

    def __call__(self, msg, path=None, line=None):
        loc = path or ""
        if line:
            loc += ":%d" % line
        self.items.append((loc + ": " if loc else "") + msg)


class Page:
    def __init__(self, pid, kind, title, meta=None, body=None, path=None, module=None):
        self.id = pid
        self.kind = kind            # home, module, detail, auto, algorithm, topic, refs
        self.title = title          # inline nodes
        self.meta = meta or {}
        self.body = body or []      # block nodes
        self.path = path            # source file (or None for generated pages)
        self.module = module
        self.blocks = []            # function blocks shown on this page (detail/auto)
        self.url = pid + ".html"
        self.numbers = {}           # label -> (kind, number)
        self.examples = []
        self.collection = None
        self.group = None
        self.headings = {}          # heading id -> title nodes

    @property
    def title_text(self):
        return plain_text(self.title)

    def summary(self):
        s = self.meta.get("summary") or self.meta.get("description")
        return parse_inline(s) if s else None


class Block:
    """A group of functions documented together (one '::: function' div)."""

    def __init__(self, kind, protos, brief, details, module, section, order, path, line):
        self.kind = kind
        self.protos = protos
        self.brief = brief
        self.details = details
        self.module = module
        self.section = section      # (id, title nodes)
        self.order = order
        self.path = path
        self.line = line
        self.page = None

    @property
    def names(self):
        return [p.name for p in self.protos]


def anchor_url(url, anchor):
    """url#anchor; bundled pages already use the fragment for the page
    (dir/f.html#page), so their anchors follow a colon."""
    if not anchor:
        return url
    return url + (":" if "#" in url else "#") + anchor


class Symbol:
    def __init__(self, name, kind, block):
        self.name, self.kind, self.block = name, kind, block

    @property
    def url(self):
        page = self.block.page
        if len(self.block.page.blocks) == 1 and len(self.block.protos) == 1:
            return page.url
        # bundled automatic pages already use the fragment ("f.html#page"),
        # so the symbol is appended after a colon
        return anchor_url(page.url, self.name)


class Site:
    def __init__(self, srcroot, flintroot, config):
        self.srcroot = srcroot
        self.flintroot = flintroot
        self.config = config
        self.pages = {}
        self.order = []
        self.symbols = {}
        self.blocks = []
        self.modules = {}           # module name -> module page
        self.bib = {}               # key -> inline nodes
        self.warn = Warn()
        self.sources = SourceIndex(os.path.join(flintroot, "src"))
        self.backlinks = {}         # page id -> list of (page, block or None)
        self.macros = {}
        self.commit = self._git_commit()

    # -- loading -----------------------------------------------------------

    def _git_commit(self):
        try:
            return subprocess.check_output(["git", "-C", self.flintroot, "rev-parse", "HEAD"],
                                           stderr=subprocess.DEVNULL).decode().strip()
        except Exception:
            # a release tarball: link to the release tag
            version = self.config.get("version", "")
            return "v" + version if version and "dev" not in version else "main"

    def load(self):
        from .progress import Progress
        self._load_macros()
        self._load_bib()
        paths = []
        for dirpath, dirnames, filenames in os.walk(self.srcroot):
            dirnames.sort()
            for fn in sorted(filenames):
                if fn.endswith(".md") and not (dirpath == self.srcroot and fn.startswith("references")):
                    paths.append(os.path.join(dirpath, fn))
        from .parallel import pmap
        pr = Progress("reading sources", len(paths))
        parsed = pmap(_parse_file, paths, pr)
        pr.done()
        self._scan_queue = []
        for path, res in zip(paths, parsed):
            self._load_page(path, res)
        # scan the C sources of the documented modules (in parallel), then
        # merge the results in the order of the pages
        mods = list(dict.fromkeys(self._scan_queue))
        pr = Progress("scanning C sources", len(mods))
        data = pmap(_collect_module, mods, pr, state=self.sources)
        for mod, d in zip(mods, data):
            self.sources.merge(mod, d)
        pr.done()
        # a failed preprocessor run loses the declarations generated by
        # macros for that configuration: say why, rather than only reporting
        # these functions as missing from the source
        pages = dict((p.meta.get("source", m), p) for m, p in self.modules.items())
        for mod, errs in self.sources.preprocess_errors.items():
            page = pages.get(mod)
            for flags, msg in errs:
                self.warn("preprocessing %s.h with '%s' failed (%s); the functions generated by macros "
                          "for these flags are not checked:\n    %s"
                          % (mod, flags, os.environ.get("CC", "cc"), msg.replace("\n", "\n    ")),
                          page.path if page else None)
        pr = Progress("cross-referencing")
        self._assign_blocks()
        self._collect_guide_mentions()
        self._number_and_collect()
        self._check_links()
        pr.done()
        pr = Progress("checking headers")
        self._check_sources()
        pr.done()

    def _load_macros(self):
        path = os.path.join(self.srcroot, "macros.tex")
        if not os.path.exists(path):
            return
        with open(path, encoding="utf-8") as f:
            for ln in f:
                m = re.match(r"\\newcommand\{(\\\w+)\}(?:\[(\d)\])?\{(.*)\}\s*$", ln.strip())
                if m:
                    self.macros[m.group(1)] = m.group(3)

    def _load_bib(self):
        """references.md plus any references-*.md files (one per contributor
        or topic), all in the format '[Key]: text'."""
        import glob
        path = os.path.join(self.srcroot, "references.md")
        paths = ([path] if os.path.exists(path) else []) + sorted(
            glob.glob(os.path.join(self.srcroot, "references-*.md")))
        meta = {}
        for p in paths:
            with open(p, encoding="utf-8") as f:
                text = f.read()
            m, body, off = split_front_matter(text)
            if p == path:
                meta = m
            for mm in re.finditer(r"^\[([\w.-]+)\]:\s*(.+?)\s*$", body, re.M):
                if mm.group(1) in self.bib:
                    self.warn("duplicate reference [%s]" % mm.group(1), p)
                self.bib[mm.group(1)] = parse_inline(mm.group(2))
        if paths:
            page = Page("references", "refs", parse_inline(meta.get("title", "References")),
                        meta, [], path)
            self._add_page(page)

    def _add_page(self, page):
        if page.id in self.pages:
            self.warn("duplicate page id %r" % page.id, page.path)
        self.pages[page.id] = page
        self.order.append(page)

    def _load_page(self, path, parsed=None):
        rel = os.path.relpath(path, self.srcroot)[:-3].replace(os.sep, "/")
        if parsed is None:
            parsed = _parse_file(path)
        if parsed[0] == "error":
            _, msg, epath, line = parsed
            self.warn(msg, epath or path, line)
            return
        _, meta, blocks = parsed
        parts = rel.split("/")
        title = meta.get("title")
        if title is None and blocks and blocks[0].kind == "heading" and blocks[0]["level"] == 1:
            title_nodes = blocks.pop(0).children
        else:
            title_nodes = parse_inline(title or parts[-1])
        collection = group = None
        if rel == "index":
            kind, module = "home", None
        elif parts[0] in COLLECTIONS:
            module = None
            collection = parts[0]
            if len(parts) == 2 and parts[1] == "index":
                kind = "collection"
            elif len(parts) == 3 and parts[2] == "index":
                kind, group = "group", parts[1]
            else:
                group = parts[1] if len(parts) == 3 else None
                if collection == "tutorials":
                    kind = "tutorial"
                elif collection in ("guide", "history"):
                    kind = "manual"
                else:
                    kind = "topic" if meta.get("kind") == "topic" else "algorithm"
        elif parts[-1] == "index":
            kind, module = "module", parts[0]
        elif meta.get("functions") or meta.get("kind") == "functions":
            kind, module = "detail", parts[0]
        else:
            # a hand-written page in a module directory that documents no
            # functions of its own: a guide, a direct child of the module
            kind, module = "guide", parts[0]
        page = Page(rel, kind, title_nodes, meta, blocks, path, module)
        page.collection, page.group = collection, group
        if kind == "module":
            page.id = module
            page.url = module + "/index.html"
            self.modules[module] = page
            if meta.get("template-headers"):
                self.sources.template_headers[meta.get("source", module)] = meta_list(meta["template-headers"])
            if meta.get("preprocess"):
                self.sources.preprocess_flags[meta.get("source", module)] = (
                    meta["preprocess"] if isinstance(meta["preprocess"], list) else [meta["preprocess"]])
            self._scan_queue.append(meta.get("source", module))
            for extra in meta_list(meta.get("extra-headers")):
                self._scan_queue.append(extra)
            self._extract_blocks(page)
        self._add_page(page)

    def _extract_blocks(self, page):
        """Pull '::: function' divs out of a module overview."""
        section = (None, [])
        subsection = None
        n = 0
        for node in page.body:
            if node.kind == "heading" and node["level"] == 2:
                sid = node["id"] or slugify(plain_text(node.children))
                node.attrs["id"] = sid
                section = (sid, node.children)
                subsection = None
            elif node.kind == "heading" and node["level"] == 3:
                sid = node["id"] or slugify(plain_text(node.children))
                node.attrs["id"] = sid
                subsection = (sid, node.children)
            if node.kind == "div" and node["kind"] in ("function", "type", "macro", "variable"):
                brief, details = [], []
                for c in node.children:
                    if c.kind == "div" and c["kind"] == "details":
                        details.extend(c.children)
                    else:
                        brief.append(c)
                try:
                    protos = [Prototype(p) for p in node["protos"]]
                except ValueError as e:
                    self.warn(str(e), page.path, node["line"])
                    continue
                b = Block(node["kind"], protos, brief, details, page.module, section, n,
                          page.path, node["line"])
                b.subsection = subsection
                n += 1
                node.attrs["block"] = b
                self.blocks.append(b)
                for p in protos:
                    if p.name in self.symbols:
                        self.warn("%s documented twice" % p.name, page.path, node["line"])
                    self.symbols[p.name] = Symbol(p.name, node["kind"], b)

    def _assign_blocks(self):
        # detail pages claim functions; everything else gets an automatic page
        for page in list(self.order):
            if page.kind != "detail":
                continue
            for name in meta_list(page.meta.get("functions")):
                sym = self.symbols.get(name)
                if sym is None:
                    self.warn("page claims unknown function %s" % name, page.path)
                    continue
                b = sym.block
                if b.page is not None and b.page is not page:
                    self.warn("%s claimed by both %s and %s" % (name, b.page.id, page.id),
                              page.path)
                    continue
                if b.page is None:
                    b.page = page
                    page.blocks.append(b)
            page.blocks.sort(key=lambda b: b.order)
            secs = []
            for b in page.blocks:
                if b.section[0] not in secs:
                    secs.append(b.section[0])
            if len(secs) > 1:
                self.warn("function page %s claims functions from several sections (%s); "
                          "a function page belongs to one section, use a guide for broader topics"
                          % (page.id, ", ".join(str(x) for x in secs)), page.path)
        for b in self.blocks:
            if b.page is None:
                pid = "%s/%s" % (b.module, b.names[0])
                title = [Node("code", [], text=", ".join(b.names))]
                page = Page(pid, "auto", title, {}, [], b.path, b.module)
                page.blocks = [b]
                b.page = page
                self._add_page(page)
        for page in self.order:
            if page.kind == "detail" and not page.blocks:
                self.warn("detail page documents no functions", page.path)

    def _collect_guide_mentions(self):
        """Symbols that each module guide mentions (in text or in its programs),
        so that function pages can link back to the guides that use them."""
        self.guide_mentions = {}
        for page in self.order:
            if page.kind != "guide":
                continue
            names = set()
            for nd in Node("root", page.body).walk():
                text = None
                if nd.kind == "code":
                    text = nd.get("text") or ""
                elif nd.kind == "div" and nd.get("kind") == "functions":
                    text = nd.get("args") or ""
                if text:
                    names.update(w for w in re.findall(r"[A-Za-z_]\w*", text) if w in self.symbols)
            self.guide_mentions[page.id] = names

    # -- collections (algorithms, tutorials) ----------------------------------

    def _sort_key(self, p):
        try:
            order = float(p.meta.get("order", 50))
        except ValueError:
            order = 50
        return (order, p.title_text.lower())

    def collection_root(self, coll):
        return self.pages.get(coll + "/index")

    def collection_groups(self, coll):
        """[(group page or None, [pages])] in display order."""
        groups = sorted((p for p in self.order if p.collection == coll and p.kind == "group"),
                        key=self._sort_key)
        out = []
        loose = sorted((p for p in self.order if p.collection == coll and p.group is None
                        and p.kind not in ("collection", "group")), key=self._sort_key)
        if loose:
            out.append((None, loose))
        for g in groups:
            pages = sorted((p for p in self.order if p.collection == coll and p.group == g.group
                            and p.kind != "group"), key=self._sort_key)
            out.append((g, pages))
        return out

    # -- cross references ---------------------------------------------------

    def bundle_auto_pages(self):
        """Serve all pages of each directory (module, or group of a
        collection) from one shell page (<dir>/f.html#<name>) instead of one
        file each, for hosts that limit the number of files. The home page,
        the references and the collection roots stay separate files."""
        for page in self.order:
            if page.kind in ("home", "refs", "collection"):
                continue
            d, _, name = page.url[:-5].rpartition("/") if page.url.endswith(".html") else ("", "", "")
            if not d:
                continue
            page.canonical_url = page.url
            page.url = "%s/f.html#%s" % (d, name)
            page._anchors = None

    def modules_in_order(self):
        """Modules grouped by the home page's module-categories, then by name."""
        home = self.pages.get("index")
        cats = home.meta.get("module-categories") if home else None
        cats = cats if isinstance(cats, list) else []
        def key(item):
            name, page = item
            c = page.meta.get("category")
            try:
                order = float(page.meta.get("order", 50))
            except ValueError:
                order = 50
            return (cats.index(c) if c in cats else len(cats), order, name)
        return sorted(self.modules.items(), key=key)

    def module_guides(self, module):
        """Guides of a module, in display order."""
        return sorted((p for p in self.order if p.kind == "guide" and p.module == module),
                      key=self._sort_key)

    def module_pages(self, module):
        """Detail/auto pages of a module in overview order."""
        seen, out = set(), []
        for b in sorted((b for b in self.blocks if b.module == module), key=lambda b: b.order):
            if b.page.id not in seen:
                seen.add(b.page.id)
                out.append(b.page)
        return out

    def page_content(self, page):
        """All block-level content rendered on a page (for numbering and link collection)."""
        content = []
        for b in page.blocks:
            content.extend(b.brief)
            content.extend(b.details)
        content.extend(page.body)
        return content

    def _number_and_collect(self):
        for page in self.order:
            counters = {"thm": 0, "alg": 0, "eq": 0}
            content = self.page_content(page)
            root = Node("root", content)
            ex = 0
            for node in root.walk():
                if node.kind == "heading":
                    hid = node["id"] or slugify(plain_text(node.children))
                    page.headings.setdefault(hid, node.children)
                if node.kind == "div" and node["kind"] in THEOREM_KINDS + ("algorithm",):
                    key = "alg" if node["kind"] == "algorithm" else "thm"
                    counters[key] += 1
                    node.attrs["number"] = counters[key]
                    if node["id"]:
                        page.numbers[node["id"]] = (node["kind"].capitalize(), counters[key])
                if node.kind == "code" and node.get("flags", {}).get("example"):
                    ex += 1
                    name = node["flags"].get("name") or "%s_%d" % (page.id.split("/")[-1], ex)
                    node.attrs["example_name"] = name
                    page.examples.append(node)
                    for problem in check_output_flags(node.get("output_flags") or {}):
                        self.warn("example %s: %s" % (name, problem), page.path, node.get("line"))
            # backlinks from function documentation to other pages (algorithms)
            if page.kind == "module":
                for node in page.body:
                    if node.kind == "div" and node.get("block"):
                        self._collect_links(node, page, node["block"])
                    else:
                        self._collect_links(node, page, None)
            else:
                for b in page.blocks:
                    for nd in b.brief + b.details:
                        self._collect_links(nd, page, b)
                for nd in page.body:
                    self._collect_links(nd, page, None)

    def _collect_links(self, node, page, block):
        for nd in node.walk():
            target = None
            if nd.kind == "link":
                target = nd["target"]
            elif nd.kind == "div" and nd["kind"] == "summary":
                target = nd["args"].split()[0] if nd["args"] else None
            if not target:
                continue
            tp = self.resolve_page(target.split("#")[0])
            if tp is not None and tp is not page and tp.kind in PROSE_KINDS:
                lst = self.backlinks.setdefault(tp.id, [])
                if block is not None:
                    entry = (block.page, block)
                else:
                    entry = (page, None)
                if entry not in lst:
                    lst.append(entry)

    def resolve_page(self, target):
        target = target.rstrip("/")
        if target in self.pages:
            return self.pages[target]
        if target + "/index" in self.pages:
            return self.pages[target + "/index"]
        if target in self.modules:
            return self.modules[target]
        if target.endswith("/index") and target[:-6] in self.modules:
            return self.modules[target[:-6]]
        return None

    def anchors(self, page):
        """Every fragment id a link may point to on a page."""
        if getattr(page, "_anchors", None) is None:
            ids = set(page.headings) | set(page.numbers)
            for b in page.blocks:
                ids.update(b.names)
            if page.kind == "module":
                ids.update(n for b in self.blocks if b.module == page.module for n in b.names)
            for nd in Node("root", self.page_content(page)).walk():
                if nd.get("id"):
                    ids.add(nd["id"])
            ids.update(("references", "used-by", "guides", "examples"))
            page._anchors = ids
        return page._anchors

    def _check_links(self):
        for page in self.order:
            for nd in Node("root", self.page_content(page)).walk():
                if nd.kind == "cite":
                    for key in nd["keys"]:
                        if key not in self.bib:
                            self.warn("unknown citation [@%s]" % key, page.path, nd.get("line"))
                    continue
                if nd.kind == "div" and nd.get("kind") == "summary":
                    t = (nd.get("args") or "").split()
                    if t and self.resolve_page(t[0]) is None:
                        self.warn("unknown page in ::: summary %s" % t[0], page.path, nd.get("line"))
                    continue
                if nd.kind != "link":
                    continue
                if self.resolve(nd["target"], page)[2] == "unresolved":
                    self.warn("unresolved link target %r" % nd["target"], page.path, nd.get("line"))
                    continue
                target = nd["target"]
                if "://" in target or target.startswith("mailto:") or "#" not in target:
                    continue
                base, anchor = target.split("#", 1)
                tp = page if not base else self.resolve_page(base)
                if tp is page and page.module in self.modules and anchor not in self.anchors(page):
                    tp = self.modules[page.module]
                if tp is not None and anchor not in self.anchors(tp):
                    self.warn("link to %s: no anchor #%s on page %s" % (target, anchor, tp.id),
                              page.path, nd.get("line"))

    def resolve(self, target, page):
        """Resolve a link target. Returns (url_relative_to_site_root or None,
        default link text nodes or None, kind). Absolute URLs are returned
        with kind 'external'."""
        if "://" in target or target.startswith("mailto:"):
            return target, None, "external"
        anchor = ""
        if "#" in target:
            target, anchor = target.split("#", 1)
        if not target and anchor not in self.anchors(page) and page.module in self.modules:
            # a same-page link in a block's text that is shown on a function
            # page: it refers to the module overview the block comes from
            mpage = self.modules[page.module]
            if anchor in self.anchors(mpage):
                page = mpage
        if not target:
            label = page.numbers.get(anchor)
            text = [Node("text", [], text="%s %d" % label)] if label else page.headings.get(anchor)
            return anchor_url(page.url, anchor), text, "anchor"
        name = target[:-2] if target.endswith("()") else target
        if name in self.symbols:
            sym = self.symbols[name]
            return sym.url, [Node("code", [], text=name)], "symbol"
        tp = self.resolve_page(target)
        if tp is not None:
            if anchor:
                label = tp.numbers.get(anchor)
                if label and tp is not page:
                    # "Theorem 2" alone would be read as this page's Theorem 2
                    text = [Node("text", [], text="%s %d of " % label)] + list(tp.title)
                elif label:
                    text = [Node("text", [], text="%s %d" % label)]
                else:
                    text = tp.headings.get(anchor)
                return anchor_url(tp.url, anchor), text, "page"
            return tp.url, tp.title, "page"
        return None, None, "unresolved"

    def source_url(self, name):
        loc = self.sources.locate(name)
        if not loc:
            return None
        rel, line = loc
        base = self.config.get("github", "https://github.com/flintlib/flint")
        return "%s/blob/%s/%s#L%d" % (base, self.commit[:12] if len(self.commit) == 40 else self.commit,
                                     rel, line)

    # -- checks ----------------------------------------------------------------

    def configured_features(self):
        """Features of the configured build, from src/flint-config.h: the
        number of bits of a limb ("64" or "32") and "fft_small"."""
        feats = set()
        try:
            with open(os.path.join(self.flintroot, "src", "flint-config.h")) as f:
                text = f.read()
        except OSError:
            return feats
        m = re.search(r"^#define FLINT_BITS (\d+)", text, re.M)
        if m:
            feats.add(m.group(1))
        if re.search(r"^#define FLINT_HAVE_FFT_SMALL 1", text, re.M):
            feats.add("fft_small")
        return feats

    def _check_sources(self):
        self.undocumented = {}
        # a module can only be checked against a configuration that has it:
        # front matter 'requires: 64' (or fft_small, ...) in its overview
        feats = self.configured_features()
        self.unchecked = {}
        for mod, page in self.modules.items():
            req = [r.strip() for r in str(page.meta.get("requires", "")).split(",") if r.strip()]
            unmet = [r for r in req if r not in feats]
            if unmet:
                self.unchecked[mod] = unmet
                self.undocumented[mod] = []
                continue
            documented = set(b_name for b in self.blocks if b.module == mod for b_name in b.names)
            for b in self.blocks:
                if b.module != mod or b.kind != "function":
                    continue
                for p in b.protos:
                    decl = self.sources.declarations.get(p.name)
                    if decl is None:
                        if self.sources.locate(p.name) is None and p.name not in self.sources.macros:
                            self.warn("%s is documented but not found in the source" % p.name,
                                      b.path, b.line)
                        continue
                    if p.params is not None and decl[2].signature_key() != p.signature_key():
                        self.warn("prototype of %s differs from %s:%d\n    doc:    %s\n    header: %s"
                                  % (p.name, decl[0], decl[1], p.text, decl[2].text), b.path, b.line)
            hdr = "src/%s.h" % page.meta.get("source", mod)
            missing = sorted(n for n, d in self.sources.declarations.items()
                             if d[0] in (hdr, hdr + ".in") and n not in documented and not n.startswith("_"))
            self.undocumented[mod] = missing


def _parse_file(path):
    """Read and parse one source file: ("ok", meta, blocks), or
    ("error", message, path, line). Runs in a worker process."""
    with open(path, encoding="utf-8") as f:
        text = f.read()
    try:
        meta, body, off = split_front_matter(text)
        return ("ok", meta, BlockParser(path).parse(body, off))
    except MarkupError as e:
        return ("error", e.msg, e.path, e.line)


def _collect_module(module):
    """SourceIndex.collect in a worker process."""
    from .parallel import state
    return state().collect(module)


def slugify(s):
    s = re.sub(r"[^\w\s-]", "", s.lower())
    return re.sub(r"[\s_]+", "-", s).strip("-")
