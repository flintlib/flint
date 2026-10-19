"""
The fdoc markup language: a small, strict Markdown dialect.

Block syntax
    # Heading {#optional-id}
    paragraphs, - bullet / 1. numbered lists, | pipe | tables |
    ```c [example] [name=foo]      fenced code (``example`` = compiled by make checkdoc)
    ```output [flags]              expected output of the preceding example
                                   (flags: radius-tol[=F], rel-tol=X, abs-tol=X,
                                   random; see outputs.py)
    $$ ... $$                      display math
    ::: kind [args] [{#id}]        fenced div (function, details, theorem, proof, ...)
    :::                            closes the innermost div

Inline syntax
    `code`  (auto-linked if it names a documented symbol)
    $math$  *em*  **strong**  [text](target)  [](target)  [@BibKey]  [@A; @B]
    bare https:// URLs, \\ escapes, -- (en dash)

Underscore emphasis is deliberately not supported, so identifiers and math
never need escaping.  The parser produces a plain tree of Node objects that
is shared by all output backends (HTML, LaTeX).
"""

import re


class Node:
    __slots__ = ("kind", "children", "attrs")

    def __init__(self, node_kind, children=None, **attrs):
        self.kind = node_kind
        self.children = children if children is not None else []
        self.attrs = attrs

    def __getitem__(self, key):
        return self.attrs[key]

    def get(self, key, default=None):
        return self.attrs.get(key, default)

    def __repr__(self):
        a = " ".join("%s=%r" % kv for kv in self.attrs.items())
        return "<%s %s %r>" % (self.kind, a, self.children)

    def walk(self):
        yield self
        for c in self.children:
            if isinstance(c, Node):
                yield from c.walk()
        for v in self.attrs.values():
            if isinstance(v, Node):
                yield from v.walk()
            elif isinstance(v, list):
                for x in v:
                    if isinstance(x, Node):
                        yield from x.walk()
                    elif isinstance(x, list):
                        for y in x:
                            if isinstance(y, Node):
                                yield from y.walk()


class MarkupError(Exception):
    def __init__(self, msg, path=None, line=None):
        self.msg, self.path, self.line = msg, path, line
        loc = "%s:%s: " % (path, line) if path else ""
        super().__init__(loc + msg)


# ---------------------------------------------------------------------------
# Front matter
# ---------------------------------------------------------------------------

def split_front_matter(text):
    """Return (dict, body, body_line_offset). Values are strings; lists are
    written either inline (a, b, c) or as indented '- item' lines."""
    meta = {}
    if not text.startswith("---\n"):
        return meta, text, 0
    end = text.find("\n---\n", 3)
    if end < 0:
        raise MarkupError("unterminated front matter")
    block = text[4:end].split("\n")
    key = None
    for ln in block:
        if not ln.strip() or ln.lstrip().startswith("#"):
            continue
        m = re.match(r"^([A-Za-z][\w-]*):\s*(.*)$", ln)
        if m:
            key = m.group(1)
            meta[key] = m.group(2).strip()
        elif key and ln.startswith((" ", "\t")):
            item = ln.strip()
            if item.startswith("- "):
                if not isinstance(meta[key], list):
                    meta[key] = [] if not meta[key] else [meta[key]]
                meta[key].append(item[2:].strip())
            else:
                meta[key] = (meta[key] + " " + item).strip()
    body = text[end + 5:]
    return meta, body, text[:end + 5].count("\n")


def meta_list(value):
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [x for x in re.split(r"[,\s]+", value) if x]


# ---------------------------------------------------------------------------
# Block parser
# ---------------------------------------------------------------------------

RE_FENCE = re.compile(r"^(`{3,})\s*(.*)$")
RE_DIV_OPEN = re.compile(r"^(:{3,})\s*([A-Za-z][\w-]*)\s*(.*?)\s*$")
RE_DIV_CLOSE = re.compile(r"^:{3,}\s*$")
RE_HEADING = re.compile(r"^(#{1,6})\s+(.*?)\s*$")
RE_ID = re.compile(r"\s*\{#([\w.:-]+)\}\s*$")
RE_ITEM = re.compile(r"^(\s*)([-*]|\d+[.)])\s+(.*)$")


def split_cells(ln):
    """Split a table row at '|', except inside $...$, `...` and after a
    backslash."""
    parts, cur, math, code, i = [], [], False, False, 0
    while i < len(ln):
        ch = ln[i]
        if ch == "\\" and i + 1 < len(ln):
            cur.append(ln[i:i + 2]); i += 2; continue
        if ch == "`" and not math:
            code = not code
        elif ch == "$" and not code:
            math = not math
        elif ch == "|" and not math and not code:
            parts.append("".join(cur)); cur = []; i += 1; continue
        cur.append(ch); i += 1
    parts.append("".join(cur))
    return parts


def in_inline_math(lines):
    """True if the paragraph lines so far end inside an unclosed $...$, so
    that a wrapped line such as '- i}}$' continues the formula."""
    text = re.sub(r"`[^`]*`", "", " ".join(lines)).replace("\\$", "")
    return text.count("$") % 2 == 1


def interrupts_para(line):
    """A list item may interrupt a paragraph only if it is a bullet or starts
    at 1 (as in CommonMark), so that a wrapped line such as '0. The next
    sentence' or '2023. In ...' stays part of the paragraph."""
    m = RE_ITEM.match(line)
    return bool(m) and (not m.group(2)[0].isdigit() or m.group(2)[:-1] == "1")
RE_TABLE_SEP = re.compile(r"^\s*\|?\s*:?-{2,}:?\s*(\|\s*:?-{2,}:?\s*)*\|?\s*$")


class BlockParser:
    def __init__(self, path=None):
        self.path = path

    def error(self, msg, line):
        raise MarkupError(msg, self.path, line)

    def parse(self, text, line0=0):
        lines = text.expandtabs(4).split("\n")
        return self.blocks(lines, 0, len(lines), line0)

    # lines[i:j] -> list of block nodes
    def blocks(self, lines, i, j, line0):
        out = []
        while i < j:
            ln = lines[i]
            s = ln.strip()
            if not s:
                i += 1
                continue

            m = RE_FENCE.match(s)
            if m:
                i = self.code_block(lines, i, j, line0, m, out)
                continue

            if s.startswith("$$"):
                i = self.math_block(lines, i, j, line0, out)
                continue

            m = RE_DIV_OPEN.match(s)
            if m:
                i = self.div(lines, i, j, line0, m, out)
                continue
            if RE_DIV_CLOSE.match(s):
                self.error("unmatched ':::'", line0 + i + 1)

            m = RE_HEADING.match(s)
            if m:
                level = len(m.group(1))
                title = m.group(2)
                hid = None
                mid = RE_ID.search(title)
                if mid:
                    hid = mid.group(1)
                    title = title[:mid.start()]
                out.append(Node("heading", parse_inline(title), level=level,
                                id=hid, line=line0 + i + 1))
                i += 1
                continue

            if s.startswith("|") and i + 1 < j and RE_TABLE_SEP.match(lines[i + 1]):
                i = self.table(lines, i, j, line0, out)
                continue

            if RE_ITEM.match(ln):
                i = self.list_block(lines, i, j, line0, out)
                continue

            # paragraph
            k = i
            buf = []
            while k < j:
                t = lines[k]
                ts = t.strip()
                if not ts:
                    break
                if k > i and (RE_FENCE.match(ts) or ts.startswith("$$") or
                              RE_DIV_OPEN.match(ts) or RE_DIV_CLOSE.match(ts) or
                              RE_HEADING.match(ts) or
                              (interrupts_para(t) and not in_inline_math(buf))):
                    break
                buf.append(ts)
                k += 1
            out.append(Node("para", parse_inline(" ".join(buf)), line=line0 + i + 1))
            i = k
        return out

    def code_block(self, lines, i, j, line0, m, out):
        ticks, info = m.group(1), m.group(2).split()
        lang = info[0] if info else ""
        flags = {}
        for w in info[1:]:
            if "=" in w:
                a, b = w.split("=", 1)
                flags[a] = b
            else:
                flags[w] = True
        indent = len(lines[i]) - len(lines[i].lstrip())
        k = i + 1
        body = []
        while k < j and lines[k].strip() != ticks:
            body.append(lines[k][indent:] if lines[k][:indent].strip() == "" else lines[k])
            k += 1
        if k >= j:
            self.error("unterminated code block", line0 + i + 1)
        node = Node("code", [], lang=lang, text="\n".join(body), flags=flags,
                    line=line0 + i + 1)
        if lang == "output" and out and out[-1].kind == "code" and out[-1]["flags"].get("example"):
            out[-1].attrs["output"] = node["text"]
            out[-1].attrs["output_flags"] = flags
        else:
            out.append(node)
        return k + 1

    def math_block(self, lines, i, j, line0, out):
        s = lines[i].strip()
        if s.endswith("$$") and len(s) > 4:
            out.append(Node("displaymath", [], tex=s[2:-2].strip(), line=line0 + i + 1))
            return i + 1
        buf = [s[2:]]
        k = i + 1
        while k < j and not lines[k].strip().endswith("$$"):
            buf.append(lines[k])
            k += 1
        if k >= j:
            self.error("unterminated $$", line0 + i + 1)
        buf.append(lines[k].strip()[:-2])
        out.append(Node("displaymath", [], tex="\n".join(x.strip() for x in buf).strip(),
                        line=line0 + i + 1))
        return k + 1

    def find_div_end(self, lines, i, j, line0):
        depth = 1
        k = i + 1
        fence = None
        while k < j:
            s = lines[k].strip()
            if fence:
                if s == fence:
                    fence = None
            else:
                fm = RE_FENCE.match(s)
                if fm:
                    fence = fm.group(1)
                elif RE_DIV_OPEN.match(s):
                    if not s.endswith(":::") or s.count(":::") < 2:
                        depth += 1
                elif RE_DIV_CLOSE.match(s):
                    depth -= 1
                    if depth == 0:
                        return k
            k += 1
        self.error("unterminated ':::' block", line0 + i + 1)

    def div(self, lines, i, j, line0, m, out):
        kind, rest = m.group(2), m.group(3)
        single = rest.endswith(":::")
        if single:
            rest = rest[:-3].strip()
        did = None
        mid = RE_ID.search(rest)
        if mid:
            did = mid.group(1)
            rest = rest[:mid.start()].strip()
        if single:
            out.append(Node("div", [], kind=kind, args=rest, id=did, line=line0 + i + 1))
            return i + 1
        k = self.find_div_end(lines, i, j, line0)
        inner_lines = lines[i + 1:k]
        node = Node("div", [], kind=kind, args=rest, id=did, line=line0 + i + 1)
        if kind in ("function", "type", "macro", "variable"):
            # leading non-blank lines are prototypes
            p = 0
            protos = []
            cur = ""
            while p < len(inner_lines) and inner_lines[p].strip():
                cur = (cur + " " + inner_lines[p].strip()).strip()
                if cur.count("(") <= cur.count(")"):
                    protos.append(cur)
                    cur = ""
                p += 1
            if cur:
                protos.append(cur)
            node.attrs["protos"] = protos
            node.children = self.blocks(lines, i + 1 + p, k, line0)
        else:
            node.children = self.blocks(lines, i + 1, k, line0)
        out.append(node)
        return k + 1

    def table(self, lines, i, j, line0, out):
        def cells(ln):
            ln = ln.strip()
            if ln.startswith("|"):
                ln = ln[1:]
            if ln.endswith("|") and not ln.endswith("\\|"):
                ln = ln[:-1]
            parts = split_cells(ln)
            return [parse_inline(p.strip().replace("\\|", "|")) for p in parts]
        head = cells(lines[i])
        aligns = []
        for c in lines[i + 1].strip().strip("|").split("|"):
            c = c.strip()
            aligns.append("center" if c.startswith(":") and c.endswith(":")
                          else "right" if c.endswith(":") else "left")
        rows = []
        k = i + 2
        while k < j and lines[k].strip().startswith("|"):
            rows.append(cells(lines[k]))
            k += 1
        out.append(Node("table", [], head=head, aligns=aligns, rows=rows, line=line0 + i + 1))
        return k

    def list_block(self, lines, i, j, line0, out):
        m = RE_ITEM.match(lines[i])
        base = len(m.group(1))
        ordered = m.group(2)[0].isdigit()
        start = int(m.group(2)[:-1]) if ordered else 1
        items = []
        k = i
        while k < j:
            m = RE_ITEM.match(lines[k])
            if not m or len(m.group(1)) != base or m.group(2)[0].isdigit() != ordered:
                break
            content_indent = len(m.group(1)) + len(m.group(2)) + 1
            item_lines = [" " * content_indent + m.group(3)]
            first_line = k
            k += 1
            while k < j:
                t = lines[k]
                if not t.strip():
                    # blank: continue only if next nonblank is indented
                    q = k + 1
                    while q < j and not lines[q].strip():
                        q += 1
                    if q < j and (len(lines[q]) - len(lines[q].lstrip())) > base:
                        item_lines.extend([""] * (q - k))
                        k = q
                        continue
                    break
                ind = len(t) - len(t.lstrip())
                if ind <= base and (RE_ITEM.match(t) or RE_DIV_CLOSE.match(t.strip())
                                    or RE_DIV_OPEN.match(t.strip())):
                    break
                if ind <= base and RE_FENCE.match(t.strip()):
                    break
                item_lines.append(t)
                k += 1
            # dedent
            ind = min((len(t) - len(t.lstrip()) for t in item_lines if t.strip()), default=0)
            ind = min(ind, content_indent)
            sub = [t[ind:] if t[:ind].strip() == "" else t.lstrip() for t in item_lines]
            items.append(self.blocks(sub, 0, len(sub), line0 + first_line))
            # skip blank lines between items
            q = k
            while q < j and not lines[q].strip():
                q += 1
            if q < j and RE_ITEM.match(lines[q]) and len(RE_ITEM.match(lines[q]).group(1)) == base:
                k = q
        out.append(Node("list", [], ordered=ordered, start=start, items=items,
                        line=line0 + i + 1))
        return k


# ---------------------------------------------------------------------------
# Inline parser
# ---------------------------------------------------------------------------

RE_URL = re.compile(r"https?://[^\s<>()\]]*[^\s<>()\].,;:!?'\"]")


def parse_inline(s):
    out = []
    buf = []

    def flush():
        if buf:
            t = "".join(buf).replace("---", "—").replace("--", "–")
            out.append(Node("text", [], text=t))
            buf.clear()

    i = 0
    n = len(s)
    while i < n:
        c = s[i]
        if c == "\\" and i + 1 < n and not s[i + 1].isalnum():
            buf.append(s[i + 1])
            i += 2
            continue
        if c == "`":
            m = re.match(r"`+", s[i:])
            ticks = m.group(0)
            end = s.find(ticks, i + len(ticks))
            if end > 0:
                flush()
                out.append(Node("code", [], text=s[i + len(ticks):end].strip()))
                i = end + len(ticks)
                continue
        if c == "$":
            k = i + 1
            while k < n and not (s[k] == "$" and s[k - 1] != "\\"):
                k += 1
            if k < n and k > i + 1:
                flush()
                out.append(Node("math", [], tex=s[i + 1:k]))
                i = k + 1
                continue
        if c == "*":
            dbl = s.startswith("**", i)
            mark = "**" if dbl else "*"
            end = _find_close(s, i + len(mark), mark)
            if end > 0:
                flush()
                out.append(Node("strong" if dbl else "em",
                                parse_inline(s[i + len(mark):end])))
                i = end + len(mark)
                continue
        if c == "[":
            if s.startswith("[@", i):
                end = s.find("]", i)
                if end > 0:
                    flush()
                    keys, locs = [], []
                    for k in s[i + 1:end].split(";"):
                        k, _, loc = k.strip().lstrip("@").partition(",")
                        keys.append(k.strip())
                        locs.append(loc.strip())
                    out.append(Node("cite", [], keys=keys, locs=locs))
                    i = end + 1
                    continue
            end = _match_bracket(s, i)
            if end > 0 and end + 1 < n and s[end + 1] == "(":
                close = s.find(")", end + 2)
                if close > 0:
                    flush()
                    target = s[end + 2:close].strip()
                    out.append(Node("link", parse_inline(s[i + 1:end]), target=target))
                    i = close + 1
                    continue
        if c == "h" and (s.startswith("http://", i) or s.startswith("https://", i)):
            m = RE_URL.match(s, i)
            if m:
                flush()
                url = m.group(0)
                out.append(Node("link", [Node("text", [], text=url)], target=url, bare=True))
                i = m.end()
                continue
        buf.append(c)
        i += 1
    flush()
    return out


def _find_close(s, i, mark):
    # closing marker must not be preceded by whitespace
    k = i
    while True:
        k = s.find(mark, k)
        if k < 0:
            return -1
        if k > i and not s[k - 1].isspace() and (mark == "**" or not s.startswith("**", k)):
            return k
        if mark == "*" and s.startswith("**", k):
            # skip over a strong span
            e = s.find("**", k + 2)
            if e < 0:
                return -1
            k = e + 2
            continue
        k += len(mark)


def _match_bracket(s, i):
    depth = 0
    k = i
    while k < len(s):
        if s[k] == "\\":
            k += 2
            continue
        if s[k] == "`":
            e = s.find("`", k + 1)
            if e < 0:
                return -1
            k = e + 1
            continue
        if s[k] == "[":
            depth += 1
        elif s[k] == "]":
            depth -= 1
            if depth == 0:
                return k
        k += 1
    return -1


def plain_text(nodes):
    """Flatten inline nodes to plain text (for search index, titles)."""
    out = []
    for nd in nodes:
        if nd.kind == "text":
            out.append(nd["text"])
        elif nd.kind == "code":
            out.append(nd["text"])
        elif nd.kind == "math":
            out.append(tex_to_text(nd["tex"]))
        elif nd.kind == "cite":
            out.append("[" + ", ".join(nd["keys"]) + "]")
        elif nd.kind == "para":
            out.append(plain_text(nd.children))
        else:
            out.append(plain_text(nd.children))
    return "".join(out)


TEX_SYMBOLS = {"bmod": "mod", "pmod": "mod", "lfloor": "\u230a", "rfloor": "\u230b",
               "lceil": "\u2308", "rceil": "\u2309", "le": "\u2264", "leq": "\u2264",
               "ge": "\u2265", "geq": "\u2265", "ne": "\u2260", "cdot": "\u00b7",
               "times": "\u00d7", "ldots": "\u2026", "dots": "\u2026", "mid": "|",
               "nmid": "\u2224", "equiv": "\u2261", "to": "\u2192", "pi": "\u03c0",
               "mu": "\u03bc", "varphi": "\u03c6", "sqrt": "\u221a", "in": "\u2208",
               "bits": "FLINT_BITS", "infty": "\u221e", "pm": "\u00b1"}


def tex_to_text(tex):
    """Rough plain-text rendering of a formula (for search snippets)."""
    def cmd(m):
        name = m.group(1)
        return TEX_SYMBOLS.get(name, "" if name in ("left", "right", "mathit", "mathrm",
                                                     "operatorname", "mathtt", "frac", "text") else name)
    t = re.sub(r"\\([A-Za-z]+)", cmd, tex)
    t = t.replace("\\_", "_").replace("\\,", " ").replace("{", "").replace("}", "")
    return " ".join(t.split())
