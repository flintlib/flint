"""
C helpers: prototype parsing, a tiny tokenizer for highlighting, and a
scanner that locates definitions/declarations in the FLINT source tree.
"""

import os
import re

C_KEYWORDS = set("""auto break case const continue default do else enum extern for goto if
inline register restrict return sizeof static struct switch typedef union volatile while
_Bool _Complex""".split())

C_TYPES = set("""void char short int long float double signed unsigned size_t
ulong slong ulong_ptr nn_ptr nn_srcptr mp_limb_t mp_size_t mp_ptr mp_srcptr
flint_bitcnt_t flint_rand_t fmpz_t fmpz fmpq_t n_factor_t n_primes_t n_primes_struct
n_ecm_t n_ecm_s nmod_t int64_t uint64_t int32_t uint32_t FILE""".split())

QUALIFIERS = {"ULONG_EXTRAS_INLINE", "FLINT_DLL", "FLINT_FORCE_INLINE", "FLINT_CONST",
              "FLINT_WARN_UNUSED", "static", "inline", "extern", "FLINT_INLINE",
              "NMOD_INLINE", "FMPZ_INLINE", "WARN_UNUSED_RESULT", "FLINT_NORETURN",
              "FLINT_MALLOC", "FLINT_RETURNS_NONNULL", "FLINT_DEPRECATED", "FLINT_WARN_UNUSED_RESULT"}

RE_TOKEN = re.compile(r"""
    (?P<comment>/\*.*?\*/|//[^\n]*)
  | (?P<pre>^[ \t]*\#[^\n]*)
  | (?P<string>"(?:[^"\\\n]|\\.)*"|'(?:[^'\\\n]|\\.)*')
  | (?P<number>\b(?:0[xX][0-9a-fA-F]+|\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)[uUlL]*\b)
  | (?P<ident>[A-Za-z_]\w*)
  | (?P<space>\s+)
  | (?P<punct>.)
""", re.S | re.X | re.M)


def tokenize(src):
    for m in RE_TOKEN.finditer(src):
        kind = m.lastgroup
        yield kind, m.group(kind)


class Prototype:
    """A parsed C prototype such as 'ulong n_gcd(ulong x, ulong y)'."""

    def __init__(self, text):
        self.text = " ".join(re.sub(r"FLINT_UNUSED\((\w+(?:\s*\[\s*\])?)\)", r"\1", text).split())
        t = self.text.rstrip(";")
        p = t.find("(")
        if p < 0:
            # a type, variable or object-like macro
            words = t.split()
            self.name = words[-1].lstrip("*")
            self.ret = " ".join(words[:-1])
            self.params = None
            return
        head = t[:p].rstrip()
        m = re.search(r"([A-Za-z_]\w*)$", head)
        if not m:
            raise ValueError("cannot parse prototype: %r" % text)
        self.name = m.group(1)
        self.ret = head[:m.start()].strip()
        inner = t[p + 1:t.rfind(")")]
        self.params = split_params(inner)

    def param_types(self):
        if self.params is None:
            return None
        return tuple(normalize_type(param_type(x)) for x in self.params if x.strip() not in ("", "void"))

    def signature_key(self):
        return (normalize_type(self.ret), self.param_types())


def split_params(inner):
    out, depth, cur = [], 0, []
    for ch in inner:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    if "".join(cur).strip():
        out.append("".join(cur).strip())
    return out


def param_type(p):
    """'const ulong * x' -> 'const ulong *'; 'ulong ** sqrt' -> 'ulong **'."""
    p = p.strip()
    if p in ("...",):
        return p
    p = re.sub(r"\[[^\]]*\]$", "*", p)
    m = re.match(r"^(.*?[\s*])([A-Za-z_]\w*)$", p)
    if m and m.group(1).strip() and m.group(2) not in C_TYPES:
        return m.group(1)
    return p


def is_qualifier(w):
    return w in QUALIFIERS or (w.isupper() and w.endswith(("_INLINE", "_DLL", "_DEF")))


def normalize_type(t):
    words = [w for w in re.findall(r"[A-Za-z_]\w*|\*", t) if not is_qualifier(w)]
    s = " ".join(words)
    s = s.replace("mp_limb_t", "ulong").replace("unsigned long", "ulong")
    s = s.replace("mp_ptr", "nn_ptr").replace("mp_srcptr", "nn_srcptr")
    return s


# ---------------------------------------------------------------------------
# Source scanning
# ---------------------------------------------------------------------------

RE_DEF = re.compile(
    r"^(?P<head>[A-Za-z_][\w \t*]*?\n?[ \t]*\**[ \t]*)?\b(?P<name>[A-Za-z_]\w*)[ \t]*\((?P<params>[^;{}()]*(?:\([^;{}()]*\)[^;{}()]*)*)\)[ \t]*\n?[ \t]*\{",
    re.M)
RE_DECL = re.compile(
    r"^(?P<head>(?:[A-Za-z_][\w]*[ \t*]+)+\**)(?P<name>[A-Za-z_]\w*)[ \t]*\((?P<params>[^;{}()]*(?:\([^;{}()]*\)[^;{}()]*)*)\)[ \t]*;",
    re.M)
RE_MACRO = re.compile(r"^[ \t]*#[ \t]*define[ \t]+(?P<name>[A-Za-z_]\w*)\b", re.M)
NOT_FUNCS = {"if", "while", "for", "switch", "return", "sizeof"}


def _strip_comments(text):
    # keep line structure so that line numbers stay correct
    return re.sub(r"/\*.*?\*/", lambda m: re.sub(r"[^\n]", " ", m.group(0)), text, flags=re.S)


class SourceIndex:
    """Locates where symbols are defined (for 'view source' links) and how
    they are declared in headers (for prototype checking)."""

    def __init__(self, srcdir):
        self.srcdir = srcdir
        self.definitions = {}   # name -> (relpath, line)
        self.declarations = {}  # name -> (relpath, line, Prototype)
        self.macros = {}
        self.generated = set()      # declared only after macro expansion
        self.preprocess_flags = {}  # module -> list of compiler flag strings
        self.template_headers = {}  # module -> extra headers whose expansion counts as the module's
        self.scanned = set()

    def scan_module(self, module):
        if module in self.scanned:
            return
        self.scanned.add(module)
        self._scan_module(module)
        self.drop_macro_names()
        self.scan_preprocessed(module, self.preprocess_flags.get(module))

    def _scan_module(self, module):
        header = os.path.join(self.srcdir, module + ".h")
        files = []
        mdir = os.path.join(self.srcdir, module)
        if os.path.isdir(mdir):
            for fn in sorted(os.listdir(mdir)):
                if fn.endswith((".c", ".h")):
                    files.append(os.path.join(mdir, fn))
        if os.path.exists(header):
            files.append(header)
        for path in files:
            rel = os.path.relpath(path, os.path.dirname(self.srcdir))
            with open(path, encoding="utf-8", errors="replace") as f:
                text = _strip_comments(f.read())
            for m in RE_DEF.finditer(text):
                name = m.group("name")
                if name in NOT_FUNCS:
                    continue
                line = text.count("\n", 0, m.start("name")) + 1
                prev = self.definitions.get(name)
                # prefer a .c file named after the function, then any .c file
                if prev is None or (prev[0].endswith(".h") and path.endswith(".c")):
                    self.definitions[name] = (rel, line)
            if path == header or path.endswith(".h"):
                for m in RE_DECL.finditer(text):
                    name = m.group("name")
                    if name in C_KEYWORDS or name in C_TYPES or m.group("head").lstrip().startswith("typedef"):
                        continue
                    line = text.count("\n", 0, m.start("name")) + 1
                    try:
                        proto = Prototype(m.group("head") + name + "(" + m.group("params") + ")")
                    except ValueError:
                        continue
                    self.declarations.setdefault(name, (rel, line, proto))
                for m in RE_DEF.finditer(text):
                    name = m.group("name")
                    if name in NOT_FUNCS or not m.group("head"):
                        continue
                    if name in C_KEYWORDS or name in C_TYPES or m.group("head").lstrip().startswith("typedef"):
                        continue
                    line = text.count("\n", 0, m.start("name")) + 1
                    try:
                        proto = Prototype(m.group("head") + name + "(" + m.group("params") + ")")
                    except ValueError:
                        continue
                    self.declarations.setdefault(name, (rel, line, proto))
                for m in RE_MACRO.finditer(text):
                    line = text.count("\n", 0, m.start("name")) + 1
                    self.macros.setdefault(m.group("name"), (rel, line))

    def scan_preprocessed(self, module, flag_sets=None):
        """Run the C preprocessor over the module header and record
        declarations that only exist after macro expansion (functions
        generated by macros).  Optional: silently skipped without a
        compiler."""
        import shutil
        import subprocess
        header = os.path.join(self.srcdir, module + ".h")
        cc = os.environ.get("CC", "cc")
        if not os.path.exists(header) or not shutil.which(cc.split()[0]):
            return
        for flags in (flag_sets or [""]):
            try:
                stubs = os.path.join(os.path.dirname(os.path.abspath(__file__)), "stubs")
                flags = flags.replace("{stubs}", stubs)
                res = subprocess.run(cc.split() + flags.split() + ["-E", "-I" + self.srcdir, header],
                                     capture_output=True, timeout=60)
            except Exception:
                continue
            if res.returncode == 0:
                # keep only the text that comes from this module's own headers
                mdir = os.path.join(self.srcdir, module) + os.sep
                extra = [os.path.normpath(os.path.join(self.srcdir, h if h.endswith(".h") else h + ".h"))
                         for h in self.template_headers.get(module, [])]
                keep, own = [], False
                for ln in res.stdout.decode(errors="replace").split("\n"):
                    m = re.match(r'^# \d+ "([^"]*)"', ln)
                    if m:
                        f = os.path.normpath(m.group(1))
                        own = (f == os.path.normpath(header) or f.startswith(os.path.normpath(mdir))
                               or f in extra)
                        continue
                    if own:
                        keep.append(ln)
                self._scan_preprocessed_text("\n".join(keep), header)

    def _scan_preprocessed_text(self, text, header):
        text = re.sub(r"__attribute__\s*\(\((?:[^()]|\([^()]*\))*\)\)", " ", text)
        text = re.sub(r"\b(__inline__|__inline|__extension__|__restrict)\b", " ", text)
        # one declaration per statement: join lines, split on ; and {
        text = re.sub(r"\s+", " ", text)
        rel = os.path.relpath(header, os.path.dirname(self.srcdir))
        for m in re.finditer(r"(?:^|(?<=[;}{]))\s*((?:[A-Za-z_]\w*[\s*]+)+)([A-Za-z_]\w*)\s*\(([^;{}()]*(?:\([^()]*\)[^;{}()]*)*)\)\s*[;{]",
                             text):
            name = m.group(2)
            if name in NOT_FUNCS or name in self.declarations:
                continue
            head = m.group(1)
            if head.strip().split()[0] in ("return", "typedef", "else"):
                continue
            try:
                proto = Prototype(head + name + "(" + m.group(3) + ")")
            except ValueError:
                continue
            self.declarations[name] = (rel, 0, proto)
            self.generated.add(name)

    def drop_macro_names(self):
        for name in list(self.declarations):
            if name in self.macros and name not in self.generated:
                del self.declarations[name]
        for name in list(self.definitions):
            if name in self.macros:
                del self.definitions[name]

    def locate(self, name):
        if name in self.definitions:
            return self.definitions[name]
        if name in self.declarations:
            d = self.declarations[name]
            if d[1] == 0:
                return None
            return d[0], d[1]
        return self.macros.get(name)
