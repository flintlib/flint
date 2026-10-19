"""
LaTeX backend: the whole site as one book-style document.

Layout: the user guide; a chapter of modules with one section per module
(its overview, with every function's short description, its guides and the
function pages that have extended content); the algorithms, the examples
and tutorials and the history; then the bibliography.  Theorem-like
blocks keep the numbers they have in the HTML output.
"""

import os
import re

from .markup import Node, plain_text
from .outputs import describe as describe_output
from .model import THEOREM_KINDS

PREAMBLE = r"""\documentclass[11pt,a4paper]{report}
\usepackage[utf8]{inputenc}
\DeclareUnicodeCharacter{03B1}{\ensuremath{\alpha}}
\DeclareUnicodeCharacter{03B2}{\ensuremath{\beta}}
\DeclareUnicodeCharacter{03B3}{\ensuremath{\gamma}}
\DeclareUnicodeCharacter{03B4}{\ensuremath{\delta}}
\DeclareUnicodeCharacter{03B5}{\ensuremath{\varepsilon}}
\DeclareUnicodeCharacter{03B6}{\ensuremath{\zeta}}
\DeclareUnicodeCharacter{03B7}{\ensuremath{\eta}}
\DeclareUnicodeCharacter{03B8}{\ensuremath{\theta}}
\DeclareUnicodeCharacter{03B9}{\ensuremath{\iota}}
\DeclareUnicodeCharacter{03BA}{\ensuremath{\kappa}}
\DeclareUnicodeCharacter{03BB}{\ensuremath{\lambda}}
\DeclareUnicodeCharacter{03BC}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{03BD}{\ensuremath{\nu}}
\DeclareUnicodeCharacter{03BE}{\ensuremath{\xi}}
\DeclareUnicodeCharacter{03C0}{\ensuremath{\pi}}
\DeclareUnicodeCharacter{03C1}{\ensuremath{\rho}}
\DeclareUnicodeCharacter{03C3}{\ensuremath{\sigma}}
\DeclareUnicodeCharacter{03C4}{\ensuremath{\tau}}
\DeclareUnicodeCharacter{03C5}{\ensuremath{\upsilon}}
\DeclareUnicodeCharacter{03C6}{\ensuremath{\varphi}}
\DeclareUnicodeCharacter{03C7}{\ensuremath{\chi}}
\DeclareUnicodeCharacter{03C8}{\ensuremath{\psi}}
\DeclareUnicodeCharacter{03C9}{\ensuremath{\omega}}
\DeclareUnicodeCharacter{0393}{\ensuremath{\Gamma}}
\DeclareUnicodeCharacter{0394}{\ensuremath{\Delta}}
\DeclareUnicodeCharacter{0398}{\ensuremath{\Theta}}
\DeclareUnicodeCharacter{039B}{\ensuremath{\Lambda}}
\DeclareUnicodeCharacter{039E}{\ensuremath{\Xi}}
\DeclareUnicodeCharacter{03A0}{\ensuremath{\Pi}}
\DeclareUnicodeCharacter{03A3}{\ensuremath{\Sigma}}
\DeclareUnicodeCharacter{03A6}{\ensuremath{\Phi}}
\DeclareUnicodeCharacter{03A8}{\ensuremath{\Psi}}
\DeclareUnicodeCharacter{03A9}{\ensuremath{\Omega}}
\DeclareUnicodeCharacter{00B1}{\ensuremath{\pm}}
\DeclareUnicodeCharacter{00D7}{\ensuremath{\times}}
\DeclareUnicodeCharacter{00B7}{\ensuremath{\cdot}}
\DeclareUnicodeCharacter{2260}{\ensuremath{\neq}}
\DeclareUnicodeCharacter{2264}{\ensuremath{\leq}}
\DeclareUnicodeCharacter{2265}{\ensuremath{\geq}}
\DeclareUnicodeCharacter{2192}{\ensuremath{\to}}
\DeclareUnicodeCharacter{2190}{\ensuremath{\leftarrow}}
\DeclareUnicodeCharacter{221E}{\ensuremath{\infty}}
\DeclareUnicodeCharacter{221A}{\ensuremath{\surd}}
\DeclareUnicodeCharacter{2211}{\ensuremath{\sum}}
\DeclareUnicodeCharacter{2208}{\ensuremath{\in}}
\DeclareUnicodeCharacter{2248}{\ensuremath{\approx}}
\DeclareUnicodeCharacter{2261}{\ensuremath{\equiv}}
\DeclareUnicodeCharacter{2223}{\ensuremath{\mid}}
\DeclareUnicodeCharacter{2308}{\ensuremath{\lceil}}
\DeclareUnicodeCharacter{2309}{\ensuremath{\rceil}}
\DeclareUnicodeCharacter{230A}{\ensuremath{\lfloor}}
\DeclareUnicodeCharacter{230B}{\ensuremath{\rfloor}}
\DeclareUnicodeCharacter{00B2}{\ensuremath{^2}}
\DeclareUnicodeCharacter{00B3}{\ensuremath{^3}}
\DeclareUnicodeCharacter{2080}{\ensuremath{_0}}
\DeclareUnicodeCharacter{2081}{\ensuremath{_1}}
\DeclareUnicodeCharacter{2082}{\ensuremath{_2}}
\DeclareUnicodeCharacter{2124}{\ensuremath{\mathbb{Z}}}
\DeclareUnicodeCharacter{211A}{\ensuremath{\mathbb{Q}}}
\DeclareUnicodeCharacter{211D}{\ensuremath{\mathbb{R}}}
\DeclareUnicodeCharacter{2102}{\ensuremath{\mathbb{C}}}
\DeclareUnicodeCharacter{2115}{\ensuremath{\mathbb{N}}}
\DeclareUnicodeCharacter{2212}{\ensuremath{-}}
\DeclareUnicodeCharacter{00B5}{\ensuremath{\mu}}
\DeclareUnicodeCharacter{00B0}{\ensuremath{^\circ}}
\DeclareUnicodeCharacter{2218}{\ensuremath{\circ}}
\DeclareUnicodeCharacter{2286}{\ensuremath{\subseteq}}
\DeclareUnicodeCharacter{2229}{\ensuremath{\cap}}
\DeclareUnicodeCharacter{222A}{\ensuremath{\cup}}
\DeclareUnicodeCharacter{2026}{\ldots}
\DeclareUnicodeCharacter{00B6}{\P}
\IfFileExists{lmodern.sty}{\usepackage[T1]{fontenc}\usepackage{lmodern}}{}
\usepackage[expansion=false]{microtype}
\usepackage{amsmath,amssymb,mathrsfs}
\usepackage[margin=2.6cm]{geometry}
\usepackage{xcolor}
\usepackage{fancyvrb}
\usepackage{tabularx,booktabs,xltabular}
\newcolumntype{Y}{>{\raggedright\arraybackslash}X}
\usepackage{enumitem}
\usepackage[hidelinks]{hyperref}
\usepackage{url}
\definecolor{fdocmuted}{gray}{0.4}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.55em}
\setlist{itemsep=0.15em,topsep=0.3em}
\newcommand{\fdocproto}[1]{\par{\raggedright\hangindent=2em\hangafter=1\noindent\ttfamily #1\par}}
\newenvironment{fdocdesc}{\par\vspace{-0.3em}\begin{list}{}{\leftmargin=2em\topsep=0pt}\item[]}{\end{list}}
\newenvironment{fdocenv}[1]{\par\medskip\noindent\textbf{#1.}\ }{\par\medskip}
\newenvironment{fdocproof}{\par\noindent\textit{Proof.}\ }{\hfill$\square$\par\medskip}
\newcommand{\fdocsep}{\par\vspace{0.2em}{\color{fdocmuted}\rule{\linewidth}{0.2pt}}\par}
\newcommand{\us}{\textunderscore\allowbreak{}}
\makeatletter
% display-style (not run-in) paragraph headings, so that a table or a list
% following a deep heading starts on its own line
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}{2.5ex \@plus 1ex \@minus .2ex}{0.6ex}{\normalfont\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\z@}{2ex \@plus .8ex \@minus .2ex}{0.5ex}{\normalfont\normalsize\itshape}}
% wider number boxes in the table of contents (section numbers such as 3.155)
\renewcommand*\l@section{\@dottedtocline{1}{1.5em}{3.4em}}
\renewcommand*\l@subsection{\@dottedtocline{2}{4.9em}{4.4em}}
\renewcommand*\l@subsubsection{\@dottedtocline{3}{9.3em}{5.4em}}
\makeatother
%MACROS%
\title{%TITLE%}
\author{The FLINT development team}
\begin{document}
\maketitle
\tableofcontents
"""

SPECIAL = {"&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_", "{": r"\{",
           "}": r"\}", "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
           "\\": r"\textbackslash{}", "<": r"\textless{}", ">": r"\textgreater{}"}


def tex_escape(s):
    return "".join(SPECIAL.get(c, c) for c in s)


def tt(s):
    out = []
    for c in s:
        if c == "_":
            out.append(r"\us ")
        elif c == " ":
            out.append(r"\ ")
        else:
            out.append(SPECIAL.get(c, c))
    return r"\texttt{" + "".join(out) + "}"


def label(s):
    return re.sub(r"[^A-Za-z0-9:.-]", "-", s)


class LatexRenderer:
    def __init__(self, builder, page, level):
        self.b = builder
        self.site = builder.site
        self.page = page
        self.level = level      # sectioning level of '##' headings on this page
        self.cites = builder.cites
        self.nest = 0           # > 0 inside lists and divs (no page-breaking tables there)

    # -- inline ---------------------------------------------------------
    def inline(self, nodes):
        out = []
        for nd in nodes:
            k = nd.kind
            if k == "text":
                out.append(tex_escape(nd["text"]))
            elif k == "code":
                out.append(self.code_span(nd["text"]))
            elif k == "math":
                out.append("$" + nd["tex"] + "$")
            elif k == "em":
                out.append(r"\emph{%s}" % self.inline(nd.children))
            elif k == "strong":
                out.append(r"\textbf{%s}" % self.inline(nd.children))
            elif k == "cite":
                parts = []
                for key, loc in zip(nd["keys"], nd.get("locs") or [""] * len(nd["keys"])):
                    if key not in self.cites:
                        self.cites.append(key)
                    parts.append(r"\cite[%s]{%s}" % (tex_escape(loc), key) if loc else r"\cite{%s}" % key)
                out.append("~".join(parts))
            elif k == "link":
                out.append(self.link(nd))
            else:
                out.append(self.inline(nd.children))
        return "".join(out)

    def code_span(self, text):
        name = text[:-2] if text.endswith("()") else text
        if name in self.site.symbols:
            return r"\hyperref[fn:%s]{%s}" % (label(name), tt(text))
        return tt(text)

    def link(self, nd):
        target = nd["target"]
        url, default, kind = self.site.resolve(target, self.page)
        if kind == "external":
            if nd.get("bare"):
                return r"\url{%s}" % url
            return r"\href{%s}{%s}" % (url.replace("%", r"\%").replace("#", r"\#"), self.inline(nd.children))
        text = self.inline(nd.children) if nd.children else (self.inline(default) if default else tex_escape(target))
        if kind == "symbol":
            name = target[:-2] if target.endswith("()") else target
            return r"\hyperref[fn:%s]{%s}" % (label(name), text)
        if kind in ("page", "anchor"):
            base, _, anchor = target.partition("#")
            tp = self.site.resolve_page(base) if base else self.page
            if not base and url:
                # resolve() may redirect a block-local anchor to the module overview
                tp = next((p for p in self.site.order if p.url == url.split("#")[0]), tp)
            if tp is None:
                return text
            lab = "pg:" + label(tp.id) + (":" + label(anchor) if anchor else "")
            return r"\hyperref[%s]{%s}" % (lab, text)
        return text

    # -- blocks -----------------------------------------------------------
    def blocks(self, nodes):
        out = []
        run = []        # consecutive '::: summary' pointers share one quote
        for nd in nodes:
            if nd.kind == "div" and nd.get("kind") == "summary":
                x = self.summary_pointer(nd)
                if x:
                    run.append(x)
                continue
            if run:
                out.append("\\begin{quote}\\small %s\\end{quote}" % "\\\\\n".join(run))
                run = []
            x = self.block(nd)
            if x:
                out.append(x)
        if run:
            out.append("\\begin{quote}\\small %s\\end{quote}" % "\\\\\n".join(run))
        return "\n\n".join(out)

    def summary_pointer(self, nd):
        """A '::: summary' card in print: a reference to the section, without
        the card's abstract."""
        tp = self.site.resolve_page(nd["args"].split()[0]) if nd["args"] else None
        if tp is None:
            return ""
        return ("$\\rightarrow$ \\hyperref[pg:%s]{\\textbf{%s}} (p.~\\pageref{pg:%s})"
                % (label(tp.id), self.inline(tp.title), label(tp.id)))

    def heading_cmd(self, level):
        cmds = ["chapter", "section", "subsection", "subsubsection", "paragraph", "subparagraph"]
        return cmds[min(self.level + level - 2, len(cmds) - 1)]

    def block(self, nd):
        k = nd.kind
        if k == "para":
            return self.inline(nd.children)
        if k == "heading":
            hid = nd["id"] or re.sub(r"\W+", "-", plain_text(nd.children).lower())
            return r"\%s{%s}\label{pg:%s:%s}" % (self.heading_cmd(nd["level"]), self.inline(nd.children),
                                                 label(self.page.id), label(hid))
        if k == "list":
            env = "enumerate" if nd["ordered"] else "itemize"
            self.nest += 1
            items = "\n".join(r"\item " + self.blocks(it) for it in nd["items"])
            self.nest -= 1
            return "\\begin{%s}\n%s\n\\end{%s}" % (env, items, env)
        if k == "displaymath":
            tex = nd["tex"]
            return "\\[\n%s\n\\]" % tex
        if k == "code":
            return self.code_block(nd)
        if k == "table":
            return self.table(nd)
        if k == "div":
            return self.div(nd)
        return ""

    def code_block(self, nd):
        text = nd["text"]
        if nd.get("example_name"):
            out = ("\\par\\noindent{\\small\\textbf{Example} \\texttt{%s.c} "
                   "{\\color{fdocmuted}(compiled and run by \\texttt{make checkdoc})}}\n"
                   "\\begin{Verbatim}[fontsize=\\small,frame=single,framesep=2mm]\n%s\n\\end{Verbatim}\n"
                   % (tex_escape(nd["example_name"]), text))
            if nd.get("output") is not None:
                note = describe_output(nd.get("output_flags"))
                out += ("\n\\noindent{\\small\\textbf{Output}%s}\n"
                        "\\begin{Verbatim}[fontsize=\\small,frame=leftline,framesep=2mm]\n%s\n\\end{Verbatim}\n"
                        % (" {\\color{fdocmuted}(%s)}" % tex_escape(note) if note else "", nd["output"]))
            return out
        return "\\begin{Verbatim}[fontsize=\\small]\n%s\n\\end{Verbatim}\n" % text

    def table(self, nd):
        spec = "@{}" + "".join(r">{\hsize=%.3f\hsize}Y" % w for w in self.column_widths(nd)) + "@{}"
        head = " & ".join(r"\textbf{%s}" % self.inline(c) for c in nd["head"]) + r" \\"
        rows = [" & ".join(self.inline(c) for c in r) + r" \\" for r in nd["rows"]]
        if self.nest == 0:
            # top-level tables may break across pages; the header is repeated
            return ("{\\small\\renewcommand{\\arraystretch}{1.25}\\begin{xltabular}{\\linewidth}{%s}\\toprule\n"
                    "%s\\midrule\\endhead\n\\bottomrule\\endlastfoot\n%s\n\\end{xltabular}}" % (spec, head, "\n".join(rows)))
        return ("{\\small\\renewcommand{\\arraystretch}{1.25}\\begin{tabularx}{\\linewidth}{%s}\\toprule\n%s\\midrule\n%s\n"
                "\\bottomrule\\end{tabularx}}" % (spec, head, "\n".join(rows)))

    @staticmethod
    def column_widths(nd):
        """Relative widths of the columns (summing to the number of columns,
        as tabularx requires): proportional to the longest cell of each
        column, but wide enough for its longest word, so that a column of
        short numbers does not take as much room as one of sentences."""
        rows = [nd["head"]] + list(nd["rows"])
        n = len(nd["head"])
        weights = []
        for j in range(n):
            texts = [plain_text(r[j]) if j < len(r) else "" for r in rows]
            longest = max((len(s) for s in texts), default=1)
            word = max((len(w) for s in texts for w in s.split()), default=1)
            weights.append(max(min(longest, 60), 1.3 * word, 3))
        total = sum(weights)
        return [n * w / total for w in weights]

    def proto(self, p, anchor):
        parts = []
        if p.ret:
            parts.append(tt(p.ret) + r"\ ")
        parts.append(r"\textbf{%s}" % tt(p.name))
        if p.params is not None:
            parts.append(tt("(") + (tt(", ") + r"\allowbreak ").join(tt(x) for x in p.params) + tt(")"))
        lab = r"\phantomsection\label{fn:%s}" % label(p.name) if anchor else ""
        return r"\fdocproto{%s%s}" % (lab, "".join(parts))

    def div(self, nd):
        self.nest += 1
        try:
            return self._div(nd)
        finally:
            self.nest -= 1

    def _div(self, nd):
        kind = nd["kind"]
        lab = r"\phantomsection\label{pg:%s:%s}" % (label(self.page.id), label(nd["id"])) if nd["id"] else ""
        if kind in THEOREM_KINDS:
            head = "%s %d" % (kind.capitalize(), nd["number"])
            title = r" (%s)" % self.inline(_parse(nd["args"])) if nd["args"] else ""
            body = self.blocks(nd.children)
            if kind in ("theorem", "lemma", "proposition", "corollary"):
                body = r"\itshape " + body
            return "%s\\begin{fdocenv}{%s%s}%s\\end{fdocenv}" % (lab, head, title, body)
        if kind == "proof":
            return "\\begin{fdocproof}%s\\end{fdocproof}" % self.blocks(nd.children)
        if kind == "algorithm":
            title = "Algorithm %d" % nd["number"] + (": " + self.inline(_parse(nd["args"])) if nd["args"] else "")
            return ("%s\\par\\medskip\\noindent\\fbox{\\begin{minipage}{\\dimexpr\\linewidth-2\\fboxsep-2\\fboxrule}\\textbf{%s}\\par\\smallskip\n%s"
                    "\\end{minipage}}\\par\\medskip" % (lab, title, self.blocks(nd.children)))
        if kind in ("note", "warning"):
            return "\\begin{quote}\\textbf{%s.} %s\\end{quote}" % (tex_escape(nd["args"] or kind.capitalize()),
                                                                  self.blocks(nd.children))
        if kind == "summary":
            x = self.summary_pointer(nd)
            return "\\begin{quote}\\small %s\\end{quote}" % x if x else ""
        if kind == "functions":
            out = []
            for name in nd["args"].split():
                sym = self.site.symbols.get(name)
                if sym:
                    p = next(p for p in sym.block.protos if p.name == name)
                    out.append(self.proto(p, False))
            return "\\begin{quote}%s\\end{quote}" % "\n".join(out)
        if kind in ("function", "type", "macro", "variable"):
            b = nd["block"]
            protos = "\n".join(self.proto(p, True) for p in b.protos)
            return "%s\n\\begin{fdocdesc}%s\\end{fdocdesc}\\fdocsep" % (protos, self.blocks(b.brief))
        return self.blocks(nd.children)


def _parse(s):
    from .markup import parse_inline
    return parse_inline(s)


class LatexBuilder:
    def __init__(self, site, outdir):
        self.site = site
        self.outdir = outdir
        self.cites = []
        os.makedirs(outdir, exist_ok=True)

    def page_label(self, page):
        return r"\label{pg:%s}" % label(page.id)


    def collection_chapters(self, colls):
        from .model import COLLECTIONS
        parts = []
        site = self.site
        for coll in colls:
            ctitle = COLLECTIONS[coll]
            root = site.collection_root(coll)
            if root is None:
                continue
            r = LatexRenderer(self, root, 2)
            sec = ["\\chapter{%s}%s\n%s" % (tex_escape(ctitle), self.page_label(root), r.blocks(root.body))]
            for g, pages in site.collection_groups(coll):
                level = 3
                if g is not None:
                    r = LatexRenderer(self, g, 3)
                    sec.append("\\section{%s}%s\n%s" % (r.inline(g.title), self.page_label(g), r.blocks(g.body)))
                    level = 4
                for p in pages:
                    r = LatexRenderer(self, p, level)
                    cmd = "subsection" if level == 4 else "section"
                    sec.append("\\%s{%s}%s\n%s" % (cmd, r.inline(p.title), self.page_label(p), r.blocks(p.body)))
            parts.append("\n\n".join(sec))
        return parts

    def build(self):
        site = self.site
        parts = []
        home = site.pages.get("index")
        if home:
            r = LatexRenderer(self, home, 1)
            parts.append("\\chapter*{Introduction}\\addcontentsline{toc}{chapter}{Introduction}\n"
                         + r.blocks([n for n in home.body if n.kind in ("para", "list")]))
        parts.extend(self.collection_chapters(["guide"]))
        # all modules form one chapter: one section per module, with its
        # overview's headings as subsections, then its guides and the
        # function pages with extended content
        mods = ["\\chapter{Modules}\\label{pg:modules}\n"
                "Each section documents one module (one header file): an overview with every "
                "function and its short description, the module's guides, and the function "
                "pages that have extended content."]
        for name, mpage in site.modules_in_order():
            r = LatexRenderer(self, mpage, 2)
            mods.append("\\section{%s}%s\n%s" % (tt(mpage.meta.get("header", name + ".h")),
                                                 self.page_label(mpage), r.blocks(mpage.body)))
            # guides: narrative pages, direct children of the module
            for g in site.module_guides(name):
                r = LatexRenderer(self, g, 3)
                mods.append("\\subsection{%s}%s\n%s" % (r.inline(g.title), self.page_label(g), r.blocks(g.body)))
            # function pages with extended content
            details = []
            for p in site.module_pages(name):
                if not p.body and not any(b.details for b in p.blocks):
                    continue
                r = LatexRenderer(self, p, 4)
                title = r.inline(p.title) if p.kind == "detail" else tt(", ".join(p.blocks[0].names))
                if len(p.blocks) > 1:
                    # keep each block's notes with its own prototypes
                    head = "\n\n".join("\n".join(r.proto(pr, False) for pr in b.protos)
                                         + ("\n\\begin{fdocdesc}%s\\end{fdocdesc}" % r.blocks(b.details)
                                            if b.details else "")
                                         for b in p.blocks)
                    text = r.blocks(p.body)
                else:
                    head = "\n".join(r.proto(pr, False) for pr in p.blocks[0].protos)
                    text = r.blocks(p.blocks[0].details + p.body)
                details.append("\\subsubsection{%s}%s\n%s\n\n%s" % (title, self.page_label(p), head, text))
            if details:
                mods.append("\\subsection{%s: function details}\n%s" % (tt(name), "\n\n".join(details)))
        parts.append("\n\n".join(mods))
        parts.extend(self.collection_chapters(["algorithms", "tutorials", "history"]))
        if self.cites:
            r = LatexRenderer(self, site.pages["references"], 1)
            items = "\n".join(r"\bibitem[%s]{%s} %s" % (k, k, r.inline(site.bib[k]))
                              for k in sorted(self.cites, key=str.lower) if k in site.bib)
            parts.append("\\begin{thebibliography}{XXXXXXXXXX}\n%s\n\\end{thebibliography}" % items)
        macros = "\n".join(r"\providecommand{%s}{}\renewcommand{%s}{%s}" % (k, k, v)
                           for k, v in site.macros.items())
        pre = PREAMBLE.replace("%MACROS%", macros).replace(
            "%TITLE%", "FLINT %s" % tex_escape(site.config.get("version", "")))
        path = os.path.join(self.outdir, "flint.tex")
        with open(path, "w", encoding="utf-8") as f:
            f.write(pre + "\n\n".join(parts) + "\n\\end{document}\n")
        return path
