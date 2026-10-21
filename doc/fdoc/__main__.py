"""
fdoc — the FLINT documentation builder.

    python3 -m fdoc build    [--out build/html]          HTML site
    python3 -m fdoc latex    [--out build/latex] [--pdf] LaTeX (and PDF)
    python3 -m fdoc examples [--cc cc] [--cflags ...] [--libs ...]
                                                         compile + run every example
    python3 -m fdoc check                                consistency report

Run from the doc/ directory (or pass --src / --flint).
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


# headers that ./configure generates and that the header checks need
CONFIGURED_HEADERS = ("flint.h", "flint-config.h", "flint-mparam.h", "gmpcompat.h")


def load_site(args):
    from .model import Site
    src = os.path.abspath(args.src)
    flint = os.path.abspath(args.flint)
    missing = [h for h in CONFIGURED_HEADERS if not os.path.exists(os.path.join(flint, "src", h))]
    if missing:
        sys.exit("fdoc: %s not found in %s: run ./configure first"
                 % (", ".join(missing), os.path.join(flint, "src")))
    try:
        with open(os.path.join(flint, "VERSION")) as f:
            version = f.read().strip()
    except OSError:
        version = ""
    # the generation date: SOURCE_DATE_EPOCH (reproducible builds) or today
    import datetime
    epoch = os.environ.get("SOURCE_DATE_EPOCH")
    date = (datetime.datetime.fromtimestamp(int(epoch), datetime.timezone.utc).date()
            if epoch else datetime.date.today())
    config = {"project": "FLINT", "version": version,
              "github": "https://github.com/flintlib/flint",
              "baseurl": "https://flintlib.org/doc/",
              "date": date.isoformat()}
    site = Site(src, flint, config)
    site.load()
    return site


def only_filter(args):
    """Path prefixes (relative to src/) selected with --only, or None."""
    only = getattr(args, "only", None)
    if not only:
        return None
    return [os.path.join(os.path.abspath(args.src), x.strip().rstrip("/")) for x in only.split(",") if x.strip()]


def selected(path, prefixes):
    if prefixes is None:
        return True
    path = os.path.abspath(path) if path else ""
    return any(path == p or path.startswith(p + os.sep) or path.startswith(p + ".") for p in prefixes)


def report(site, strict=False, args=None):
    prefixes = only_filter(args) if args is not None else None
    items = [w for w in site.warn.items
             if prefixes is None or any(p in w for p in prefixes)]
    for w in items:
        print("warning:", w, file=sys.stderr)
    return 1 if (strict and items) else 0


def cmd_build(args):
    from .html import HtmlBuilder
    t0 = time.time()
    site = load_site(args)
    if getattr(args, "bundle", False):
        site.bundle_auto_pages()
    n = HtmlBuilder(site, os.path.abspath(args.out)).build()
    if getattr(args, "redirects", False):
        from .redirects import write_redirects
        r = write_redirects(site, os.path.abspath(args.out))
        print("wrote %d redirect pages for the URLs of the old documentation" % r)
    dt = time.time() - t0
    print("built %d pages in %.2f s -> %s" % (n, dt, args.out))
    return report(site, args.strict, args)


def cmd_latex(args):
    from .latex import LatexBuilder
    from .progress import Progress
    site = load_site(args)
    out = os.path.abspath(args.out)
    pr = Progress("writing LaTeX")
    path = LatexBuilder(site, out).build()
    pr.done()
    print("wrote", path)
    rc = report(site, args.strict, args)
    if args.pdf:
        if shutil.which("pdflatex") is None:
            print("pdflatex not found: install a TeX distribution to build the PDF")
            return 1
        for i in range(2):
            pr = Progress("pdflatex, pass %d of 2" % (i + 1))
            res = subprocess.run(["pdflatex", "-interaction=nonstopmode", "-halt-on-error",
                                  os.path.basename(path)], cwd=out, capture_output=True)
            pr.done()
            if res.returncode != 0:
                break
        if res.returncode != 0:
            print(res.stdout.decode()[-3000:])
            return 1
        print("wrote", path[:-4] + ".pdf")
    return rc


def cmd_examples(args):
    """Extract every '```c example' block, compile it against FLINT, run it,
    and compare stdout with the '```output' block that follows it."""
    from . import outputs
    site = load_site(args)
    outdir = os.path.abspath(args.out)
    os.makedirs(outdir, exist_ok=True)
    cc = args.cc.split()
    flint = os.path.abspath(args.flint)
    if args.cflags is None:
        # in-tree build: programs include <flint/xxx.h>, so provide a
        # directory in which flint/ points to the source headers
        inc = os.path.join(outdir, "_include")
        os.makedirs(inc, exist_ok=True)
        link = os.path.join(inc, "flint")
        if not os.path.islink(link):
            os.symlink(os.path.join(flint, "src"), link)
        args.cflags = "-std=c11 -Wall -I%s -I%s -I%s" % (inc, os.path.join(flint, "src"), flint)
    if args.libs is None:
        args.libs = "-L%s -lflint -lmpfr -lgmp -lm -Wl,-rpath,%s" % (flint, flint)
    cflags = args.cflags.split()
    libs = args.libs.split()
    # features of the configured build, for examples marked requires=...
    # (64, 32, fft_small; several separated by commas)
    probe = os.path.join(outdir, "_features_probe.c")
    with open(probe, "w") as f:
        f.write('#include <stdio.h>\n#include <flint/flint.h>\n'
                'int main(void)\n{\n    printf("%d", FLINT_BITS);\n'
                '#if FLINT_HAVE_FFT_SMALL\n    printf(" fft_small");\n#endif\n'
                '    printf("\\n");\n    return 0;\n}\n')
    res = subprocess.run(cc + cflags + ["-o", probe[:-2], probe] + libs, capture_output=True)
    if res.returncode != 0:
        print("fdoc: cannot compile a test program against FLINT:\n" + res.stderr.decode())
        return 1
    features = set(subprocess.run([probe[:-2]], capture_output=True).stdout.decode().split())
    # requirements of every example, set by the documentation in the front
    # matter of the home page ('example-requires'): the expected outputs are
    # those of a build with these features
    home = site.pages.get("index")
    default_req = str(home.meta.get("example-requires", "")) if home else ""
    prefixes = only_filter(args)
    jobs = []
    for page in site.order:
        if not selected(page.path, prefixes):
            continue
        for ex in page.examples:
            jobs.append((page, ex))

    def run(job):
        """Compile and run one example; returns (status, name, where, message)."""
        page, ex = job
        name = ex["example_name"]
        where = "%s:%d" % (os.path.relpath(page.path, site.srcroot), ex["line"])
        req = list(dict.fromkeys(r.strip() for r in (default_req + "," + str(ex["flags"].get("requires", "")))
                                 .split(",") if r.strip()))
        unmet = [r for r in req if r not in features]
        if unmet:
            return "skip", name, where, "needs a build with " + ", ".join(
                ("%s-bit limbs" % r) if r.isdigit() else r for r in unmet)
        src = os.path.join(outdir, name + ".c")
        exe = os.path.join(outdir, name)
        with open(src, "w") as f:
            f.write(ex["text"].rstrip() + "\n")
        res = subprocess.run(cc + cflags + ["-o", exe, src] + libs, capture_output=True)
        if res.returncode != 0:
            return "fail", name, where, "compile\n" + res.stderr.decode()
        try:
            res = subprocess.run([exe], capture_output=True, timeout=args.timeout)
        except subprocess.TimeoutExpired:
            return "fail", name, where, "timeout"
        got = res.stdout.decode().rstrip()
        if res.returncode != 0:
            return "fail", name, where, "exit status %d" % res.returncode
        want = ex.get("output")
        if want is None:
            return "ran", name, where, "no expected output"
        if not outputs.compare(want, got, ex.get("output_flags")):
            return "fail", name, where, "output\n--- expected\n%s\n--- got\n%s" % (want.rstrip(), got)
        return "ok", name, where, ""

    from concurrent.futures import ThreadPoolExecutor
    from .progress import Progress
    from .parallel import jobs as _jobs
    njobs = _jobs()
    pr = Progress("compiling and running examples", len(jobs))
    results = []
    with ThreadPoolExecutor(max_workers=njobs) as pool:
        for r in pool.map(run, jobs):
            results.append(r)
            pr.step()
    failures = sum(1 for r in results if r[0] == "fail")
    skipped = sum(1 for r in results if r[0] == "skip")
    pr.done("%d failed" % failures if failures else "")
    reasons = {}
    for status, name, where, msg in results:
        if status == "skip":
            reasons[msg] = reasons.get(msg, 0) + 1
    for msg, n in sorted(reasons.items()):
        print("skipped %d example%s: %s" % (n, "s" if n > 1 else "", msg))
    for status, name, where, msg in results:
        if status == "fail":
            print("FAIL %s (%s): %s" % (name, where, msg))
        elif args.verbose:
            print("%-4s %-28s %s%s" % (status, name, where, " (%s)" % msg if msg else ""))
    ok = sum(1 for r in results if r[0] == "ok")
    ran = sum(1 for r in results if r[0] == "ran")
    print("%d examples: %d compiled, ran and printed the expected output, %d compiled and ran "
          "(no expected output given), %d failed, %d skipped"
          % (len(results), ok, ran, failures, skipped))
    return 1 if failures else 0


def cmd_check(args):
    site = load_site(args)
    prefixes = only_filter(args)
    for mod, missing in sorted(site.undocumented.items()):
        if prefixes is not None and not selected(os.path.join(site.srcroot, mod), prefixes):
            continue
        # a function declared in one module's header may be documented in
        # another module (e.g. nmod_mat_charpoly is declared in nmod_poly.h)
        missing = [name for name in missing if name not in site.symbols]
        if missing:
            print("%s: %d public functions declared in the header are not documented:" % (mod, len(missing)))
            for name in missing:
                print("    " + name)
    for mod, unmet in sorted(getattr(site, "unchecked", {}).items()):
        print("%s: not checked against the headers: needs a build with %s"
              % (mod, ", ".join(("%s-bit limbs" % r) if r.isdigit() else r for r in unmet)))
    nblocks = len(site.blocks)
    ndetail = sum(1 for p in site.order if p.kind == "detail")
    nex = sum(len(p.examples) for p in site.order)
    nalg = sum(1 for p in site.order if p.kind in ("algorithm", "topic"))
    print("%d symbols in %d blocks; %d hand-written function pages; %d examples; %d algorithm/topic pages"
          % (len(site.symbols), nblocks, ndetail, nex, nalg))
    return report(site, True, args)


def main(argv=None):
    # progress goes to stderr: keep stdout line-buffered so that reports
    # and progress lines appear in order when both go to the same log
    try:
        sys.stdout.reconfigure(line_buffering=True)
    except (AttributeError, ValueError):
        pass
    ap = argparse.ArgumentParser(prog="fdoc")
    ap.add_argument("--src", default=os.path.join(os.path.dirname(HERE), "src"))
    ap.add_argument("--flint", default=os.path.dirname(os.path.dirname(HERE)))
    # the general options are accepted before and after the command, so that
    # e.g. 'make doc DOCFLAGS=--strict' works (the defaults are set here; the
    # subcommands' copies only override them when given)
    def general(parser, suppress):
        d = (lambda v: argparse.SUPPRESS) if suppress else (lambda v: v)
        parser.add_argument("--strict", action="store_true", default=d(False),
                            help="treat warnings as errors")
        parser.add_argument("-q", "--quiet", action="store_true", default=d(False),
                            help="no progress output")
        parser.add_argument("-j", "--jobs", type=int, default=d(0),
                            help="parallel jobs (default: the number of CPUs)")
    general(ap, False)
    common = argparse.ArgumentParser(add_help=False)
    general(common, True)
    sub = ap.add_subparsers(dest="cmd", required=True)
    _add_parser = sub.add_parser
    sub.add_parser = lambda name, **kw: _add_parser(name, parents=[common], **kw)
    p = sub.add_parser("build")
    p.add_argument("--out", default="build/html")
    p.add_argument("--redirects", action="store_true",
                   help="also write <old>.html pages redirecting the URLs of the old (Sphinx) documentation")
    p.add_argument("--bundle", action="store_true",
                   help="serve each module's automatic function pages from one shell page")
    p = sub.add_parser("latex")
    p.add_argument("--out", default="build/latex")
    p.add_argument("--pdf", action="store_true")
    p = sub.add_parser("examples")
    p.add_argument("--out", default="build/examples")
    p.add_argument("--cc", default=os.environ.get("CC", "cc"))
    p.add_argument("--cflags", default=None,
                   help="default: the in-tree headers and build of --flint")
    p.add_argument("--libs", default=None,
                   help="default: the in-tree libflint of --flint")
    p.add_argument("--timeout", type=float, default=20)
    p.add_argument("-v", "--verbose", action="store_true", help="also list the examples that pass")
    p.add_argument("--only", help="comma-separated paths under src/ (e.g. fmpz_poly,algorithms/polynomials)")
    p = sub.add_parser("check")
    p.add_argument("--only", help="comma-separated paths under src/: report only these")
    args = ap.parse_args(argv)
    if args.quiet:
        from . import progress
        progress.QUIET = True
    from . import parallel
    parallel.JOBS = args.jobs
    return {"build": cmd_build, "latex": cmd_latex, "examples": cmd_examples,
            "check": cmd_check}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
