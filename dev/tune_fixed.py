#!/usr/bin/env python3
"""Tune and emit the fully specialized per-size implementations of the
fixed module (run from the top-level FLINT directory, after building
libflint):

    python3 dev/tune_fixed.py FUNC N [options]

FUNC is one of exp, log1p, atan, trig (trig covers sin_cos and tan,
which share the half-angle reconstruction).  The script

  1. generates an out-of-tree source file containing one fully
     specialized candidate per reduction parameter r in [--rmin,
     --rmax] -- each candidate is the SAME body the production file
     will use, with the hand series built for exactly that r;
  2. builds it against the in-tree libflint and runs it: every
     candidate is validated against MPFR over --points adversarial
     inputs and timed (best of five);
  3. selects the fastest candidate whose measured error stays within
     --margin of the documented budget.  The margin exists because a
     sweep-sized sample UNDERESTIMATES the maximum error: candidates
     at half the budget have later blown it at scale (see the r table
     comment in tan_bitwise_rs.c);
  4. with --emit, writes src/fixed/FUNC_opt_N.c containing the chosen
     specialization under its public name (fixed_exp_opt_3 and so on).
     --pin R skips the sweep and emits at exactly r = R: running with
     --pin at the shipped r must reproduce the shipped file byte for
     byte, which is the reproducibility check.

Example:  python3 dev/tune_fixed.py exp 4 --pin 16 --emit
          python3 dev/tune_fixed.py trig 9 --rmin 14 --rmax 26 --emit
"""

import argparse
import importlib.util
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)

LICENSE = """/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
"""


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    sys.modules[name] = m
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------------
# per-function description: how to emit a series, how to emit a body,
# the public name, the error budget, and the MPFR reference
# ---------------------------------------------------------------------

class Func:
    def series(self, o, n, r, name):
        raise NotImplementedError

    def body(self, o, n, r, fn_name, series_name, static):
        raise NotImplementedError

    def includes(self):
        return ['#include "flint.h"', '#include "longlong.h"',
                '#include "mpn_extras.h"', '#include "fixed.h"',
                '#include "impl.h"']


class ExpFunc(Func):
    name = "exp"
    outputs = 1
    def __init__(self):
        self.g = load(os.path.join(HERE, "gen_fixed_exp_bitwise_small.py"),
                      "g_exp_body")
        self.hs = load(os.path.join(HERE, "gen_fixed_exp_hand_mixed.py"),
                       "g_exp_series")
    HAND_R = {1: 12, 2: 16}    # hand-written series, not regenerable
    def budget(self, r):
        return 9 * r + 100
    def series(self, o, n, r, name):
        if n in self.HAND_R:
            if r != self.HAND_R[n]:
                raise ValueError("exp n <= 2 series are hand-written "
                                 "for r = %d only" % self.HAND_R[n])
            src = open(os.path.join(HERE,
                                    "exp_series_hand.inc")).read()
            fn = "_fixed_exp_rs_opt_%d" % n
            i = src.rindex("/*", 0, src.index("%s(nn_ptr res" % fn))
            j = src.index("\n}\n", i) + 3
            body = src[i:j].replace(fn, name).replace(
                "void\n", "static void\n", 1)
            o("#ifndef FIXED_EXP_E2")
            for line in src.splitlines():
                if line.startswith("#define FIXED_EXP_"):
                    o(line)
            o("#endif")
            for line in body.rstrip("\n").split("\n"):
                o(line)
            return
        self.hs.emit(o, n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        body = self.g.gen(n, opt=True, fn_name=fn_name,
                series_name=series_name, r_const=r, static=static)
        if isinstance(body, str):
            body = body.split("\n")
        for line in body:
            o(line)


class GenOneFunc(Func):
    """atan and log1p: bodies from gen_one(o, n, opt=...), series from
    gen_fixed_trig_hand_mixed."""
    def __init__(self, body_gen, series_func):
        self.g = load(os.path.join(HERE, body_gen), "g_" + self.name)
        self.hs = load(os.path.join(HERE, "gen_fixed_trig_hand_mixed.py"),
                       "g_trig_series")
        self.series_func = series_func
    def series(self, o, n, r, name):
        self.hs.emit(o, self.series_func, n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        for line in self.g.gen_one_named(n, r, fn_name, series_name,
                static):
            o(line)


class AtanFunc(GenOneFunc):
    name = "atan"
    outputs = 1
    def __init__(self):
        GenOneFunc.__init__(self, "gen_fixed_atan_bitwise_small.py", "atan")
    def budget(self, r):
        return 4 * r + 64


class Log1pFunc(GenOneFunc):
    name = "log1p"
    outputs = 1
    def __init__(self):
        GenOneFunc.__init__(self, "gen_fixed_log1p_bitwise_small.py",
                            "atanh")
    def budget(self, r):
        return 3 * r + 64
    def body(self, o, n, r, fn_name, series_name, static):
        if n > 2:
            GenOneFunc.body(self, o, n, r, fn_name, series_name,
                            static)
            return

        # n <= 2 predate the body generator: the hand-written register
        # bodies live in dev/log1p_hand.inc and are carried verbatim
        if r != 16:
            raise SystemExit("log1p n <= 2 exists only at r = 16 "
                             "(hand-written; see dev/log1p_hand.inc)")

        src = open(os.path.join(HERE, "log1p_hand.inc")).read()
        i = src.index("/* n = %d:" % n)
        j = src.index("/* n = %d:" % (n + 1)) if n == 1 else len(src)
        body = src[i:j].rstrip("\n")
        body = body.replace("@FN@", fn_name)
        body = body.replace("@SERIES@", series_name)
        if static:
            body = body.replace("void\n" + fn_name,
                                "static void\n" + fn_name)
        for line in body.split("\n"):
            o(line)

    def series(self, o, n, r, name):
        # the n = 2 hand body carries its own inline atanh series
        if n == 2:
            self.hs.emit(o, "atanh", 2, 16, name)
        elif n >= 3:
            GenOneFunc.series(self, o, n, r, name)


class TrigFunc(Func):
    """sin_cos and tan through the shared half-angle reconstruction:
    n <= 4 fully in registers, larger n through the exported
    _fixed_tan_halfangle_mid helper."""
    name = "trig"
    outputs = 3
    def __init__(self):
        self.g = load(os.path.join(HERE,
            "gen_fixed_tan_halfangle_small.py"), "g_trig_body")
        self.hs = load(os.path.join(HERE, "gen_fixed_trig_hand_mixed.py"),
                       "g_trig_series")
    def budget(self, r):
        return 6 * r + 128        # sin/cos; tan's 8r + 256 is looser
    def budget_tan(self, r):
        return 8 * r + 256
    def series(self, o, n, r, name):
        self.hs.emit(o, "tan", n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        if n == 1:
            src = open(os.path.join(HERE, "trig_hand.inc")).read()
            body = src[src.index("/* n = 1"):]
            body = body.replace("@FN@", fn_name).replace("@R@", str(r))
            body = body.replace("@SERIES@", series_name)
            if static:
                body = body.replace("void\n" + fn_name,
                                    "static void\n" + fn_name)
            for line in body.rstrip("\n").split("\n"):
                o(line)
        elif n <= 4:
            for line in self.g.emit_named(n, r, fn_name, series_name,
                    static):
                o(line)
        else:
            o("static void" if static else "void")
            o("%s(nn_ptr ysin, nn_ptr ycos, nn_ptr ytan," % fn_name)
            o("    nn_srcptr x)")
            o("{")
            o("    _fixed_tan_halfangle_mid(ysin, ycos, ytan, x, %d, %d,"
              % (n, r))
            o("        %s);" % series_name)
            o("}")

FUNCS = {"exp": ExpFunc, "log1p": Log1pFunc, "atan": AtanFunc,
         "trig": TrigFunc}


# ---------------------------------------------------------------------
# harness generation
# ---------------------------------------------------------------------

def harness_source(F, n, rs, points):
    L = []
    o = L.append
    o("/* generated by dev/tune_fixed.py: %s, n = %d -- build out of"
      % (F.name, n))
    o("   tree, run, discard */")
    for inc in F.includes():
        o(inc)
    o('#include <stdio.h>')
    o('#include <time.h>')
    o('#include <mpfr.h>')
    o('#include "hand_mulhi.inc"')
    if F.name == "trig" and n <= 4:
        o('#include "tan_rotate.inc"')
    o("")
    for r in rs:
        F.series(o, n, r, "ts_%d" % r)
        F.body(o, n, r, "cand_%d" % r, "ts_%d" % r, static=True)
    o("")
    if F.outputs == 3:
        o("typedef void (*cand_fn)(nn_ptr, nn_ptr, nn_ptr, nn_srcptr);")
    else:
        o("typedef void (*cand_fn)(nn_ptr, nn_srcptr);")
    o("static const struct { cand_fn f; int r; } cands[] = {")
    for r in rs:
        o("    { cand_%d, %d }," % (r, r))
    o("};")
    o("""
static double now(void)
{ struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t);
  return t.tv_sec * 1e9 + t.tv_nsec; }

int main(void)
{
    flint_rand_t state;
    slong i, j, k;
    flint_rand_init(state);
""")
    o("    for (k = 0; k < (slong) (sizeof(cands)/sizeof(cands[0])); k++)")
    o("    {")
    o("        ulong x[%d], y[%d], yc[%d], yt[%d], xs[64][%d];"
      % (n, n + 2, n + 2, n + 2, n))
    o("        double e = 0, ec = 0, et = 0, te = 1e30, t0;")
    o("        mpfr_t xm, fm, f2, f3, vm, dm;")
    o("        mpz_t z;")
    o("        mpz_init(z);")
    o("        mpfr_init2(xm, %d); mpfr_init2(fm, %d);" % (64*n, 64*n+64))
    o("        mpfr_init2(f2, %d); mpfr_init2(f3, %d);" % (64*n+64, 64*n+64))
    o("        mpfr_init2(vm, %d); mpfr_init2(dm, %d);" % (64*n+64, 64*n+64))
    o("        for (i = 0; i < %d; i++)" % points)
    o("        {")
    o("            for (j = 0; j < %d; j++) x[j] = n_randlimb(state);" % n)
    o("            if (i %% 7 == 0) x[%d - 1] |= UWORD(3) << 62;" % n)
    o("            if (i %% 11 == 0) { for (j = 0; j < %d; j++)" % n)
    o("                x[j] = ~UWORD(0); }")
    o("            if (i %% 13 == 5) { for (j = 0; j < %d; j++)" % n)
    o("                x[j] = 0; x[0] = i; }")
    if F.outputs == 3:
        o("            cands[k].f(y, yc, yt, x);")
    else:
        o("            cands[k].f(y, x);")
    o("            mpz_import(z, %d, -1, sizeof(ulong), 0, 0, x);" % n)
    o("            mpfr_set_z_2exp(xm, z, -%d, MPFR_RNDN);" % (64*n))
    if F.name == "exp":
        o("            mpfr_exp(fm, xm, MPFR_RNDN);")
        vn = n + 1
    elif F.name == "log1p":
        o("            mpfr_log1p(fm, xm, MPFR_RNDN);")
        vn = n
    elif F.name == "atan":
        o("            mpfr_atan(fm, xm, MPFR_RNDN);")
        vn = n
    else:
        o("            mpfr_sin_cos(fm, f2, xm, MPFR_RNDN);")
        o("            mpfr_tan(f3, xm, MPFR_RNDN);")
        vn = n + 1
    def check(dst, ref, acc):
        o("            mpz_import(z, %d, -1, sizeof(ulong), 0, 0, %s);"
          % (vn, dst))
        o("            mpfr_set_z_2exp(vm, z, -%d, MPFR_RNDN);" % (64*n))
        o("            mpfr_sub(dm, vm, %s, MPFR_RNDA);" % ref)
        o("            mpfr_abs(dm, dm, MPFR_RNDU);")
        o("            mpfr_mul_2si(dm, dm, %d, MPFR_RNDU);" % (64*n))
        o("            { double u = mpfr_get_d(dm, MPFR_RNDU);")
        o("              %s = FLINT_MAX(%s, u); }" % (acc, acc))
    check("y", "fm", "e")
    if F.outputs == 3:
        check("yc", "f2", "ec")
        check("yt", "f3", "et")
    o("        }")
    o("        mpfr_clears(xm, fm, f2, f3, vm, dm, (mpfr_ptr) 0);")
    o("        mpz_clear(z);")
    o("        for (j = 0; j < 64; j++)")
    o("            for (i = 0; i < %d; i++) xs[j][i] = n_randlimb(state);"
      % n)
    o("        for (j = 0; j < 5; j++)")
    o("        {")
    o("            t0 = now();")
    if F.outputs == 3:
        o("            for (i = 0; i < 120000; i++)")
        o("                cands[k].f(y, yc, NULL, xs[i & 63]);")
    else:
        o("            for (i = 0; i < 120000; i++)")
        o("                cands[k].f(y, xs[i & 63]);")
    o("            te = FLINT_MIN(te, (now() - t0) / 120000);")
    o("        }")
    o('        printf("r %d  ns %.1f  err %.0f %.0f %.0f\\n",')
    o("            cands[k].r, te, e, ec, et);")
    o("        fflush(stdout);")
    o("    }")
    o("    flint_rand_clear(state);")
    o("    return 0;")
    o("}")
    return "\n".join(L) + "\n"


def emit_production(F, n, r, path):
    L = []
    o = L.append
    o(LICENSE)
    o("/* GENERATED by dev/tune_fixed.py (%s, n = %d, r = %d) -- do not"
      % (F.name, n, r))
    o("   edit directly; retune or re-emit with")
    o("       python3 dev/tune_fixed.py %s %d --pin %d --emit */"
      % (F.name, n, r))
    o("")
    for inc in F.includes():
        o(inc)
    o('#include "hand_mulhi.inc"')
    if F.name == "trig" and n <= 4:
        o('#include "tan_rotate.inc"')
    o("")
    o("#if FLINT_BITS == 64")
    o("")
    if F.name == "trig":
        F.series(o, n, r, "_fixed_tan_series_%d" % n)
        F.body(o, n, r, "_fixed_trig_opt_%d" % n,
               "_fixed_tan_series_%d" % n, static=False)
        o("")
        o("void")
        o("fixed_sin_cos_opt_%d(nn_ptr ysin, nn_ptr ycos, nn_srcptr x)" % n)
        o("{")
        o("    _fixed_trig_opt_%d(ysin, ycos, NULL, x);" % n)
        o("}")
        o("")
        o("void")
        o("fixed_tan_opt_%d(nn_ptr res, nn_srcptr x)" % n)
        o("{")
        o("    _fixed_trig_opt_%d(NULL, NULL, res, x);" % n)
        o("}")
    else:
        F.series(o, n, r, "_fixed_%s_series_%d" % (F.name, n))
        F.body(o, n, r, "fixed_%s_opt_%d" % (F.name, n),
               "_fixed_%s_series_%d" % (F.name, n), static=False)
    o("")
    o("#endif /* FLINT_BITS == 64 */")
    open(path, "w").write("\n".join(L).rstrip("\n") + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("func", choices=sorted(FUNCS))
    ap.add_argument("n", type=int)
    ap.add_argument("--rmin", type=int, default=4)
    ap.add_argument("--rmax", type=int, default=34)
    ap.add_argument("--points", type=int, default=6000)
    ap.add_argument("--margin", type=float, default=0.5,
        help="accept only errors below margin * budget")
    ap.add_argument("--pin", type=int, default=None)
    ap.add_argument("--emit", action="store_true")
    ap.add_argument("--no-run", action="store_true",
        help="with --pin --emit: emit without building the harness "
             "(bootstrap only, e.g. while the library itself does not "
             "link yet)")
    a = ap.parse_args()

    F = FUNCS[a.func]()
    n = a.n

    if a.pin is not None:
        rs = [a.pin]
    else:
        rs = list(range(a.rmin, a.rmax + 1))

    # drop candidates the series generator cannot build
    ok_rs = []
    for r in rs:
        try:
            F.series([].append, n, r, "probe")
            ok_rs.append(r)
        except Exception:
            pass
    rs = ok_rs
    if not rs:
        sys.exit("no buildable candidates in the given r range")

    if a.no_run:
        if a.pin is None or not a.emit:
            sys.exit("--no-run needs --pin and --emit")
        path = os.path.join(TOP, "src", "fixed",
                            "%s_opt_%d.c" % (a.func, n))
        emit_production(F, n, a.pin, path)
        print("emitted", path, "(unvalidated: rerun without --no-run "
              "once the library links)")
        return

    src = os.path.join("/tmp", "tune_%s_%d.c" % (a.func, n))
    binp = os.path.join("/tmp", "tune_%s_%d" % (a.func, n))
    open(src, "w").write(harness_source(F, n, rs, a.points))

    cmd = ["cc", "-O2", "-march=native", "-I", "src", "-I", ".",
           "-I", "src/fixed", "-o", binp, src,
           "-L", TOP, "-Wl,-rpath," + TOP,
           "-lflint", "-lmpfr", "-lgmp", "-lm"]
    subprocess.run(cmd, cwd=TOP, check=True)
    out = subprocess.run([binp], cwd=TOP, check=True,
                         capture_output=True, text=True).stdout
    sys.stderr.write(out)

    accepted = []
    for line in out.splitlines():
        w = line.split()
        r, ns = int(w[1]), float(w[3])
        errs = [float(x) for x in w[5:8]]
        bud = F.budget(r)
        lim = [a.margin * bud] * 3
        if a.func == "trig":
            lim[2] = a.margin * F.budget_tan(r)
        if any(e > l for e, l in zip(errs, lim)):
            continue
        accepted.append((r, ns, max(e / l for e, l in zip(errs, lim)
                                    if l > 0)))
    if not accepted:
        sys.exit("every candidate exceeded margin * budget: rerun with "
                 "more points to see, or widen --margin knowingly")
    # among near-ties (within 1.5% of the fastest accepted), take the
    # cleanest error: sweep-sized samples UNDERESTIMATE the maximum,
    # and candidates near the margin have blown the budget at scale
    fastest = min(ns for _, ns, _ in accepted)
    ties = [c for c in accepted if c[1] <= 1.015 * fastest]
    r, ns, _ = min(ties, key=lambda c: c[2])
    print("selected r = %d (%.1f ns)" % (r, ns))

    if a.emit:
        path = os.path.join(TOP, "src", "fixed",
                            "%s_opt_%d.c" % (a.func, n))
        emit_production(F, n, r, path)
        print("emitted", path)


if __name__ == "__main__":
    main()
