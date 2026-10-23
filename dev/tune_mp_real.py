#!/usr/bin/env python3
"""Tune and emit the fully specialized per-size implementations of the
mp_real module (run from the top-level FLINT directory, after building
libflint):

    python3 dev/tune_mp_real.py FUNC N [options]

FUNC is one of exp, log1p, atan and sin_cos, and N in 4 .. 10: the
dense-table kernels of dev/gen_mp_real_dense_bodies.py, whose candidates
are r = 7 L - 1 for L reduction levels with the hand series, and (exp,
atan, sin_cos) r = 32 resp. (exp) r = 63 for the hardcoded series of the
2^-32 resp. 2^-64 tier after five levels resp. nine and a top level
(tables: dev/gen_mp_real_dense_tables.py).  (The smaller sizes have the
table-driven kernels of src/mp_real/small.h, and tan takes the sine and
cosine kernels and a division.)  The script

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
  4. with --emit, writes src/mp_real/FUNC_opt_N.c containing the chosen
     specialization under its public name (_mp_real_exp_opt_4 and so on).
     --pin R skips the sweep and emits at exactly r = R: running with
     --pin at the shipped r must reproduce the shipped file byte for
     byte, which is the reproducibility check.

Example:  python3 dev/tune_mp_real.py exp 4 --pin 34 --emit
          python3 dev/tune_mp_real.py sin_cos 9 --rmin 27 --rmax 41 --emit
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
                '#include "mpn_extras.h"', '#include "mp_real.h"',
                '#include "impl.h"']


class ExpFunc(Func):
    name = "exp"
    outputs = 1
    def __init__(self):
        self.hs = load(os.path.join(HERE, "gen_mp_real_exp_hand_mixed.py"),
                       "g_exp_series")
        self.d = load(os.path.join(HERE, "gen_mp_real_dense_bodies.py"),
                      "g_dense_body")
    # the sizes of the dense-table bodies (gen_mp_real_dense_bodies.py):
    # r = 7 L - 1 for L reduction levels of seven bits
    DENSE = range(4, 11)
    def budget(self, r):
        return 9 * r + 100
    def series(self, o, n, r, name):
        L, top, rs = self.d.dense_levels(r)
        if rs:
            self.d.emit_rs_series(o, "exp", n, name)
            return
        self.hs.emit(o, n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        L, top, rs = self.d.dense_levels(r)
        for line in self.d.gen_exp(n, L, fn_name, series_name, static,
                                   top=top):
            o(line)


class AtanFunc(Func):
    name = "atan"
    outputs = 1
    # the sizes of the dense-table bodies, r = 7 L - 1
    DENSE = range(4, 11)
    def __init__(self):
        self.hs = load(os.path.join(HERE, "gen_mp_real_trig_hand_mixed.py"),
                       "g_trig_series")
        self.d = load(os.path.join(HERE, "gen_mp_real_dense_bodies.py"),
                      "g_dense_body")
    def series(self, o, n, r, name):
        L, top, rs = self.d.dense_levels(r)
        if top:
            raise ValueError("no top level for atan")
        if rs:
            self.d.emit_rs_series(o, "atan", n, name)
            return
        self.hs.emit(o, "atan", n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        L, top, rs = self.d.dense_levels(r)
        for line in self.d.gen_atan(n, L, fn_name, series_name, static):
            o(line)
    def budget(self, r):
        return 4 * r + 64


class Log1pFunc(Func):
    name = "log1p"
    outputs = 1
    # the sizes of the dense-table bodies, r = 7 L - 1
    DENSE = range(4, 11)
    def __init__(self):
        self.d = load(os.path.join(HERE, "gen_mp_real_dense_bodies.py"),
                      "g_dense_body")
    def budget(self, r):
        return 3 * r + 64
    def body(self, o, n, r, fn_name, series_name, static):
        for line in self.d.gen_log1p(n, (r + 1) // 7, fn_name,
                                     series_name, static):
            o(line)
    def series(self, o, n, r, name):
        if (r + 1) % 7 != 0:
            raise ValueError("dense log1p kernels take r = 7 L - 1")
        self.d.emit_log1p_series(o, n, r, name)


class SinCosFunc(Func):
    """sin and cos: the dense-table kernels (gen_mp_real_dense_bodies.py),
    emitted as sin_cos_opt_<n>.c"""
    name = "sin_cos"
    outputs = 2
    DENSE = range(4, 11)
    def __init__(self):
        self.hs = load(os.path.join(HERE, "gen_mp_real_trig_hand_mixed.py"),
                       "g_trig_series")
        self.d = load(os.path.join(HERE, "gen_mp_real_dense_bodies.py"),
                      "g_dense_body")
    def budget(self, r):
        return 6 * r + 128
    def series(self, o, n, r, name):
        if n not in self.DENSE:
            raise ValueError("sin_cos kernels: n in 4..10")
        L, top, rs = self.d.dense_levels(r)
        if top:
            raise ValueError("no top level for sin_cos")
        if rs:
            self.d.emit_rs_series(o, "sin_cos", n, name)
            return
        self.hs.emit_sincos(o, n, r, name)
    def body(self, o, n, r, fn_name, series_name, static):
        L, top, rs = self.d.dense_levels(r)
        for line in self.d.gen_sin_cos(n, L, fn_name, series_name, static):
            o(line)


FUNCS = {"exp": ExpFunc, "log1p": Log1pFunc, "atan": AtanFunc,
         "sin_cos": SinCosFunc}


# ---------------------------------------------------------------------
# harness generation
# ---------------------------------------------------------------------

def harness_source(F, n, rs, points):
    L = []
    o = L.append
    o("/* generated by dev/tune_mp_real.py: %s, n = %d -- build out of"
      % (F.name, n))
    o("   tree, run, discard */")
    for inc in F.includes():
        o(inc)
    o('#include <stdio.h>')
    o('#include <time.h>')
    o('#include <mpfr.h>')
    o('#include "hand_mulhi.inc"')
    o("")
    for r in rs:
        F.series(o, n, r, "ts_%d" % r)
        F.body(o, n, r, "cand_%d" % r, "ts_%d" % r, static=True)
    o("")
    if F.outputs == 2:
        o("typedef void (*cand_fn)(nn_ptr, nn_ptr, nn_srcptr);")
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
    if F.outputs == 2:
        o("            cands[k].f(y, yc, x);")
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
    if F.outputs == 2:
        check("yc", "f2", "ec")
    o("        }")
    o("        mpfr_clears(xm, fm, f2, f3, vm, dm, (mpfr_ptr) 0);")
    o("        mpz_clear(z);")
    o("        for (j = 0; j < 64; j++)")
    o("            for (i = 0; i < %d; i++) xs[j][i] = n_randlimb(state);"
      % n)
    o("        for (j = 0; j < 5; j++)")
    o("        {")
    o("            t0 = now();")
    if F.outputs == 2:
        o("            for (i = 0; i < 120000; i++)")
        o("                cands[k].f(y, yc, xs[i & 63]);")
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
    o("/* GENERATED by dev/tune_mp_real.py (%s, n = %d, r = %d) -- do not"
      % (F.name, n, r))
    o("   edit directly; retune or re-emit with")
    o("       python3 dev/tune_mp_real.py %s %d --pin %d --emit */"
      % (F.name, n, r))
    o("")
    for inc in F.includes():
        o(inc)
    o('#include "hand_mulhi.inc"')
    o("")
    o("#if FLINT_BITS == 64")
    o("")
    if F.name == "sin_cos":
        F.series(o, n, r, "_mp_real_sin_cos_series_%d" % n)
        F.body(o, n, r, "_mp_real_sin_cos_opt_%d" % n,
               "_mp_real_sin_cos_series_%d" % n, static=False)
    else:
        F.series(o, n, r, "_mp_real_%s_series_%d" % (F.name, n))
        F.body(o, n, r, "_mp_real_%s_opt_%d" % (F.name, n),
               "_mp_real_%s_series_%d" % (F.name, n), static=False)
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
    if n not in F.DENSE:
        sys.exit("the per-size kernels take n in 4 .. 10 (smaller sizes: "
                 "src/mp_real/small.h)")

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
        path = os.path.join(TOP, "src", "mp_real",
                            "%s_opt_%d.c" % (a.func, n))
        emit_production(F, n, a.pin, path)
        print("emitted", path, "(unvalidated: rerun without --no-run "
              "once the library links)")
        return

    src = os.path.join("/tmp", "tune_%s_%d.c" % (a.func, n))
    binp = os.path.join("/tmp", "tune_%s_%d" % (a.func, n))
    open(src, "w").write(harness_source(F, n, rs, a.points))

    cmd = ["cc", "-O2", "-march=native", "-I", "src", "-I", ".",
           "-I", "src/mp_real", "-o", binp, src,
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
        path = os.path.join(TOP, "src", "mp_real",
                            "%s_opt_%d.c" % (a.func, n))
        emit_production(F, n, r, path)
        print("emitted", path)


if __name__ == "__main__":
    main()
