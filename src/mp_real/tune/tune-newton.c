/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Timing of _mp_real_neglog_newton and _mp_real_atan_newton against the
   forward functions they are built on and against the bitwise
   inverses, with a scan over the number of terms N and the row
   width m.

   usage: tune-newton n [reps] [Nmin Nmax] [scan_m]

   Prints, for the size n (limbs), the time of the forward functions,
   the bitwise inverses, and of each Newton variant with its worst
   error in ulps against arb. */

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "arb.h"
#include "mp_real.h"
#include "mp_real/impl.h"

static double
_ulps(nn_srcptr y, slong n, const arb_t ref)
{
    arb_t v, d;
    fmpz_t f;
    mag_t m;
    double u;

    arb_init(v); arb_init(d); fmpz_init(f); mag_init(m);
    fmpz_set_ui_array(f, y, n);
    arb_set_fmpz(v, f);
    arb_mul_2exp_si(v, v, -FLINT_BITS * n);
    arb_sub(d, v, ref, FLINT_BITS * n + 64);
    arb_get_mag(m, d);
    mag_mul_2exp_si(m, m, FLINT_BITS * n);
    u = mag_get_d(m);
    arb_clear(v); arb_clear(d); fmpz_clear(f); mag_clear(m);
    return u;
}

#include <sys/time.h>
static double
wall(void)
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + 1e-6 * tv.tv_usec;
}

/* best-of-reps time of f applied to the inputs in xs */
#define TIME(expr, reps, out) do { double _b = 1e300; int _r; \
    for (_r = 0; _r < (reps); _r++) { double _t0 = wall(); expr; \
    _b = FLINT_MIN(_b, wall() - _t0); } out = _b; } while (0)

int
main(int argc, char ** argv)
{
    flint_rand_t st;
    slong n, N, Nmin, Nmax;
    int reps;
    nn_ptr x, xl, y, y2;
    arb_t xa, xla, rl, ra;
    fmpz_t f;
    double tf, tb, t0, t1;

    flint_rand_init(st);
    n = (argc > 1) ? atol(argv[1]) : 1000;
    reps = (argc > 2) ? atoi(argv[2]) : 3;
    Nmin = (argc > 3) ? atol(argv[3]) : 0;
    Nmax = (argc > 4) ? atol(argv[4]) : 0;

    x = flint_malloc(n * sizeof(ulong));
    xl = flint_malloc(n * sizeof(ulong));
    y = flint_malloc((n + 2) * sizeof(ulong));
    y2 = flint_malloc((n + 2) * sizeof(ulong));
    arb_init(xa); arb_init(xla); arb_init(rl); arb_init(ra); fmpz_init(f);

    flint_mpn_urandomb(x, st, FLINT_BITS * n);
    x[n - 1] &= ~(UWORD(1) << (FLINT_BITS - 1));       /* atan input in [0, 1/2) */
    flint_mpn_copyi(xl, x, n);
    xl[n - 1] |= UWORD(1) << (FLINT_BITS - 1);         /* log input in [1/2, 1) */

    fmpz_set_ui_array(f, x, n); arb_set_fmpz(xa, f); arb_mul_2exp_si(xa, xa, -FLINT_BITS * n);
    fmpz_set_ui_array(f, xl, n); arb_set_fmpz(xla, f); arb_mul_2exp_si(xla, xla, -FLINT_BITS * n);
    arb_log(rl, xla, FLINT_BITS * n + 64); arb_neg(rl, rl);
    arb_atan(ra, xa, FLINT_BITS * n + 64);

    /* warm the caches (tables, constants) */
    _mp_real_neglog_newton(y, NULL, xl, n);
    _mp_real_atan_newton(y, NULL, x, n);

    flint_printf("n = %wd limbs (%wd bits)\n", n, FLINT_BITS * n);

    /* -log */
    TIME(_mp_real_exp_diophantine(y, NULL, x, n + 1), reps, tf);
    flint_printf("  exp_diophantine(n+1)   %10.3f ms\n", tf * 1e3);
    if (n <= 40000)
    {
        TIME(_mp_real_log1p_bitwise_rs(y, NULL, x, n, 0), reps, tb);
        flint_printf("  log1p_bitwise          %10.3f ms  (%.2fx exp_dioph)\n", tb * 1e3, tb / tf);
    }
    TIME(_mp_real_neglog_newton(y, NULL, xl, n), reps, t0);
    flint_printf("  neglog_newton default  %10.3f ms  (%.3fx exp_dioph) %.2f ulps\n",
        t0 * 1e3, t0 / tf, _ulps(y, n, rl));
    for (N = Nmin; N <= Nmax; N++)
    {
        TIME(_mp_real_neglog_newton_tune(y, NULL, xl, n, 0, N), reps, t1);
        flint_printf("    N = %2wd              %10.3f ms  (%.3fx exp_dioph) %.2f ulps\n",
            N, t1 * 1e3, t1 / tf, _ulps(y, n, rl));
    }

    /* atan */
    TIME(_mp_real_sin_cos_diophantine(y, y2, NULL, x, n + 1), reps, tf);
    flint_printf("  sin_cos_diophantine    %10.3f ms\n", tf * 1e3);
    if (n <= 40000)
    {
        TIME(_mp_real_atan_bitwise_rs(y, NULL, x, n, 0), reps, tb);
        flint_printf("  atan_bitwise           %10.3f ms  (%.2fx sin_cos_dioph)\n", tb * 1e3, tb / tf);
    }
    TIME(_mp_real_atan_newton(y, NULL, x, n), reps, t0);
    flint_printf("  atan_newton default    %10.3f ms  (%.3fx sin_cos_dioph) %.2f ulps\n",
        t0 * 1e3, t0 / tf, _ulps(y, n, ra));
    for (N = Nmin; N <= Nmax; N++)
    {
        TIME(_mp_real_atan_newton_tune(y, NULL, x, n, 0, N), reps, t1);
        flint_printf("    N = %2wd              %10.3f ms  (%.3fx sin_cos_dioph) %.2f ulps\n",
            N, t1 * 1e3, t1 / tf, _ulps(y, n, ra));
    }

    flint_free(x); flint_free(xl); flint_free(y); flint_free(y2);
    arb_clear(xa); arb_clear(xla); arb_clear(rl); arb_clear(ra); fmpz_clear(f);
    flint_rand_clear(st);
    flint_cleanup();
    return 0;
}
