/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Two products with a common operand, z1 = a x and z2 = b x:

      sep     two flint_mpn_mul calls
      shared  _fft_small_mpn_mul_shared (two-prime up to
              FLINT_MPN_MUL_NP2_MAX_AN limbs, wide beyond)
      cplx    flint_mpn_mul_complex_fft_small on (x + 0 i)(a + b i),
              which computes the same pair (by the shared code)
      tuned   flint_mpn_mul_shared, which chooses between shared and
              sep (the column "rule" shows its choice)

    printing the wall times and the ratios to sep: the best of five
    rounds, each timing the four variants in turn (TIMEIT), so that
    drift in the machine's speed reaches all of them alike. Every
    variant is checked against sep once per size. Each shape ends, for
    the two-prime and the wide sizes, with the geometric means of
    shared/sep and cplx/sep over the sizes where the rule chooses the
    transformed path, and the geometric mean and the maximum of
    tuned/sep over all sizes.

    Usage: p-mpn_mul_shared [threads [shape [nmin [nmax]]]]

    shape (xn = n throughout):
      equal  an = bn = n
      q      an = 1.17 n, bn = 0.86 n  (T1 Q2, Q1 Q2 in the merges of
                                        binary splitting for pi)
      r      an = 3.85 n, bn = 1.42 n  (R1 T2, R1 R2 there)
      short  an = n, bn = n / 8
      all    each of the above (default)

    n runs over powers of two times 1, 1.25, 1.5, 1.75 from nmin
    (default 128) to nmax (default 2^20; the r shape then takes about
    350 MB, a full run of all shapes some tens of minutes).
*/

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "mpn_extras.h"

#if FLINT_HAVE_FFT_SMALL

#include "fft_small.h"
#include "profiler.h"

static void
_mul(nn_ptr z, nn_srcptr a, slong an, nn_srcptr b, slong bn)
{
    if (an >= bn)
        flint_mpn_mul(z, a, an, b, bn);
    else
        flint_mpn_mul(z, b, bn, a, an);
}

#define ROUNDS 5
#define NV 4

/* one TIMEIT measurement (wall time per call) of variant v */
static double
_time_once(int v, nn_ptr z1, slong * l1, nn_ptr z2, slong * l2,
           nn_srcptr a, slong an, nn_srcptr b, slong bn,
           nn_srcptr x, slong xn)
{
    double FLINT_SET_BUT_UNUSED(tcpu), twall;

    TIMEIT_START;
    if (v == 0)
    {
        _mul(z1, a, an, x, xn);
        _mul(z2, b, bn, x, xn);
    }
    else if (v == 1)
    {
        _fft_small_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);
    }
    else if (v == 2)
    {
        if (!flint_mpn_mul_complex_fft_small(z1, l1, z2, l2,
                    x, xn, 0, x, 0, 0, a, an, 0, b, bn, 0))
            flint_abort();
    }
    else
    {
        flint_mpn_mul_shared(z1, a, an, z2, b, bn, x, xn);
    }
    TIMEIT_STOP_VALUES(tcpu, twall);

    return twall;
}

/* the outputs of variant v against the references r1, r2 */
static int
_check(int v, nn_ptr z1, slong l1, nn_ptr z2, slong l2, nn_srcptr r1,
       slong n1, nn_srcptr r2, slong n2)
{
    /* the complex product returns signed lengths of nonnegative parts
       here: zero-extend the magnitudes */
    if (v == 2)
    {
        if (l1 < 0 || l2 < 0 || l1 > n1 || l2 > n2)
            return 0;
        flint_mpn_zero(z1 + l1, n1 - l1);
        flint_mpn_zero(z2 + l2, n2 - l2);
    }

    return mpn_cmp(z1, r1, n1) == 0 && mpn_cmp(z2, r2, n2) == 0;
}

static void
_run_shape(const char * name, double fa, double fb, slong nmin, slong nmax,
           nn_srcptr data_a, nn_srcptr data_b, nn_srcptr data_x,
           nn_ptr z1, nn_ptr z2, nn_ptr r1, nn_ptr r2)
{
    static const double steps[4] = { 1.0, 1.25, 1.5, 1.75 };
    /* per geometry: log sums of shared, cplx over the rule's sizes, of
       tuned over all sizes; counts; the worst tuned/sep */
    double lsum[2][3] = { { 0, 0, 0 }, { 0, 0, 0 } }, worst[2] = { 0, 0 };
    slong cnt[2][2] = { { 0, 0 }, { 0, 0 } };
    slong p, j, g;

    flint_printf("shape %s (an = %.2f n, bn = %.3f n, xn = n), %wd thread(s)\n",
                 name, fa, fb, flint_get_num_threads());
    flint_printf("%9s %9s %9s | %9s %9s %9s %9s | %7s %7s %7s | %-4s %s\n",
                 "an", "bn", "xn", "sep", "shared", "cplx", "tuned",
                 "shr/sep", "cpx/sep", "tun/sep", "rule", "geometry");

    for (p = 1; p <= nmax; p *= 2)
    {
        for (j = 0; j < 4; j++)
        {
            slong n = (slong) (p * steps[j]);
            slong an = FLINT_MAX(1, (slong) (fa * n));
            slong bn = FLINT_MAX(1, (slong) (fb * n));
            slong xn = n;
            double t[NV];
            slong l1 = 0, l2 = 0;
            int good = 1, rule, k, v;

            if (n < nmin || n > nmax)
                continue;

            /* the rounds; sep writes the references r1, r2, the others
               z1, z2, checked after their last round */
            for (v = 0; v < NV; v++)
                t[v] = 1e300;
            for (k = 0; k < ROUNDS; k++)
            {
                for (v = 0; v < NV; v++)
                {
                    t[v] = FLINT_MIN(t[v], _time_once(v,
                            v ? z1 : r1, &l1, v ? z2 : r2, &l2,
                            data_a, an, data_b, bn, data_x, xn));
                    if (v != 0 && k == ROUNDS - 1)
                        good = good && _check(v, z1, l1, z2, l2,
                                              r1, an + xn, r2, bn + xn);
                }
            }

            rule = _flint_mpn_mul_shared_use_fft(an, bn, xn);
            g = (FLINT_MAX(FLINT_MAX(an, bn), xn) > FLINT_MPN_MUL_NP2_MAX_AN);

            flint_printf("%9wd %9wd %9wd | %9.3g %9.3g %9.3g %9.3g | "
                         "%7.3f %7.3f %7.3f | %-4s %s%s\n",
                         an, bn, xn, t[0], t[1], t[2], t[3],
                         t[1] / t[0], t[2] / t[0], t[3] / t[0],
                         rule ? "fft" : "sep", g ? "wide" : "two-prime",
                         good ? "" : "   WRONG RESULT");

            if (rule)
            {
                lsum[g][0] += log(t[1] / t[0]);
                lsum[g][1] += log(t[2] / t[0]);
                cnt[g][0]++;
            }
            lsum[g][2] += log(t[3] / t[0]);
            cnt[g][1]++;
            worst[g] = FLINT_MAX(worst[g], t[3] / t[0]);
        }
    }

    for (g = 0; g < 2; g++)
    {
        if (cnt[g][1] == 0)
            continue;
        flint_printf("%s sizes:", g ? "wide" : "two-prime");
        if (cnt[g][0] != 0)
            flint_printf(" where the rule takes fft, shared/sep %.3f, "
                "cplx/sep %.3f;", exp(lsum[g][0] / cnt[g][0]),
                exp(lsum[g][1] / cnt[g][0]));
        flint_printf(" tuned/sep geometric mean %.3f, max %.3f\n",
            exp(lsum[g][2] / cnt[g][1]), worst[g]);
    }

    flint_printf("\n");
}

int
main(int argc, char ** argv)
{
    static const char * names[4] = { "equal", "q", "r", "short" };
    static const double fa[4] = { 1.0, 1.17, 3.85, 1.0 };
    static const double fb[4] = { 1.0, 0.86, 1.42, 0.125 };
    slong nmin = 128, nmax = WORD(1) << 20, len, k;
    double fmax = 1.0;
    const char * shape = "all";
    nn_ptr a, b, x, z1, z2, r1, r2;
    flint_rand_t state;

    if (argc > 1)
        flint_set_num_threads(atoi(argv[1]));
    if (argc > 2)
        shape = argv[2];
    if (argc > 3)
        nmin = atol(argv[3]);
    if (argc > 4)
        nmax = atol(argv[4]);

    /* room for the largest operand of the shapes run */
    for (k = 0; k < 4; k++)
        if (!strcmp(shape, "all") || !strcmp(shape, names[k]))
            fmax = FLINT_MAX(fmax, FLINT_MAX(fa[k], fb[k]));
    len = (slong) (fmax * nmax) + 2;

    flint_rand_init(state);
    a = flint_malloc(len * sizeof(ulong));
    b = flint_malloc(len * sizeof(ulong));
    x = flint_malloc(len * sizeof(ulong));
    z1 = flint_malloc(2 * len * sizeof(ulong));
    z2 = flint_malloc(2 * len * sizeof(ulong));
    r1 = flint_malloc(2 * len * sizeof(ulong));
    r2 = flint_malloc(2 * len * sizeof(ulong));
    flint_mpn_urandomb(a, state, len * FLINT_BITS);
    flint_mpn_urandomb(b, state, len * FLINT_BITS);
    flint_mpn_urandomb(x, state, len * FLINT_BITS);

    for (k = 0; k < 4; k++)
        if (!strcmp(shape, "all") || !strcmp(shape, names[k]))
            _run_shape(names[k], fa[k], fb[k], nmin, nmax, a, b, x,
                       z1, z2, r1, r2);

    flint_free(a);
    flint_free(b);
    flint_free(x);
    flint_free(z1);
    flint_free(z2);
    flint_free(r1);
    flint_free(r2);
    flint_rand_clear(state);
    flint_cleanup_master();
    return 0;
}

#else

int main(void)
{
    return 0;
}

#endif
