/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "profiler.h"
#include "longlong.h"
#include "ulong_extras.h"
#include "mpn_extras.h"

/* fallback for older FLINT versions (mirrors longlong.h) */
#if !defined(sub_ddddmmmmssss)
# define sub_ddddmmmmssss(s3, s2, s1, s0, a3, a2, a1, a0, b3, b2, b1, b0) \
  do { \
    ulong __t3, __t4; \
    sub_dddmmmsss(__t3, s1, s0, (ulong) 0, a1, a0, (ulong) 0, b1, b0); \
    sub_ddmmss(__t4, s2, (ulong) 0, a2, (ulong) 0, b2); \
    sub_ddmmss(s3, s2, (a3) - (b3), s2, -__t4, -__t3); \
  } while (0)
#endif
#include "gr.h"
#include "gr_dft.h"

/* Microbenchmarks of the fixed-point arithmetic primitives, for
   machine tuning of the implementation choices in nfixed.c:

   A. magnitude addition, 1-4 limbs: plain-C carry chains (as used by
      the sized methods) against the inline-assembly chains of
      longlong.h (as used by nfloat);
   B. signed (sign-magnitude) subtraction: mpn_cmp followed by an
      ordered mpn_sub_n (as in flint_mpn_signed_sub_n / nfloat)
      against an unconditional mpn_sub_n followed by conditional
      two's-complement negation, the latter both with a hand-written
      loop (as currently in nfixed.c) and with mpn_neg; and the same
      comparison for 1-4 limb specialized code;
   C. the cost of the saturation checks: saturating against
      non-saturating signed addition on random non-saturating inputs,
      for the 1-4 limb specialized code and generic sizes;
   D. classical against Karatsuba complex multiplication, for tuning
      GR_DFT_NFIXED_KARATSUBA_CUTOFF.

   Operands are random with magnitudes below 1/4, so the saturation
   and sign branches are exercised but never taken pathologically;
   per-operation times include the shared loop and store overhead. */

#define VEC 1024
#define MAXN 128

static ulong buf_a[VEC][MAXN + 1];
static ulong buf_b[VEC][MAXN + 1];
static ulong buf_r[VEC][MAXN + 1];

static void
fill(flint_rand_t state, slong n, int signs)
{
    slong i, j;
    for (i = 0; i < VEC; i++)
    {
        buf_a[i][0] = signs ? n_randint(state, 2) : 0;
        buf_b[i][0] = signs ? n_randint(state, 2) : 0;
        for (j = 1; j <= n; j++)
        {
            buf_a[i][j] = n_randlimb(state);
            buf_b[i][j] = n_randlimb(state);
        }
        buf_a[i][n] >>= 2;      /* magnitudes below 1/4 */
        buf_b[i][n] >>= 2;
    }
}

typedef void (*op_fn)(nn_ptr, nn_srcptr, nn_srcptr, slong);

static double
time_op(op_fn f, slong n)
{
    timeit_t tm;
    slong i, r, reps = 4;
    double t;

    for (i = 0; i < VEC; i++)
        f(buf_r[i], buf_a[i], buf_b[i], n);

    for (;;)
    {
        timeit_start(tm);
        for (r = 0; r < reps; r++)
            for (i = 0; i < VEC; i++)
                f(buf_r[i], buf_a[i], buf_b[i], n);
        timeit_stop(tm);
        if (tm->wall >= 100)
        {
            t = 1e6 * tm->wall / ((double) reps * VEC);
            break;
        }
        reps *= 4;
    }
    return t;
}

/* A. magnitude addition (no sign limb; a, b, r point at magnitudes) */

static void addC_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    a++; b++; r++;
    if (n == 1)
    {
        r[0] = a[0] + b[0];
    }
    else if (n == 2)
    {
        ulong c;
        r[0] = a[0] + b[0]; c = r[0] < a[0];
        r[1] = a[1] + b[1] + c;
    }
    else if (n == 3)
    {
        ulong c;
        r[0] = a[0] + b[0]; c = r[0] < a[0];
        r[1] = a[1] + b[1] + c; c = (r[1] < a[1]) | (c & (r[1] == a[1]));
        r[2] = a[2] + b[2] + c;
    }
    else
    {
        ulong c;
        r[0] = a[0] + b[0]; c = r[0] < a[0];
        r[1] = a[1] + b[1] + c; c = (r[1] < a[1]) | (c & (r[1] == a[1]));
        r[2] = a[2] + b[2] + c; c = (r[2] < a[2]) | (c & (r[2] == a[2]));
        r[3] = a[3] + b[3] + c;
    }
}

static void addA_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    a++; b++; r++;
    if (n == 1)
    {
        r[0] = a[0] + b[0];
    }
    else if (n == 2)
    {
        add_ssaaaa(r[1], r[0], a[1], a[0], b[1], b[0]);
    }
    else if (n == 3)
    {
        add_sssaaaaaa(r[2], r[1], r[0], a[2], a[1], a[0], b[2], b[1], b[0]);
    }
    else
    {
        add_ssssaaaaaaaa(r[3], r[2], r[1], r[0],
                a[3], a[2], a[1], a[0], b[3], b[2], b[1], b[0]);
    }
}

/* B. signed subtraction of magnitudes: r = |a| - |b| with returned
   sign; here the sign is stored in r[-1]... we use the leading slot */

static void ssub_cmp_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (mpn_cmp(a + 1, b + 1, n) >= 0)
    {
        mpn_sub_n(r + 1, a + 1, b + 1, n);
        r[0] = 0;
    }
    else
    {
        mpn_sub_n(r + 1, b + 1, a + 1, n);
        r[0] = 1;
    }
}

static void ssub_negloop_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (mpn_sub_n(r + 1, a + 1, b + 1, n))
    {
        slong i;
        ulong cy = 1;
        for (i = 1; i <= n; i++)
        {
            r[i] = ~r[i] + cy;
            cy = cy & (r[i] == 0);
        }
        r[0] = 1;
    }
    else
        r[0] = 0;
}

static void ssub_mpnneg_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (mpn_sub_n(r + 1, a + 1, b + 1, n))
    {
        mpn_neg(r + 1, r + 1, n);
        r[0] = 1;
    }
    else
        r[0] = 0;
}

/* sized: compare-first */
static void ssub_cmp_sized_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    int swap;
    nn_srcptr u, v;
    slong i;

    swap = 0;
    for (i = n; i >= 1; i--)
    {
        if (a[i] != b[i])
        {
            swap = a[i] < b[i];
            break;
        }
    }
    u = swap ? b : a;
    v = swap ? a : b;
    r[0] = swap;

    if (n == 1)
    {
        r[1] = u[1] - v[1];
    }
    else if (n == 2)
    {
        sub_ddmmss(r[2], r[1], u[2], u[1], v[2], v[1]);
    }
    else if (n == 3)
    {
        sub_dddmmmsss(r[3], r[2], r[1], u[3], u[2], u[1], v[3], v[2], v[1]);
    }
    else
    {
        sub_ddddmmmmssss(r[4], r[3], r[2], r[1],
                u[4], u[3], u[2], u[1], v[4], v[3], v[2], v[1]);
    }
}

/* sized: subtract, then conditional negation (plain C, mirroring the
   MAG_SSUB macros of nfixed.c) */
static void ssub_neg_sized_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    ulong bw;

    if (n == 1)
    {
        r[1] = a[1] - b[1];
        bw = a[1] < b[1];
        if (bw)
            r[1] = -r[1];
    }
    else if (n == 2)
    {
        ulong d0, d1, b0;
        d0 = a[1] - b[1]; b0 = a[1] < b[1];
        d1 = a[2] - b[2]; bw = a[2] < b[2];
        bw |= (b0 & (d1 == 0)); d1 -= b0;
        r[1] = d0; r[2] = d1;
        if (bw)
        {
            ulong cy;
            r[1] = -r[1]; cy = (r[1] == 0);
            r[2] = ~r[2] + cy;
        }
    }
    else if (n == 3)
    {
        sub_dddmmmsss(r[3], r[2], r[1], a[3], a[2], a[1], b[3], b[2], b[1]);
        bw = (a[3] < b[3]) || (a[3] == b[3] &&
                (a[2] < b[2] || (a[2] == b[2] && a[1] < b[1])));
        if (bw)
        {
            ulong cy;
            r[1] = -r[1]; cy = (r[1] == 0);
            r[2] = ~r[2] + cy; cy = cy & (r[2] == 0);
            r[3] = ~r[3] + cy;
        }
    }
    else
    {
        sub_ddddmmmmssss(r[4], r[3], r[2], r[1],
                a[4], a[3], a[2], a[1], b[4], b[3], b[2], b[1]);
        bw = r[4] >> (FLINT_BITS - 1);
        if (bw)
        {
            ulong cy;
            r[1] = -r[1]; cy = (r[1] == 0);
            r[2] = ~r[2] + cy; cy = cy & (r[2] == 0);
            r[3] = ~r[3] + cy; cy = cy & (r[3] == 0);
            r[4] = ~r[4] + cy;
        }
    }
    r[0] = bw;
}

/* C. full signed addition, saturating vs not (generic via mpn) */

static void sadd_sat_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (a[0] == b[0])
    {
        r[0] = a[0];
        if (mpn_add_n(r + 1, a + 1, b + 1, n))
        {
            slong i;
            for (i = 1; i <= n; i++)
                r[i] = ~UWORD(0);
        }
    }
    else
    {
        r[0] = a[0] ^ (ulong) (mpn_cmp(a + 1, b + 1, n) < 0);
        if (mpn_cmp(a + 1, b + 1, n) >= 0)
            mpn_sub_n(r + 1, a + 1, b + 1, n);
        else
            mpn_sub_n(r + 1, b + 1, a + 1, n);
    }
}

static void sadd_nosat_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (a[0] == b[0])
    {
        r[0] = a[0];
        mpn_add_n(r + 1, a + 1, b + 1, n);
    }
    else
    {
        r[0] = a[0] ^ (ulong) (mpn_cmp(a + 1, b + 1, n) < 0);
        if (mpn_cmp(a + 1, b + 1, n) >= 0)
            mpn_sub_n(r + 1, a + 1, b + 1, n);
        else
            mpn_sub_n(r + 1, b + 1, a + 1, n);
    }
}

/* sized signed addition with/without saturation (2 limbs shown as the
   representative case; dispatched below) */
static void sadd_sat_sized_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (a[0] == b[0])
    {
        ulong cy;
        r[0] = a[0];
        if (n == 1)
        {
            r[1] = a[1] + b[1]; cy = r[1] < a[1];
        }
        else if (n == 2)
        {
            ulong t;
            add_ssaaaa(t, r[1], a[2], a[1], b[2], b[1]);
            cy = t < a[2]; r[2] = t;
        }
        else if (n == 3)
        {
            ulong t;
            add_sssaaaaaa(t, r[2], r[1], a[3], a[2], a[1], b[3], b[2], b[1]);
            cy = t < a[3]; r[3] = t;
        }
        else
        {
            ulong t;
            add_ssssaaaaaaaa(t, r[3], r[2], r[1],
                    a[4], a[3], a[2], a[1], b[4], b[3], b[2], b[1]);
            cy = t < a[4]; r[4] = t;
        }
        if (cy)
        {
            slong i;
            for (i = 1; i <= n; i++)
                r[i] = ~UWORD(0);
        }
    }
    else
        ssub_neg_sized_op(r, a, b, n), r[0] ^= a[0];
}

static void sadd_nosat_sized_op(nn_ptr r, nn_srcptr a, nn_srcptr b, slong n)
{
    if (a[0] == b[0])
    {
        r[0] = a[0];
        if (n == 1)
            r[1] = a[1] + b[1];
        else if (n == 2)
            add_ssaaaa(r[2], r[1], a[2], a[1], b[2], b[1]);
        else if (n == 3)
            add_sssaaaaaa(r[3], r[2], r[1], a[3], a[2], a[1], b[3], b[2], b[1]);
        else
            add_ssssaaaaaaaa(r[4], r[3], r[2], r[1],
                    a[4], a[3], a[2], a[1], b[4], b[3], b[2], b[1]);
    }
    else
        ssub_neg_sized_op(r, a, b, n), r[0] ^= a[0];
}

int
main(void)
{
    flint_rand_t state;
    slong n, i;

    flint_rand_init(state);

    flint_printf("nfixed primitive microbenchmarks (ns per operation,\n");
    flint_printf("including loop and store overhead)\n\n");

    flint_printf("A. magnitude addition, plain-C carry chains vs "
            "longlong chains\n");
    flint_printf("%8s %10s %10s\n", "limbs", "plain C", "longlong");
    for (n = 1; n <= 4; n++)
    {
        fill(state, n, 0);
        flint_printf("%8wd %10.2f %10.2f\n", n,
                time_op(addC_op, n), time_op(addA_op, n));
    }

    flint_printf("\nB. signed subtraction of magnitudes\n");
    flint_printf("%8s %10s %10s %10s   (sized: cmp-first vs sub+neg)\n",
            "limbs", "cmp+sub", "sub+neg", "");
    for (n = 1; n <= 4; n++)
    {
        fill(state, n, 0);
        flint_printf("%8wd %10.2f %10.2f\n", n,
                time_op(ssub_cmp_sized_op, n),
                time_op(ssub_neg_sized_op, n));
    }
    flint_printf("%8s %10s %10s %10s   (generic mpn)\n",
            "limbs", "cmp+sub", "sub+loop", "sub+mpn_neg");
    for (n = 4; n <= MAXN; n *= 2)
    {
        fill(state, n, 0);
        flint_printf("%8wd %10.2f %10.2f %10.2f\n", n,
                time_op(ssub_cmp_op, n),
                time_op(ssub_negloop_op, n),
                time_op(ssub_mpnneg_op, n));
    }

    flint_printf("\nC. signed addition, saturating vs non-saturating\n");
    flint_printf("%8s %10s %10s   (sized)\n", "limbs", "sat", "no sat");
    for (n = 1; n <= 4; n++)
    {
        fill(state, n, 1);
        flint_printf("%8wd %10.2f %10.2f\n", n,
                time_op(sadd_sat_sized_op, n),
                time_op(sadd_nosat_sized_op, n));
    }
    flint_printf("%8s %10s %10s   (generic mpn)\n", "limbs", "sat", "no sat");
    for (n = 4; n <= MAXN; n *= 2)
    {
        fill(state, n, 1);
        flint_printf("%8wd %10.2f %10.2f\n", n,
                time_op(sadd_sat_op, n), time_op(sadd_nosat_op, n));
    }

    flint_printf("\nD. complex multiplication, schoolbook vs Karatsuba\n");
    flint_printf("%8s %12s %12s %8s\n", "limbs", "schoolbook", "karatsuba",
            "ratio");
    {
        slong nls[] = { 2, 4, 8, 12, 16, 20, 24, 32, 48, 64, 96, 128, 0 };
        slong t;

        for (t = 0; nls[t] != 0; t++)
        {
            gr_ctx_t cctx;
            gr_ptr x, y, r;
            timeit_t tm;
            slong reps, rr;
            double ts, tk;
            int status = GR_SUCCESS;

            n = nls[t];
            if (gr_dft_ctx_init_nfixed_complex(cctx, n) != GR_SUCCESS)
                continue;

            x = gr_heap_init(cctx);
            y = gr_heap_init(cctx);
            r = gr_heap_init(cctx);
            status |= gr_randtest(x, state, cctx);
            status |= gr_randtest(y, state, cctx);

            status |= _gr_dft_nfixed_cmul_schoolbook(r, x, y, cctx);
            for (reps = 4; ; reps *= 4)
            {
                timeit_start(tm);
                for (rr = 0; rr < reps; rr++)
                    status |= _gr_dft_nfixed_cmul_schoolbook(r, x, y, cctx);
                timeit_stop(tm);
                if (tm->wall >= 100)
                    break;
            }
            ts = 1e6 * tm->wall / reps;

            status |= _gr_dft_nfixed_cmul_karatsuba(r, x, y, cctx);
            for (reps = 4; ; reps *= 4)
            {
                timeit_start(tm);
                for (rr = 0; rr < reps; rr++)
                    status |= _gr_dft_nfixed_cmul_karatsuba(r, x, y, cctx);
                timeit_stop(tm);
                if (tm->wall >= 100)
                    break;
            }
            tk = 1e6 * tm->wall / reps;

            flint_printf("%8wd %12.1f %12.1f %8.2f  (status %d)\n",
                    n, ts, tk, ts / tk, status);

            gr_heap_clear(x, cctx);
            gr_heap_clear(y, cctx);
            gr_heap_clear(r, cctx);
            gr_ctx_clear(cctx);
        }
    }

    /* keep the results alive */
    {
        ulong chk = 0;
        for (i = 0; i < VEC; i++)
            chk ^= buf_r[i][1];
        if (chk == UWORD(0x123456789abcdef))
            flint_printf("(unlikely)\n");
    }

    flint_rand_clear(state);
    return 0;
}
