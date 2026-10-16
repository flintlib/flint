/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "mpn_extras.h"
#include "mp_real.h"
#include "impl.h"

/* atan(t) for t in [0, 1) by greedy vectoring with limb arithmetic:
   the path of 32-bit machines, and of explicit r below
   MP_REAL_ATAN_BITWISE_VEC_MIN limbs on 64-bit ones, where r = 0 takes
   the per-size kernels and larger sizes the signed-digit vector
   reduction of bitwise_vec.c.  (The trigonometric analog of the
   _mp_real_log1p_bitwise_rs reduction): the vector
   (X, Y) = (1, t) is rotated towards the real axis by the factors
   (1 - i 2^-i), i = 1..r, applying a factor -- i.e.
   (X, Y) <- (X + Y 2^-i, Y - X 2^-i), two shifts and two
   add/subtracts on the old values -- whenever it keeps Y >= 0, that
   is whenever Y >= X 2^-i.  This is a greedy on the angle
   atan(Y/X) with steps atan(2^-i); concavity of atan gives
   atan(2^-(i-1)) < 2 atan(2^-i), so each factor is used at most
   once, one extra conditional step at i = r absorbs the truncation
   creep, and afterwards atan(Y/X) < atan(2^-r) < 2^-r.  The angle
   is scale invariant, so the growth of |Z| (bounded by
   prod sqrt(1 + 4^-i) < 1.1644, keeping X in [1, 2) with its units
   limb fixed at 1) needs no compensation.

   Decisions run one FLINT_BITS window at a time on a single-limb
   model as in log: h is Y's window limb and
   lt = MPN_RIGHT_SHIFT_LOW(X[n], X[n-1], i mod FLINT_BITS) the
   window limb of X >> i, computable from X's top two limbs alone;
   accepted factors update Y immediately, so the model resyncs
   exactly and the ambiguity band degenerates to h == lt.  As in
   log, the residual obeys only Y < X 2^-i < 2^(1-i), so a stray bit
   above the window forces a certain accept.

   The used indices are recorded and the tabulated angles
   A_i = atan(2^-i) -- the same table the rotation reduction of
   sin/cos uses (there against x/2) -- are simply summed afterwards;
   the total is below sum_{i>=1} A_i ~ 0.8980 < 1, so an n-limb
   accumulator suffices and no rescaling is needed.  The residual
   needs ONE division t' = Y/X < 2^-r, meeting the atan series
   contract.

   Error accounting (output ulp 2^(-FLINT_BITS n)): at most r + 1
   table entries, each short by < 1 ulp; each accepted
   factor truncates X and Y by < 1 ulp each, perturbing the residual
   angle by < 2 ulp per accept; the halving, the division (within one
   ulp) and the series error (<= 15 ulp two-sided) are constants.  Total
   below the returned bound *err = 4 r + 64. */

#define ATANS_A(i) (tab + (i) * nc)

/* smallest n at which each reduction parameter becomes optimal:
   interleaved timings of every ladder value at 26 sizes from 10 to
   640 limbs and at 800, 1000, 1300 and 1700 limbs (x86-64), the
   optimum being flat to a few percent over one or two ladder steps
   beyond ~100 limbs; the earlier tables from
   src/mp_real/tune/tune-bitwise-r.c, whose binary searches on single
   timings drifted with machine noise, chose r one to three steps too
   large from ~40 limbs on (5-25% slower) */
static const int _mp_real_atan_bitwise_rs_r_tab[] =
    {32, 64, 128, 192, 256, 320};
static const short _mp_real_atan_bitwise_rs_n_tab[] =
    {1, 18, 100, 270, 450, 620};

#if FLINT_BITS == 64
/* the same for _mp_real_atan_bitwise_vec (its own sweep, x86-64 AVX2:
   every ladder value at 27 sizes from 4 to 600 limbs, the minimum of
   three interleaved runs, the optimum flat to a few percent over one or
   two steps) */
static const int _mp_real_atan_bitwise_vec_r_tab[] =
    {32, 48, 64, 96, 192, 320, 384, 512};
static const short _mp_real_atan_bitwise_vec_n_tab[] =
    {1, 18, 30, 44, 88, 200, 240, 380};
#endif

int
_mp_real_atan_bitwise_rs_default_r(slong n)
{
    slong j;

#if FLINT_BITS == 64
    /* the reduction parameters of the per-size kernels that r = 0
       selects at these sizes (see impl.h) */
    static const int opt_r[] = MP_REAL_ATAN_OPT_R;
    if (n <= MP_REAL_ATAN_OPT_R_MAX)
        return opt_r[n];

    if (n >= MP_REAL_ATAN_BITWISE_VEC_MIN)
    {
        for (j = 0; j + 1 < (slong) (sizeof(_mp_real_atan_bitwise_vec_n_tab)
                / sizeof(short))
                && n >= _mp_real_atan_bitwise_vec_n_tab[j + 1]; j++)
            ;
        return _mp_real_atan_bitwise_vec_r_tab[j];
    }
#endif

    for (j = 0; j + 1 < (slong) (sizeof(_mp_real_atan_bitwise_rs_n_tab) / sizeof(short))
            && n >= _mp_real_atan_bitwise_rs_n_tab[j + 1]; j++)
        ;
    return _mp_real_atan_bitwise_rs_r_tab[j];
}

FLINT_FORCE_INLINE void
_mp_real_atan_bitwise_rs_noerr(nn_ptr res, nn_srcptr x, slong n, int r)
{
    slong i, c, wn, nc, num, j;
    nn_ptr X, Y, vx, vy, acc, t, y;
    slong * used;
    TMP_INIT;

    FLINT_ASSERT(n >= 1);
    FLINT_ASSERT(r == 0 || r >= 32);
    /* n = 1 is unsupported on 32-bit limbs: the clamp below would
       drop r under the r >= 32 contract of the residual series */
    FLINT_ASSERT(FLINT_BITS == 64 || n >= 2);

    if (r == 0)
        r = _mp_real_atan_bitwise_rs_default_r(n);

    /* the residual series contract is r >= 32 */
    r = FLINT_MAX(r, 32);
    r = FLINT_MIN((slong) r, FLINT_BITS * n - 16);

#if FLINT_BITS == 64
    if (n >= MP_REAL_ATAN_BITWISE_VEC_MIN)
    {
        _mp_real_atan_bitwise_vec(res, x, n, r);
        return;
    }
#endif

    wn = n + 1;

    nn_srcptr tab;
    tab = _mp_real_atans_tab(n, r, &nc);

    TMP_START;
    X = TMP_ALLOC((wn + n + wn + wn + wn + n + n) * sizeof(ulong));
    Y = X + wn;
    vx = Y + n;
    vy = vx + wn;
    acc = vy + wn;
    t = acc + wn;
    used = TMP_ALLOC(MP_REAL_BITWISE_REDUCE_USED_ALLOC(r) * sizeof(slong));

    /* X = 1, Y = t */
    flint_mpn_zero(X, n);
    X[n] = 1;
    flint_mpn_copyi(Y, x, n);
    num = 0;

/* materialize v = trunc(X >> ii) into vx (the leading unit of X
   lands inside the fraction for b == 0 boundaries) and
   trunc(Y >> ii) into vy, both as n-limb fractions */
#define MAKE_V(ii) \
    do { \
        slong q_ = (ii) / FLINT_BITS, b_ = (ii) - q_ * FLINT_BITS; \
        flint_mpn_zero(vx, n); \
        flint_mpn_zero(vy, n); \
        if (b_ != 0) \
        { \
            mpn_rshift(vx, X + q_, wn - q_, (int) b_); \
            if (n - q_ > 0) \
                mpn_rshift(vy, Y + q_, n - q_, (int) b_); \
        } \
        else \
        { \
            flint_mpn_copyi(vx, X + q_, wn - q_); \
            if (n - q_ > 0) \
                flint_mpn_copyi(vy, Y + q_, n - q_); \
        } \
    } while (0)

/* accept factor ii (v's already made): X += vy, Y -= vx */
#define APPLY(ii) \
    do { \
        mpn_add(X, X, wn, vy, n); \
        mpn_sub_n(Y, Y, vx, n); \
        used[num++] = (ii); \
    } while (0)

/* exact compare-accept: Y >= trunc(X >> ii)? */
#define EXACT_STEP(ii) \
    do { \
        MAKE_V(ii); \
        if (mpn_cmp(Y, vx, n) >= 0) \
            APPLY(ii); \
    } while (0)

    for (c = 0; FLINT_BITS * c <= r; c++)
    {
        slong i0 = (c == 0) ? 1 : FLINT_BITS * c;
        slong i1 = FLINT_MIN((slong) r, FLINT_BITS * (c + 1) - 1);
        ulong h;

        /* i0 <= r always: for c = 0 because r >= 32, and for
           c >= 1 by the loop condition */
        EXACT_STEP(i0);

        h = Y[n - 1 - c];

        for (i = i0 + 1; i <= i1; i++)
        {
            slong b = i - FLINT_BITS * c;
            ulong lt;

            /* stray bit above the window: certain accept, since
               then Y >= 2^(1 - i) b-shifted > X 2^-i */
            if (c > 0 && Y[n - c] != 0)
            {
                MAKE_V(i);
                APPLY(i);
                h = Y[n - 1 - c];
                continue;
            }

            lt = MPN_RIGHT_SHIFT_LOW(X[n], X[n - 1], (int) b);

            if (h > lt)
            {
                MAKE_V(i);
                APPLY(i);
                h = Y[n - 1 - c];
            }
            else if (h == lt)
            {
                EXACT_STEP(i);
                h = Y[n - 1 - c];
            }
            /* else certain reject */
        }
    }

    /* extra steps at i = r absorb the truncation creep, restoring
       Y < trunc(X >> r) <= X 2^-r, i.e. the residual angle below
       atan(2^-r): the series contract.  X grows with each accepted
       factor, so the comparison is recomputed each time. */
    for (;;)
    {
        slong num_before = num;

        EXACT_STEP(r);
        if (num == num_before)
            break;
    }

#undef MAKE_V
#undef APPLY
#undef EXACT_STEP

    /* sum the table entries: the total is below 0.898 < 1, so no
       carry limb is needed */
    flint_mpn_zero(acc, n);
    for (j = 0; j < num; j++)
        mpn_add_n(acc, acc, ATANS_A(used[j]) + (nc - n), n);

    /* t' = Y / X < 2^-r: one division, X in [1, 2) with the top
       limb X[n] = 1 */
    _mp_real_div_fraction(t, Y, X, n);

    /* res = acc + atan(t') (the series does not support aliased
       output) */
    y = t + n;
    _mp_real_atan_rs(y, NULL, t, n);
    mpn_add_n(res, acc, y, n);

    TMP_END;
}

void
_mp_real_atan_bitwise_rs(nn_ptr res, ulong * err, nn_srcptr x, slong n, int r)
{
#if FLINT_BITS == 64
    /* r = 0 leaves the reduction to us: at the sizes of the per-size
       kernels (small.h, the dense tables) those serve the call (any size
       they decline takes the bitwise path below) */
    if (r == 0 && n <= MP_REAL_ATAN_OPT_R_MAX)
    {
        ulong e;
        if (_mp_real_atan_opt(res, &e, x, n))
        {
            if (err != NULL)
                *err = e;
            return;
        }
    }
#endif

    _mp_real_atan_bitwise_rs_noerr(res, x, n, r);
    if (err != NULL)
        *err = 4 * (r ? (ulong) r : (ulong) _mp_real_atan_bitwise_rs_default_r(n)) + 64;
}
