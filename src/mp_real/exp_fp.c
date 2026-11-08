/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "impl.h"

/*
    Exponentials of floating-point mantissas (mp_real.h), the kernels of
    the exponentials of nfloat, arf and arb.

    exp(x) for |x| < 2^(FLINT_BITS - 3) is reduced as |x| = q log 2 + t,
    t in [0, log 2), in fixed point with one integral and nk + 1 fraction
    limbs, nk the kernel limbs: q is a lower approximation of |x| / log 2
    from two limbs of 1/log 2 (exact or one low, then corrected; the
    reduction step is that of exp.c), and log 2 is the floor L of
    log 2 B^(nk + 1), so that the reduction errs by
    (q + 1)(log 2 - L) < 2^(FLINT_BITS - 2) B^-(nk + 1), a quarter ulp of
    B^-nk.  For x < 0, exp(x) = 2^-(q+1) exp(L - t'), t' = |x| - q log 2.
    Then exp(x) = 2^k exp(t) with the kernel on t at nk fraction limbs.
    Positive x < 1 go to the kernel directly, and small negative x to
    the alternating series of exp(-|x|).

    Beyond 2^(FLINT_BITS - 3) the quotient is a bignum: with
    N = nk + 1 + ceil((e + 3) / FLINT_BITS) fraction limbs and
    X = floor(|x| B^N), q is a lower approximation of |x| / log 2 (exact
    or one low) from the top limbs of X times the tabulated 1/log 2, as
    in the trigonometric reduction of sin_cos.c, and t = X - q L by one
    exact product, L the floor of log 2 B^N; the reduction errs by
    (q + 1)(log 2 - L) < 2^(e + 1) B^-N, below 2^-(FLINT_BITS nk + 66).

    expm1 subtracts 1 in the same fixed point, exp2(x) = 2^k exp(f log 2)
    with k = floor(x), f = x - k in [0, 1) (the reduction exact apart
    from the truncation of f at nk + 1 fraction limbs and the product
    f L, L the floor of log 2 B^(nk + 1)).
*/

slong
_mp_real_exp_fp(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb)
{
    slong nk = _mp_real_elem_limbs(pb);
    nn_ptr X;
    nn_srcptr L;
    ulong q;
    slong N, k;
    int trunc;
    TMP_INIT;

    FLINT_ASSERT(e <= MP_REAL_EXP_FP_MAX_EXP);

    TMP_START;
    X = TMP_ALLOC((nk + 3) * sizeof(ulong));

    if (e <= 0)
    {
        if (!sgnbit)
        {
            /* x in (0, 1): truncating x costs one ulp, times exp' < e */
            trunc = _mp_real_fp_get_fixed(X, nk, nk, d, n, e);
            _mp_real_exp_kernel(y, err, X, nk);
            *err += 3 * trunc;
            TMP_END;
            return 0;
        }

        /* the series of exp(-|x|) where it is faster (exp' <= 1) */
        trunc = _mp_real_fp_get_fixed(X, nk, nk, d, n, e);
        if (_mp_real_exp_neg_series(y, err, X, nk))
        {
            *err += trunc;
            TMP_END;
            return 0;
        }
    }

    /* X = |x| with one integral and N fraction limbs (truncated) */
    N = nk + 1;
    _mp_real_fp_get_fixed(X, N + 1, N, d, n, e);
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);

    /* q <= |x| / log 2, exact or one less, and t = X - q L resp.
       (q + 1) L - |x| from the two's complement (impl.h) */
    q = _mp_real_exp_quotient(X[N], X[N - 1]);
    if (sgnbit)
        mpn_neg(X, X, N + 1);
    k = _mp_real_exp_reduce_step(X, L, N, q, sgnbit);

    FLINT_ASSERT(X[N] == 0);

    /* t within two ulps of B^-nk (the reduction, the truncation of x and
       of the low fraction limb), times exp' < 2 */
    _mp_real_exp_kernel(y, err, X + 1, nk);
    *err += 4;
    FLINT_ASSERT(y[nk] <= 2);
    TMP_END;
    return k;
}

int
_mp_real_exp_fp_huge(nn_ptr y, ulong * err, nn_ptr q, slong * qn, nn_srcptr d, slong n, slong e, int sgnbit, slong pb)
{
    slong nk = _mp_real_elem_limbs(pb);
    nn_ptr X, pr, T;
    nn_srcptr L;
    slong m, N, xl, k, ql;
    TMP_INIT;

    FLINT_ASSERT(e > MP_REAL_EXP_FP_MAX_EXP);

    /* m integral limbs; N fraction limbs for log 2 and X, where
       (q + 1) < 2^(e + 1) times (log 2 - L) < B^-N stays below
       2^-(FLINT_BITS nk + 66) */
    m = (e + FLINT_BITS - 1) / FLINT_BITS;
    N = nk + 1 + (e + 3 + FLINT_BITS - 1) / FLINT_BITS;
    xl = m + N;
    k = m + 2;

    TMP_START;
    X = TMP_ALLOC((xl + 2 * m + 4 + 2 * m + N + 4) * sizeof(ulong));
    pr = X + xl;
    T = pr + 2 * m + 4;

    /* X = floor(|x| B^N), the error below one ulp of B^-N */
    _mp_real_fp_get_fixed(X, xl, N, d, n, e);
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);

    /* q = floor of a lower approximation of |x| / log 2 = |x| (1 + f):
       the top m + 1 limbs of X (m integral, one fraction) times the
       floor of f B^k, plus the same limbs shifted by k, the limbs k + 1
       up of the sum.  The truncations of X and f drop less than
       2/B + 2^(e - FLINT_BITS k) < 1 before the floor, so q is exact or
       one low. */
    flint_mpn_mul(pr, _mp_real_const_ptr(MP_REAL_CONST_ID_INV_LOG2_FRAC, k), k, X + N - 1, m + 1);
    pr[2 * m + 3] = mpn_add_n(pr + k, pr + k, X + N - 1, m + 1);
    flint_mpn_copyi(q, pr + k + 1, m + 1);
    ql = m + 1;
    while (ql > 1 && q[ql - 1] == 0)
        ql--;

    /* X = |x| - q log 2 in [0, 2 log 2): q L as ql + N limbs (its top
       limb zero when that exceeds xl) */
    flint_mpn_mul(T, L, N, q, ql);
    if (ql + N > xl)
    {
        FLINT_ASSERT(T[xl] == 0);
        mpn_sub_n(X, X, T, xl);
    }
    else
        mpn_sub(X, X, xl, T, ql + N);

    /* one correction when q was one low */
    if (!flint_mpn_zero_p(X + N, m) || mpn_cmp(X, L, N) >= 0)
    {
        mpn_sub(X, X, xl, L, N);
        q[ql] = mpn_add_1(q, q, ql, 1);
        ql += (q[ql] != 0);
    }
    FLINT_ASSERT(flint_mpn_zero_p(X + N, m));

    if (!sgnbit)
    {
        /* exp(x) = 2^q exp(t), t = X B^-N, the top nk limbs of X within
           two ulps of B^-nk (the reduction, the truncation of x and of
           the dropped limbs), times exp' < 2 */
        *qn = ql;
        _mp_real_exp_kernel(y, err, X + (N - nk), nk);
        *err += 4;
        TMP_END;
        return 0;
    }
    else
    {
        /* exp(x) = 2^-(q + 1) exp(L - t) */
        mpn_sub_n(X, L, X, N);
        q[ql] = mpn_add_1(q, q, ql, 1);
        *qn = ql + (q[ql] != 0);
        _mp_real_exp_kernel(y, err, X + (N - nk), nk);
        *err += 4;
        TMP_END;
        return 1;
    }
}

slong
_mp_real_expm1_fp(nn_ptr y, ulong * elo, ulong * ehi, nn_srcptr d, slong n, slong e, int sgnbit, slong pb)
{
    slong nk = _mp_real_elem_limbs(pb);
    ulong err;
    slong k;

    k = _mp_real_exp_fp(y, &err, d, n, e, sgnbit, pb);

    if (e <= 0)
    {
        if (!sgnbit)
        {
            /* exp(x) - 1 in [0, e - 1) */
            FLINT_ASSERT(k == 0);
            y[nk] -= 1;
            *elo = *ehi = err;
            return 0;
        }
        else
        {
            /* 1 - Y 2^k B^-nk, k in {0, -1, -2}: shift Y down by -k bits
               (one more unit for the dropped bits) */
            if (k != 0)
            {
                int dropped;
                FLINT_ASSERT(k >= -2);
                dropped = (y[0] & ((UWORD(1) << (-k)) - 1)) != 0;
                mpn_rshift(y, y, nk + 1, (unsigned int) (-k));
                err = (err >> (-k)) + 1 + dropped;
            }
            /* exp(x) < 1: a Y at or above B^nk (by its error, as for a
               tiny x truncated to zero) is clamped to B^nk, which only
               moves it toward the value */
            if (y[nk] != 0)
            {
                _flint_mpn_zero_small(y, nk);
                y[nk] = 1;
            }
            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            *elo = *ehi = err;
            return 0;
        }
    }
    else
    {
        if (!sgnbit)
        {
            /* Y 2^(k - FLINT_BITS nk) - 1 with k >= 1: subtract the bit
               FLINT_BITS nk - k of Y, or a fraction of a unit */
            slong b = FLINT_BITS * nk - k;

            *elo = *ehi = err;
            if (b >= 0)
                mpn_sub_1(y + b / FLINT_BITS, y + b / FLINT_BITS,
                    nk + 1 - b / FLINT_BITS, UWORD(1) << (b % FLINT_BITS));
            else
                *ehi = err + 1;

            return k;
        }
        else
        {
            /* 1 - exp(x), exp(x) < 1/2 (k <= -1): Y 2^k at nk fraction
               limbs (one more unit for the dropped bits) */
            slong s = -k;

            FLINT_ASSERT(k <= -1);

            if (s >= FLINT_BITS * (nk + 1))
            {
                /* exp(x) below a unit: 1 - exp(x) just below 1 */
                _flint_mpn_zero_small(y, nk);
                y[nk] = 1;
                *elo = 1;
                *ehi = 0;
                return 0;
            }
            else
            {
                slong ls = s / FLINT_BITS, i;
                unsigned int bs = s % FLINT_BITS;
                int dropped = 0;

                for (i = 0; i < ls; i++)
                    dropped |= (y[i] != 0);
                if (bs != 0)
                {
                    dropped |= ((y[ls] << (FLINT_BITS - bs)) != 0);
                    mpn_rshift(y, y + ls, nk + 1 - ls, bs);
                }
                else
                    flint_mpn_copyi(y, y + ls, nk + 1 - ls);
                _flint_mpn_zero_small(y + nk + 1 - ls, ls);
                err = (s >= FLINT_BITS) ? 1 : (err >> s) + 1;
                err += dropped;
            }

            mpn_neg(y, y, nk + 1);
            y[nk] += 1;
            *elo = *ehi = err;
            return 0;
        }
    }
}

slong
_mp_real_exp2_fp(nn_ptr y, ulong * err, nn_srcptr d, slong n, slong e, int sgnbit, slong pb)
{
    slong nk = _mp_real_elem_limbs(pb);
    nn_ptr X, T;
    nn_srcptr L;
    slong N, k;
    int trunc;
    TMP_INIT;

    FLINT_ASSERT(e <= FLINT_BITS - 3);

    N = nk + 1;
    TMP_START;
    X = TMP_ALLOC((2 * N + 3) * sizeof(ulong));
    T = X + N + 2;

    /* X = x with one integral limb in two's complement and N fraction
       limbs, truncated towards -inf */
    X[0] = X[N] = 0;    /* (written first only to keep gcc from warning) */
    trunc = _mp_real_fp_get_fixed(X, N + 1, N, d, n, e);

    if (sgnbit)
    {
        mpn_neg(X, X, N + 1);
        if (trunc)
            mpn_sub_1(X, X, N + 1, 1);
    }

    k = (slong) X[N];

    /* an integer: 2^k exactly */
    if (!trunc && flint_mpn_zero_p(X, N))
    {
        _flint_mpn_zero_small(y, nk);
        y[nk] = 1;
        *err = 0;
        TMP_END;
        return k;
    }

    /* t = f L in [0, log 2): the product truncated (2 ulps of B^-N),
       f's truncation (1 ulp) times log 2 */
    L = _mp_real_const_ptr(MP_REAL_CONST_ID_LOG2, N);
    flint_mpn_mulhigh_n(T, X, L, N);

    /* t within two ulps of B^-nk, times exp' < 2 */
    _mp_real_exp_kernel(y, err, T + 1, nk);
    *err += 4;
    TMP_END;
    return k;
}
