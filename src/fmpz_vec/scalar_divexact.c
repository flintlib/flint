/*
    Copyright 2000-2003, 2005, 2013 Free Software Foundation, Inc.
    Copyright 1991-2018, 2021, 2022 Free Software Foundation, Inc.
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* This file includes code adapted from GMP (gmp-impl.h, mpn/generic/dive_1.c). */

#include <gmp.h>
#include "longlong.h"
#include "fmpz.h"
#include "fmpz_vec.h"

/* Helper functions that could go in longlong.h, ulong_extras.h
   and mpn_extras.h as soon as they find a use elsewhere. */

/* Tuned for Zen 3; with assembly code we would not need GMP. */
#define DIVEXACT_1_ODD_GMP_CUTOFF 35
#define DIVEXACT_1_EVEN_GMP_CUTOFF 25

/* The table is slightly faster but requires GMP internals. */
#define USE_BINVERT_TABLE 0

#if USE_BINVERT_TABLE
extern const unsigned char  __gmp_binvert_limb_table[128];

/* Invert n mod 2^FLINT_BITS */
static ulong n_binvert(ulong n)
{
    ulong r;

    FLINT_ASSERT(n % 2);

    r = __gmp_binvert_limb_table[(n / 2) & 0x7F];
    r = 2 * r - r * r * n;          /* 8 bits */
    r = 2 * r - r * r * n;          /* 16 bits */
    r = 2 * r - r * r * n;          /* 32 bits */
#if FLINT_BITS == 64
    r = 2 * r - r * r * n;          /* 64 bits */
#endif

    return r;
}
#else

/* Use Hurchalla's algorithm https://arxiv.org/abs/2204.04342 */
static ulong n_binvert(ulong a)
{
    ulong r, y;

    r = (3 * a) ^ 2;    /* 5 bits */
    y = 1 - a * r;
    r = r * (1 + y);    /* 10 bits */
    y *= y;
    r = r * (1 + y);    /* 20 bits */
    y *= y;
    r = r * (1 + y);    /* 40 bits */
#if FLINT_BITS == 64
    y *= y;
    r = r * (1 + y);    /* 80 bits */
#endif
    return r;
}
#endif

FLINT_FORCE_INLINE
ulong n_mulhi(ulong a, ulong b)
{
    ulong h, l;
    umul_ppmm(h, l, a, b);
    return h;
}

FLINT_FORCE_INLINE
ulong n_subc(ulong * c, ulong a, ulong b)
{
    ulong d = a - b;
    *c = d > a;     /* the compiler should recognize this idiom */
    return d;
}

static
void nn_divexact_2_1_odd(nn_ptr res, nn_srcptr x, ulong d, ulong dinv)
{
    ulong s, l, c;

    s = x[0];
    l = s * dinv;
    res[0] = l;
    c = n_mulhi(l, d);
    s = x[1];
    l = s - c;
    l = l * dinv;
    res[1] = l;
}

static
void nn_divexact_2_1_even(nn_ptr res, nn_srcptr x, ulong d, ulong dinv, unsigned int norm)
{
    ulong s, l, c;
    ulong ls, s_next;
    ulong h;

    d >>= norm;
    c = 0;
    s = x[0];
    s_next = x[1];
    ls = (s >> norm) | (s_next << (FLINT_BITS - norm));
    s = s_next;
    l = n_subc(&c, ls, c);
    l = l * dinv;
	res[0] = l;
    h = n_mulhi(l, d);
    c += h;
    ls = s >> norm;
    l = ls - c;
    l = l * dinv;
    res[1] = l;
}

static
void nn_divexact_1_odd(nn_ptr res, nn_srcptr x, slong xn, ulong d, ulong dinv)
{
    ulong s, l, c;
    slong i;

    if (xn >= DIVEXACT_1_ODD_GMP_CUTOFF)
    {
        mpn_divexact_1(res, x, xn, d);
        return;
    }

    s = x[0];
    l = s * dinv;
    res[0] = l;
    c = 0;

    for (i = 1; i < xn; i++)
    {
        c += n_mulhi(l, d);
        s = x[i];
        l = n_subc(&c, s, c);
        l = l * dinv;
        res[i] = l;
    }
}

static
void nn_divexact_1_even(nn_ptr res, nn_srcptr x, slong xn, ulong d, ulong dinv, unsigned int norm)
{
    ulong s, l, c, h;
    slong i;
    ulong ls, s_next;

    if (xn >= DIVEXACT_1_EVEN_GMP_CUTOFF)
    {
        mpn_divexact_1(res, x, xn, d);
        return;
    }

    d >>= norm;
    c = 0;
    s = x[0];

    for (i = 1; i < xn; i++)
    {
        s_next = x[i];
        ls = (s >> norm) | (s_next << (FLINT_BITS - norm));
        s = s_next;
        l = n_subc(&c, ls, c);
        l = l * dinv;
        res[i - 1] = l;
        h = n_mulhi(l, d);
        c += h;
    }

    ls = s >> norm;
    l = ls - c;
    l = l * dinv;
    res[xn - 1] = l;
}

#define FLINT_MPZ_GET_MPN(zd, zn, zsgnbit, z) \
    do { \
        (zsgnbit) = (z)->_mp_size < 0; \
        (zn) = FLINT_ABS((z)->_mp_size); \
        (zd) = (z)->_mp_d; \
    } \
    while (0)

/* Disable for a small speedup, which however risks creating corrupted fmpzs
   if the user somehow did not provide exact quotients. */
#define SAFE_VERSION 1

static void
_fmpz_vec_divexact_ui(fmpz * res, const fmpz * x, slong len, ulong c, int negc)
{
    ulong cinv;
    slong i;
    unsigned int norm;
    mpz_ptr z;
    slong xn;
    int negative;
    nn_srcptr xd;
    ulong cinv_signed;

    FLINT_ASSERT(c != 0);

    if (c == 1)
    {
        if (negc)
            _fmpz_vec_neg(res, x, len);
        else
            _fmpz_vec_set(res, x, len);
        return;
    }

    if (c & 1)
    {
        cinv = n_binvert(c);
        cinv_signed = negc ? -cinv : cinv;

        for (i = 0; i < len; i++)
        {
            slong v = x[i];

            if (!COEFF_IS_MPZ(v))
            {
#if SAFE_VERSION
                fmpz_set_si(res + i, v * cinv_signed);
#else
                _fmpz_demote(res + i);
                res[i] = v * cinv_signed;
#endif
            }
            else
            {
                z = COEFF_TO_PTR(v);
                FLINT_MPZ_GET_MPN(xd, xn, negative, z);
                negative ^= negc;

                if (xn == 1)
                {
                    /* cannot overflow if c >= 2 */
                    fmpz_set_si(res + i, negative ? -xd[0] * cinv : xd[0] * cinv);
                }
                else if (xn == 2)
                {
                    ulong t[2];
                    nn_divexact_2_1_odd(t, xd, c, cinv);

                    if (negative)
                        fmpz_neg_uiui(res + i, t[1], t[0]);
                    else
                        fmpz_set_uiui(res + i, t[1], t[0]);
                }
                else
                {
                    mpz_ptr rp = _fmpz_promote(res + i);
                    nn_ptr zd = FLINT_MPZ_REALLOC(rp, xn);

                    nn_divexact_1_odd(zd, xd, xn, c, cinv);

                    xn -= (zd[xn - 1] == 0);
                    FLINT_ASSERT(zd[xn - 1] != 0);
                    rp->_mp_size = negative ? -xn : xn;
#if SAFE_VERSION
                    _fmpz_demote_val(res + i);
#endif
                }
            }
        }
    }
    else
    {
        norm = flint_ctz(c);

        /* Special case for powers of 2 */
        if ((c & (c - 1)) == 0)
        {
            if (negc)
            {
                for (i = 0; i < len; i++)
                {
                    slong v = x[i];

                    if (!COEFF_IS_MPZ(v))
                    {
                        _fmpz_demote(res + i);
                        res[i] = (-v) >> norm;
                    }
                    else
                    {
                        fmpz_tdiv_q_2exp(res + i, x + i, norm);
                        fmpz_neg(res + i, res + i);
                    }
                }
            }
            else
            {
                for (i = 0; i < len; i++)
                {
                    slong v = x[i];

                    if (!COEFF_IS_MPZ(v))
                    {
                        _fmpz_demote(res + i);
                        res[i] = v >> norm;
                    }
                    else
                        fmpz_tdiv_q_2exp(res + i, x + i, norm);
                }
            }

            return;
        }

        cinv = n_binvert(c >> norm);
        cinv_signed = negc ? -cinv : cinv;

        for (i = 0; i < len; i++)
        {
            slong v = x[i];

            if (!COEFF_IS_MPZ(v))
            {
#if SAFE_VERSION
                fmpz_set_si(res + i, ((slong) (v * cinv_signed)) >> norm);
#else
                _fmpz_demote(res + i);
                res[i] = ((slong) (v * cinv_signed)) >> norm;
#endif
            }
            else
            {
                z = COEFF_TO_PTR(v);
                FLINT_MPZ_GET_MPN(xd, xn, negative, z);
                negative ^= negc;

                if (xn == 1)
                {
                    if (negative)
                        fmpz_neg_ui(res + i, (xd[0] * cinv) >> norm);
                    else
                        fmpz_set_ui(res + i, (xd[0] * cinv) >> norm);
                }
                else if (xn == 2)
                {
                    ulong t[2];
                    nn_divexact_2_1_even(t, xd, c, cinv, norm);

                    if (negative)
                        fmpz_neg_uiui(res + i, t[1], t[0]);
                    else
                        fmpz_set_uiui(res + i, t[1], t[0]);
                }
                else
                {
                    mpz_ptr rp = _fmpz_promote(res + i);
                    nn_ptr zd = FLINT_MPZ_REALLOC(rp, xn);

                    nn_divexact_1_even(zd, xd, xn, c, cinv, norm);
                    xn -= (zd[xn - 1] == 0);
                    FLINT_ASSERT(zd[xn - 1] != 0);
                    rp->_mp_size = negative ? -xn : xn;
#if SAFE_VERSION
                    _fmpz_demote_val(res + i);
#endif
                }
            }
        }
    }
}

void
_fmpz_vec_scalar_divexact_si(fmpz * vec1, const fmpz * vec2, slong len2, slong c)
{
    if (len2 == 1)
    {
        fmpz_divexact_si(vec1, vec2, c);
        return;
    }

    if (c > 0)
        _fmpz_vec_divexact_ui(vec1, vec2, len2, c, 0);
    else
        _fmpz_vec_divexact_ui(vec1, vec2, len2, -(ulong) c, 1);
}

void
_fmpz_vec_scalar_divexact_ui(fmpz * vec1, const fmpz * vec2,
                             slong len2, ulong c)
{
    if (len2 == 1)
    {
        fmpz_divexact_ui(vec1, vec2, c);
        return;
    }

    _fmpz_vec_divexact_ui(vec1, vec2, len2, c, 0);
}

void
_fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, const fmpz * vec2,
                               slong len2, const fmpz_t x)
{
    fmpz c = *x;

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 1)
            _fmpz_vec_set(vec1, vec2, len2);
        else if (c == -1)
            _fmpz_vec_neg(vec1, vec2, len2);
        else
            _fmpz_vec_scalar_divexact_si(vec1, vec2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
            fmpz_divexact(vec1 + i, vec2 + i, x);
    }
}

