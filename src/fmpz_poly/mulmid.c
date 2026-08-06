/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#if FLINT_HAVE_FFT_SMALL
#include "fft_small.h"
#endif

static void
_fmpz_poly_mulmid_tiny1(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong i, top1, top2, start, stop;

    if (poly1 == poly2 && len1 == len2)
    {
        for (i = nlo; i < nhi; i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            slong j, s = 0;

            for (j = 0; j < stop - start + 1; j++)
            {
                s += poly1[start + j] * poly1[i - stop + (stop - start + 1) - 1 - j];
            }

            s *= 2;
            if (i % 2 == 0)
                s += poly1[i / 2] * poly1[i / 2];

            fmpz_clear(res + i - nlo);
            res[i - nlo] = s;
        }
    }
    else
    {
        for (i = nlo; i < nhi; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            slong n = top1 + top2 - i + 1;

            slong j, s;

            s = 0;
            for (j = 0; j < n; j++)
                s += poly1[i - top2 + j] * poly2[i - top1 + n - 1 - j];

            fmpz_clear(res + i - nlo);
            res[i - nlo] = s;
        }
    }
}

static void
_fmpz_poly_mulmid_tiny2(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong i, top1, top2, start, stop;

    if (poly1 == poly2 && len1 == len2)
    {
        for (i = nlo; i < nhi; i++)
        {
            start = FLINT_MAX(0, i - len1 + 1);
            stop = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            slong j;
            ulong s0, s1, t0, t1;

            s0 = s1 = 0;
            for (j = 0; j < stop - start + 1; j++)
            {
                smul_ppmm(t1, t0, poly1[start + j], poly1[i - stop + (stop - start + 1) - 1 - j]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }

            s1 = (s1 << 1) | (s0 >> (FLINT_BITS - 1));
            s0 <<= 1;

            if (i % 2 == 0)
            {
                smul_ppmm(t1, t0, poly1[i / 2], poly1[i / 2]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }

            if (((slong) s1) >= 0)
            {
                fmpz_set_uiui(res + i - nlo, s1, s0);
            }
            else
            {
                sub_ddmmss(s1, s0, 0, 0, s1, s0);
                fmpz_neg_uiui(res + i - nlo, s1, s0);
            }
        }
    }
    else
    {
        for (i = nlo; i < nhi; i++)
        {
            top1 = FLINT_MIN(len1 - 1, i);
            top2 = FLINT_MIN(len2 - 1, i);
            slong n = top1 + top2 - i + 1;

            slong j;
            ulong s0, s1, t0, t1;

            s0 = s1 = 0;
            for (j = 0; j < n; j++)
            {
                smul_ppmm(t1, t0, poly1[i - top2 + j], poly2[i - top1 + n - 1 - j]);
                add_ssaaaa(s1, s0, s1, s0, t1, t0);
            }

            if (((slong) s1) >= 0)
            {
                fmpz_set_uiui(res + i - nlo, s1, s0);
            }
            else
            {
                sub_ddmmss(s1, s0, 0, 0, s1, s0);
                fmpz_neg_uiui(res + i - nlo, s1, s0);
            }
        }
    }
}

void
_fmpz_poly_mulmid(fmpz * res, const fmpz * poly1, slong len1,
                                const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong len, bits1, bits2, rbits;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nhi != 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    /* Low truncation of inputs */
    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    /* High truncation of inputs */
    slong nlo2 = (len1 + len2 - 1) - nlo;

    if (len1 > nlo2)
    {
        slong trunc = len1 - nlo2;
        poly1 += trunc;
        len1 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    if (len2 > nlo2)
    {
        slong trunc = len2 - nlo2;
        poly2 += trunc;
        len2 -= trunc;
        nlo -= trunc;
        nhi -= trunc;
    }

    len = nhi - nlo;

    if (len1 < len2)
    {
        FLINT_SWAP(const fmpz *, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    if (len2 == 1)
    {
        _fmpz_vec_scalar_mul_fmpz(res, poly1 + nlo, len, poly2);
        return;
    }

    bits1 = _fmpz_vec_max_bits(poly1, len1);
    bits2 = (poly1 == poly2 && len1 == len2) ? bits1 : _fmpz_vec_max_bits(poly2, len2);
    bits1 = FLINT_ABS(bits1);
    bits2 = FLINT_ABS(bits2);

#if FLINT_HAVE_FFT_SMALL
    if (len2 >= 100 && (bits1 + bits2 <= 40 || bits1 + bits2 >= 128 || len2 >= 200))
        if (_fmpz_poly_mul_mid_default_mpn_ctx(res, nlo, nhi, poly1, len1, poly2, len2))
            return;
#endif

    if (bits1 <= SMALL_FMPZ_BITCOUNT_MAX && bits2 <= SMALL_FMPZ_BITCOUNT_MAX &&
        (len2 < 50 || (4 * len2 >= 3 * len && len < 150 + bits1 + bits2)))
    {
        rbits = bits1 + bits2 + FLINT_BIT_COUNT(len2);

        if (rbits <= SMALL_FMPZ_BITCOUNT_MAX)
        {
            _fmpz_poly_mulmid_tiny1(res, poly1, len1, poly2, len2, nlo, nhi);
            return;
        }
        else if (rbits <= 2 * FLINT_BITS - 1)
        {
            _fmpz_poly_mulmid_tiny2(res, poly1, len1, poly2, len2, nlo, nhi);
            return;
        }
    }

#if FLINT_HAVE_FFT_SMALL

    /* same as in mul.c */
    if (nhi <= 8 || (len2 <= 6 && FLINT_MIN(bits1, bits2) <= 5000) || len <= 3)
        _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
    else if (nlo == 0 && (len2 <= 4 || (len2 <= 8 && bits1 + bits2 >= 1500 && bits1 + bits2 <= 10000)))
        _fmpz_poly_mullow_karatsuba(res, poly1, len1, poly2, len2, nhi);
    else if
        ((len2 >= 8 && len2 <= 75 && bits1 + bits2 >= 800 && bits1 + bits2 <= 4000) ||
            (len1 + len2 >= 5000 && bits1 + bits2 >= 5000 + (len1 + len2) / 10 && flint_get_num_threads() >= 4))
        _fmpz_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
    else
        _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);

#else

    if (len <= 8 || (len2 < 7))
    {
        _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
    }
    else
    {
        slong limbs1, limbs2;

        limbs1 = (bits1 + FLINT_BITS - 1) / FLINT_BITS;
        limbs2 = (bits2 + FLINT_BITS - 1) / FLINT_BITS;

        /* todo: mulmid_karatsuba */
        if (nlo == 0 && nhi < 16 && (limbs1 > 12 || limbs2 > 12))
            _fmpz_poly_mullow_karatsuba(res, poly1, len1, poly2, len2, nhi);
        else if (limbs1 + limbs2 <= 8)
            _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
        else if ((limbs1+limbs2)/2048 > len1 + len2)
            _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
        else if ((limbs1 + limbs2)*FLINT_BITS*4 < len1 + len2)
            _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
        else
            _fmpz_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
    }

#endif
}

void
fmpz_poly_mulmid(fmpz_poly_t res, const fmpz_poly_t poly1,
                                            const fmpz_poly_t poly2, slong nlo, slong nhi)
{
    slong len1 = fmpz_poly_length(poly1);
    slong len2 = fmpz_poly_length(poly2);
    slong len;

    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nhi >= 0);

    if (len1 == 0 || len2 == 0 || nlo >= FLINT_MIN(nhi, len1 + len2 - 1))
    {
        fmpz_poly_zero(res);
        return;
    }

    nhi = FLINT_MIN(nhi, len1 + len2 - 1);
    len = nhi - nlo;

    if (res == poly1 || res == poly2)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, len);
        _fmpz_poly_mulmid(t->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, nlo, nhi);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(res, len);
        _fmpz_poly_mulmid(res->coeffs, poly1->coeffs, poly1->length,
                                    poly2->coeffs, poly2->length, nlo, nhi);
    }

    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}

