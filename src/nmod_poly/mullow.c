/*
    Copyright (C) 2008-2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"

#ifdef FLINT_HAVE_FFT_SMALL

#include "fft_small.h"

/* todo: separate squaring table */
/* todo: check unbalanced cutoffs */

static const short fft_mullow_tab[] = {1115, 1115, 597, 569, 407, 321, 306, 279, 191,
182, 166, 159, 152, 145, 139, 89, 85, 78, 75, 75, 69, 174, 174, 166, 159,
152, 152, 152, 97, 101, 106, 111, 101, 101, 101, 139, 145, 145, 139, 145,
145, 139, 145, 145, 145, 182, 182, 182, 182, 182, 182, 191, 200, 220, 210,
200, 210, 210, 210, 210, 191, 182, 182, 174, };

#endif


void _nmod_poly_mullow(mp_ptr res, mp_srcptr poly1, slong len1,
                             mp_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong bits;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (len2 <= 5)
    {
        _nmod_poly_mullow_classical(res, poly1, len1, poly2, len2, n, mod);
        return;
    }

    if (n == len1 + len2 - 1)
    {
        _nmod_poly_mul(res, poly1, len1, poly2, len2, mod);
        return;
    }

    bits = NMOD_BITS(mod);

#ifdef FLINT_HAVE_FFT_SMALL

    if (len2 >= fft_mullow_tab[bits - 1])
    {
        _nmod_poly_mul_mid_default_mpn_ctx(res, 0, n, poly1, len1, poly2, len2, mod);
        return;
    }

#endif

    if (n < 10 + bits * bits / 10)
        _nmod_poly_mullow_classical(res, poly1, len1, poly2, len2, n, mod);
    else
        _nmod_poly_mullow_KS(res, poly1, len1, poly2, len2, 0, n, mod);
}

void nmod_poly_mullow(nmod_poly_t res,
              const nmod_poly_t poly1, const nmod_poly_t poly2, slong trunc)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;

    if (len1 == 0 || len2 == 0 || trunc == 0)
    {
        nmod_poly_zero(res);

        return;
    }

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, poly1->mod);
        else
            _nmod_poly_mullow(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, poly1->mod);
        else
            _nmod_poly_mullow(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, poly1->mod);
    }

    res->length = trunc;
    _nmod_poly_normalise(res);
}

/* Assumes poly1 and poly2 are not length 0 and 0 < n <= len1 + len2 - 1 */
void
_nmod_poly_mullow_classical(mp_ptr res, mp_srcptr poly1, slong len1,
                            mp_srcptr poly2, slong len2, slong n, nmod_t mod)
{
    slong i, j, bits, log_len, nlimbs, n1, n2;
    int squaring;
    mp_limb_t c;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    if (n == 1)
    {
        res[0] = nmod_mul(poly1[0], poly2[0], mod);
        return;
    }

    if (len2 == 1)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        return;
    }

    squaring = (poly1 == poly2 && len1 == len2);

    log_len = FLINT_BIT_COUNT(len2);
    bits = FLINT_BITS - (slong) mod.norm;
    bits = 2 * bits + log_len;

    if (bits <= FLINT_BITS)
    {
        flint_mpn_zero(res, n);

        if (squaring)
        {
            for (i = 0; i < len1; i++)
            {
                c = poly1[i];

                if (2 * i < n)
                    res[2 * i] += c * c;

                c *= 2;

                for (j = i + 1; j < FLINT_MIN(len1, n - i); j++)
                    res[i + j] += poly1[j] * c;
            }
        }
        else
        {
            for (i = 0; i < len1; i++)
            {
                mp_limb_t c = poly1[i];

                for (j = 0; j < FLINT_MIN(len2, n - i); j++)
                    res[i + j] += c * poly2[j];
            }
        }

        _nmod_vec_reduce(res, res, n, mod);
        return;
    }

    if (len2 == 2)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        _nmod_vec_scalar_addmul_nmod(res + 1, poly1, len1 - 1, poly2[1], mod);
        if (n == len1 + len2 - 1)
            res[len1 + len2 - 2] = nmod_mul(poly1[len1 - 1], poly2[len2 - 1], mod);
        return;
    }

    if (bits <= 2 * FLINT_BITS)
        nlimbs = 2;
    else
        nlimbs = 3;

    if (squaring)
    {
        for (i = 0; i < n; i++)
        {
            n1 = FLINT_MAX(0, i - len1 + 1);
            n2 = FLINT_MIN(len1 - 1, (i + 1) / 2 - 1);

            c = _nmod_vec_dot_rev(poly1 + n1, poly1 + i - n2, n2 - n1 + 1, mod, nlimbs);
            c = nmod_add(c, c, mod);

            if (i % 2 == 0 && i / 2 < len1)
                NMOD_ADDMUL(c, poly1[i / 2], poly1[i / 2], mod);

            res[i] = c;
        }
    }
    else
    {
        for (i = 0; i < n; i++)
        {
            n1 = FLINT_MIN(len1 - 1, i);
            n2 = FLINT_MIN(len2 - 1, i);

            res[i] = _nmod_vec_dot_rev(poly1 + i - n2,
                                       poly2 + i - n1,
                                       n1 + n2 - i + 1, mod, nlimbs);
        }
    }
}

void
nmod_poly_mullow_classical(nmod_poly_t res,
                           const nmod_poly_t poly1, const nmod_poly_t poly2,
                           slong trunc)
{
    slong len_out;

    if (poly1->length == 0 || poly2->length == 0 || trunc == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, trunc);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_classical(temp->coeffs, poly1->coeffs,
                                        poly1->length, poly2->coeffs,
                                        poly2->length, trunc, poly1->mod);
        else
            _nmod_poly_mullow_classical(temp->coeffs, poly2->coeffs,
                                        poly2->length, poly1->coeffs,
                                        poly1->length, trunc, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, trunc);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_classical(res->coeffs, poly1->coeffs,
                                        poly1->length, poly2->coeffs,
                                        poly2->length, trunc, poly1->mod);
        else
            _nmod_poly_mullow_classical(res->coeffs, poly2->coeffs,
                                        poly2->length, poly1->coeffs,
                                        poly1->length, trunc, poly1->mod);
    }

    res->length = trunc;
    _nmod_poly_normalise(res);
}

void
_nmod_poly_mullow_KS(mp_ptr out, mp_srcptr in1, slong len1,
            mp_srcptr in2, slong len2, flint_bitcnt_t bits, slong n, nmod_t mod)
{
    slong limbs1, limbs2;
    mp_ptr tmp, mpn1, mpn2, res;
    int squaring;
    TMP_INIT;

    len1 = FLINT_MIN(len1, n);
    len2 = FLINT_MIN(len2, n);

    squaring = (in1 == in2 && len1 == len2);

    if (bits == 0)
    {
        flint_bitcnt_t bits1, bits2, loglen;

        /* Look at the actual bits of the input? This slows down the generic
        case. Are there situations where we care enough about special input? */
#if 0
        bits1  = _nmod_vec_max_bits2(in1, len1);
        bits2  = squaring ? bits1 : _nmod_vec_max_bits2(in2, len2);
#else
        bits1 = FLINT_BITS - (slong) mod.norm;
        bits2 = bits1;
#endif
        loglen = FLINT_BIT_COUNT(len2);
        bits = bits1 + bits2 + loglen;
    }

    limbs1 = (len1 * bits - 1) / FLINT_BITS + 1;
    limbs2 = (len2 * bits - 1) / FLINT_BITS + 1;

    TMP_START;
    tmp = TMP_ALLOC(sizeof(mp_limb_t) * (limbs1 + limbs2 + limbs1 + (squaring ? 0 : limbs2)));
    res = tmp;
    mpn1 = tmp + limbs1 + limbs2;
    mpn2 = squaring ? mpn1 : (mpn1 + limbs1);

    _nmod_poly_bit_pack(mpn1, in1, len1, bits);
    if (!squaring)
        _nmod_poly_bit_pack(mpn2, in2, len2, bits);

    if (squaring)
        flint_mpn_sqr(res, mpn1, limbs1);
    else
        flint_mpn_mul(res, mpn1, limbs1, mpn2, limbs2);

    _nmod_poly_bit_unpack(out, n, res, bits, mod);

    TMP_END;
}

void
nmod_poly_mullow_KS(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2,
                 flint_bitcnt_t bits, slong n)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0) || n == 0)
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;
    if (n > len_out)
        n = len_out;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_KS(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              n, poly1->mod);
        else
            _nmod_poly_mullow_KS(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              n, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mullow_KS(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              n, poly1->mod);
        else
            _nmod_poly_mullow_KS(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              n, poly1->mod);
    }

    res->length = n;
    _nmod_poly_normalise(res);
}
