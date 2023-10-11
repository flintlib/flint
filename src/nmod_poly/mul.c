/*
    Copyright (C) 2007, 2008 David Harvey (zn_poly)
    Copyright (C) 2010, 2013 William Hart
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

/* todo: check unbalanced cutoffs */

static const short fft_mul_tab[] = {1326, 1326, 1095, 802, 674, 537, 330, 306, 290,
274, 200, 192, 182, 173, 163, 99, 97, 93, 90, 82, 80, 438, 414, 324, 393,
298, 298, 268, 187, 185, 176, 176, 168, 167, 158, 158, 97, 96, 93, 92, 89,
89, 85, 85, 80, 81, 177, 172, 163, 162, 164, 176, 171, 167, 167, 164, 163,
163, 160, 165, 95, 96, 90, 94, };

static const short fft_sqr_tab[] = {1420, 1420, 1353, 964, 689, 569, 407, 353, 321,
321, 292, 279, 200, 182, 182, 159, 159, 152, 145, 139, 723, 626, 626, 569,
597, 448, 542, 292, 292, 200, 191, 191, 182, 182, 166, 166, 166, 159, 159,
159, 152, 152, 145, 145, 93, 200, 191, 182, 182, 182, 182, 191, 191, 191,
182, 182, 174, 182, 182, 182, 152, 152, 152, 145, };

#endif

void _nmod_poly_mul(mp_ptr res, mp_srcptr poly1, slong len1,
                             mp_srcptr poly2, slong len2, nmod_t mod)
{
    slong bits, cutoff_len;

    if (len2 <= 5)
    {
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
        return;
    }

    bits = NMOD_BITS(mod);
    cutoff_len = FLINT_MIN(len1, 2 * len2);

#ifdef FLINT_HAVE_FFT_SMALL

    if (poly1 == poly2 && len1 == len2)
    {
        if (cutoff_len >= fft_sqr_tab[bits - 1])
        {
            _nmod_poly_mul_mid_default_mpn_ctx(res, 0, len1 + len2 - 1, poly1, len1, poly2, len2, mod);
            return;
        }
    }
    else
    {
        if (cutoff_len >= fft_mul_tab[bits - 1])
        {
            _nmod_poly_mul_mid_default_mpn_ctx(res, 0, len1 + len2 - 1, poly1, len1, poly2, len2, mod);
            return;
        }
    }

#endif

    if (3 * cutoff_len < 2 * FLINT_MAX(bits, 10))
        _nmod_poly_mul_classical(res, poly1, len1, poly2, len2, mod);
    else if (cutoff_len * bits < 800)
        _nmod_poly_mul_KS(res, poly1, len1, poly2, len2, 0, mod);
    else if (cutoff_len * (bits + 1) * (bits + 1) < 100000)
        _nmod_poly_mul_KS2(res, poly1, len1, poly2, len2, mod);
    else
        _nmod_poly_mul_KS4(res, poly1, len1, poly2, len2, mod);
}

void nmod_poly_mul(nmod_poly_t res, const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        nmod_poly_zero(res);

        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;

        nmod_poly_init2(temp, poly1->mod.n, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);

        nmod_poly_swap(temp, res);
        nmod_poly_clear(temp);
    } else
    {
        nmod_poly_fit_length(res, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, poly1->mod);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

/* Assumes poly1 and poly2 are not length 0. */
void
_nmod_poly_mul_classical(mp_ptr res, mp_srcptr poly1,
                         slong len1, mp_srcptr poly2, slong len2, nmod_t mod)
{
    slong i, j, bits, log_len, nlimbs, n1, n2;
    int squaring;
    mp_limb_t c;

    if (len1 == 1)
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
        flint_mpn_zero(res, len1 + len2 - 1);

        if (squaring)
        {
            for (i = 0; i < len1; i++)
            {
                c = poly1[i];

                res[2 * i] += c * c;
                c *= 2;

                for (j = i + 1; j < len1; j++)
                    res[i + j] += poly1[j] * c;
            }
        }
        else
        {
            for (i = 0; i < len1; i++)
            {
                mp_limb_t c = poly1[i];

                for (j = 0; j < len2; j++)
                    res[i + j] += c * poly2[j];
            }
        }

        _nmod_vec_reduce(res, res, len1 + len2 - 1, mod);
        return;
    }

    if (len2 == 2)
    {
        _nmod_vec_scalar_mul_nmod(res, poly1, len1, poly2[0], mod);
        _nmod_vec_scalar_addmul_nmod(res + 1, poly1, len1 - 1, poly2[1], mod);
        res[len1 + len2 - 2] = nmod_mul(poly1[len1 - 1], poly2[len2 - 1], mod);
        return;
    }

    if (bits <= 2 * FLINT_BITS)
        nlimbs = 2;
    else
        nlimbs = 3;

    if (squaring)
    {
        for (i = 0; i < 2 * len1 - 1; i++)
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
        for (i = 0; i < len1 + len2 - 1; i++)
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
nmod_poly_mul_classical(nmod_poly_t res,
                        const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_classical(temp->coeffs, poly1->coeffs,
                                     poly1->length, poly2->coeffs,
                                     poly2->length, poly1->mod);
        else
            _nmod_poly_mul_classical(temp->coeffs, poly2->coeffs,
                                     poly2->length, poly1->coeffs,
                                     poly1->length, poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_classical(res->coeffs, poly1->coeffs, poly1->length,
                                     poly2->coeffs, poly2->length, poly1->mod);
        else
            _nmod_poly_mul_classical(res->coeffs, poly2->coeffs, poly2->length,
                                     poly1->coeffs, poly1->length, poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

void
_nmod_poly_mul_KS(mp_ptr out, mp_srcptr in1, slong len1,
                  mp_srcptr in2, slong len2, flint_bitcnt_t bits, nmod_t mod)
{
    slong len_out = len1 + len2 - 1, limbs1, limbs2;
    mp_ptr tmp, mpn1, mpn2, res;
    int squaring;
    TMP_INIT;

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

    _nmod_poly_bit_unpack(out, len_out, res, bits, mod);

    TMP_END;
}

void
nmod_poly_mul_KS(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2,
                 flint_bitcnt_t bits)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_KS(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length, bits,
                              poly1->mod);
        else
            _nmod_poly_mul_KS(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length, bits,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

/*
   Multiplication/squaring using Kronecker substitution at 2^b and -2^b.
*/
void
_nmod_poly_mul_KS2(mp_ptr res, mp_srcptr op1, slong n1,
                  mp_srcptr op2, slong n2, nmod_t mod)
{
   int sqr, v3m_neg;
   ulong bits, b, w;
   slong n1o, n1e, n2o, n2e, n3o, n3e, n3, k1, k2, k3;
   mp_ptr v1_buf0, v2_buf0, v1_buf1, v2_buf1, v1_buf2, v2_buf2;
   mp_ptr v1o, v1e, v1p, v1m, v2o, v2e, v2p, v2m, v3o, v3e, v3p, v3m;
   mp_ptr z;
   TMP_INIT;

   if (n2 == 1)
   {
      /* code below needs n2 > 1, so fall back on scalar multiplication */
      _nmod_vec_scalar_mul_nmod(res, op1, n1, op2[0], mod);
      return;
   }

   TMP_START;

   sqr = (op1 == op2 && n1 == n2);

   /* bits in each output coefficient */
   bits = 2 * (FLINT_BITS - mod.norm) + FLINT_CLOG2(n2);

   /* we're evaluating at x = B and -B, where B = 2^b, and b = ceil(bits / 2) */
   b = (bits + 1) / 2;

   /* number of ulongs required to store each output coefficient */
   w = (2*b - 1)/FLINT_BITS + 1;

   /*
      Write f1(x) = f1e(x^2) + x * f1o(x^2)
            f2(x) = f2e(x^2) + x * f2o(x^2)
             h(x) =  he(x^2) + x *  ho(x^2)
      "e" = even, "o" = odd
   */

   n1o = n1 / 2;
   n1e = n1 - n1o;

   n2o = n2 / 2;
   n2e = n2 - n2o;

   n3 = n1 + n2 - 1;    /* length of h */
   n3o = n3 / 2;
   n3e = n3 - n3o;

   /*
      f1(B) and |f1(-B)| are at most ((n1 - 1) * b + mod->bits) bits long.
      However, when evaluating f1e(B^2) and B * f1o(B^2) the bitpacking
      routine needs room for the last chunk of 2b bits. Therefore we need to
      allow room for (n1 + 1) * b bits. Ditto for f2.
   */
   k1 = ((n1 + 1)*b - 1)/FLINT_BITS + 1;
   k2 = ((n2 + 1)*b - 1)/FLINT_BITS + 1;
   k3 = k1 + k2;

   /* allocate space */
   v1_buf0 = TMP_ALLOC(sizeof(mp_limb_t) * 3 * k3); /* k1 limbs */
   v2_buf0 = v1_buf0 + k1;         /* k2 limbs */
   v1_buf1 = v2_buf0 + k2;         /* k1 limbs */
   v2_buf1 = v1_buf1 + k1;         /* k2 limbs */
   v1_buf2 = v2_buf1 + k2;         /* k1 limbs */
   v2_buf2 = v1_buf2 + k1;         /* k2 limbs */

   /*
      arrange overlapping buffers to minimise memory use
      "p" = plus, "m" = minus
   */
   v1e = v1_buf0;
   v2e = v2_buf0;
   v1o = v1_buf1;
   v2o = v2_buf1;
   v1p = v1_buf2;
   v2p = v2_buf2;
   v1m = v1_buf0;
   v2m = v2_buf0;
   v3m = v1_buf1;
   v3p = v1_buf0;
   v3e = v1_buf2;
   v3o = v1_buf0;

   z = TMP_ALLOC(sizeof(mp_limb_t) * w * n3e);

   if (!sqr)
   {
      /* multiplication version */

      /* evaluate f1e(B^2) and B * f1o(B^2) */
      _nmod_poly_KS2_pack(v1e, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      /* evaluate f2e(B^2) and B * f2o(B^2) */
      _nmod_poly_KS2_pack(v2e, op2, n2e, 2, 2 * b, 0, k2);
      _nmod_poly_KS2_pack(v2o, op2 + 1, n2o, 2, 2 * b, b, k2);

      /*
         compute f1(B) = f1e(B^2) + B * f1o(B^2)
             and f2(B) = f2e(B^2) + B * f2o(B^2)
      */
      mpn_add_n(v1p, v1e, v1o, k1);
      mpn_add_n(v2p, v2e, v2o, k2);

      /*
         compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
             and |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
      */
      v3m_neg  = signed_mpn_sub_n(v1m, v1e, v1o, k1);
      v3m_neg ^= signed_mpn_sub_n(v2m, v2e, v2o, k2);

      /*
         compute  h(B)   =  f1(B)   *  f2(B)
         compute |h(-B)| = |f1(-B)| * |f2(-B)|
         v3m_neg is set if h(-B) is negative
      */
      flint_mpn_mul(v3m, v1m, k1, v2m, k2);
      flint_mpn_mul(v3p, v1p, k1, v2p, k2);
   }
   else
   {
      /* squaring version */

      /* evaluate f1e(B^2) and B * f1o(B^2) */
      _nmod_poly_KS2_pack(v1e, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1o, op1 + 1, n1o, 2, 2 * b, b, k1);

      /* compute f1(B) = f1e(B^2) + B * f1o(B^2) */
      mpn_add_n(v1p, v1e, v1o, k1);

      /* compute |f1(-B)| = |f1e(B^2) - B * f1o(B^2)| */
      signed_mpn_sub_n(v1m, v1e, v1o, k1);

      /*
         compute h(B)  = f1(B)^2
         compute h(-B) = f1(-B)^2
         v3m_neg is cleared (since f1(-B)^2 is never negative)
      */
      flint_mpn_sqr(v3m, v1m, k1);
      flint_mpn_sqr(v3p, v1p, k1);
      v3m_neg = 0;
   }

   /*
      he(B^2) and B * ho(B^2) are both at most b * (n3 + 1) bits long (since
      the coefficients don't overlap). The buffers used below are at least
      b * (n1 + n2 + 2) = b * (n3 + 3) bits long. So we definitely have
      enough room for 2 * he(B^2) and 2 * B * ho(B^2).
   */

   /* compute 2 * he(B^2) = h(B) + h(-B) */
   if (v3m_neg)
      mpn_sub_n(v3e, v3p, v3m, k3);
   else
      mpn_add_n(v3e, v3p, v3m, k3);

   /* unpack coefficients of he, and reduce mod m */
   _nmod_poly_KS2_unpack(z, v3e, n3e, 2 * b, 1);
   _nmod_poly_KS2_reduce(res, 2, z, n3e, w, mod);

   /* compute 2 * b * ho(B^2) = h(B) - h(-B) */
   if (v3m_neg)
      mpn_add_n(v3o, v3p, v3m, k3);
   else
      mpn_sub_n(v3o, v3p, v3m, k3);

   /* unpack coefficients of ho, and reduce mod m */
   _nmod_poly_KS2_unpack(z, v3o, n3o, 2 * b, b + 1);
   _nmod_poly_KS2_reduce(res + 1, 2, z, n3o, w, mod);

   TMP_END;
}

void
nmod_poly_mul_KS2(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS2(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS2(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS2(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS2(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}

/*
   Multiplication/squaring using Kronecker substitution at 2^b, -2^b,
   2^(-b) and -2^(-b).
*/
void
_nmod_poly_mul_KS4(mp_ptr res, mp_srcptr op1, slong n1,
                  mp_srcptr op2, slong n2, nmod_t mod)
{
   int sqr, v3m_neg;
   ulong bits, b, w, a1, a2, a3;
   slong n1o, n1e, n2o, n2e, n3o, n3e, n3, k1, k2, k3;
   mp_ptr v1_buf0, v2_buf0, v1_buf1, v2_buf1, v1_buf2, v2_buf2, v1_buf3, v2_buf3, v1_buf4, v2_buf4;
   mp_ptr v1on, v1en, v1pn, v1mn, v2on, v2en, v2pn, v2mn, v3on, v3en, v3pn, v3mn;
   mp_ptr v1or, v1er, v1pr, v1mr, v2or, v2er, v2pr, v2mr, v3or, v3er, v3pr, v3mr;
   mp_ptr z, zn, zr;
   TMP_INIT;

   if (n2 == 1)
   {
      /* code below needs n2 > 1, so fall back on scalar multiplication */
      _nmod_vec_scalar_mul_nmod(res, op1, n1, op2[0], mod);
      return;
   }

   TMP_START;

   sqr = (op1 == op2 && n1 == n2);

   /* bits in each output coefficient */
   bits = 2 * (FLINT_BITS - mod.norm) + FLINT_CLOG2(n2);

   /*
      we're evaluating at x = B, -B, 1/B, -1/B,
      where B = 2^b, and b = ceil(bits / 4)
   */
   b = (bits + 3) / 4;

   /* number of ulongs required to store each base-B^2 digit */
   w = (2*b - 1)/FLINT_BITS + 1;

   /*
      Write f1(x) = f1e(x^2) + x * f1o(x^2)
            f2(x) = f2e(x^2) + x * f2o(x^2)
             h(x) =  he(x^2) + x *  ho(x^2)
      "e" = even, "o" = odd
   */

   n1o = n1 / 2;
   n1e = n1 - n1o;

   n2o = n2 / 2;
   n2e = n2 - n2o;

   n3 = n1 + n2 - 1;   /* length of h */
   n3o = n3 / 2;
   n3e = n3 - n3o;

   /*
      Put k1 = number of limbs needed to store f1(B) and |f1(-B)|.
      In f1(B), the leading coefficient starts at bit position b * (n1 - 1)
      and has length 2b, and the coefficients overlap so we need an extra bit
      for the carry: this gives (n1 + 1) * b + 1 bits. Ditto for f2.
   */
   k1 = ((n1 + 1) * b)/FLINT_BITS + 1;
   k2 = ((n2 + 1) * b)/FLINT_BITS + 1;
   k3 = k1 + k2;

   /* allocate space */
   v1_buf0 = TMP_ALLOC(sizeof(mp_limb_t) * 5 * k3); /* k1 limbs */
   v2_buf0 = v1_buf0 + k1;         /* k2 limbs */
   v1_buf1 = v2_buf0 + k2;         /* k1 limbs */
   v2_buf1 = v1_buf1 + k1;         /* k2 limbs */
   v1_buf2 = v2_buf1 + k2;         /* k1 limbs */
   v2_buf2 = v1_buf2 + k1;         /* k2 limbs */
   v1_buf3 = v2_buf2 + k2;         /* k1 limbs */
   v2_buf3 = v1_buf3 + k1;         /* k2 limbs */
   v1_buf4 = v2_buf3 + k2;         /* k1 limbs */
   v2_buf4 = v1_buf4 + k1;         /* k2 limbs */

   /*
      arrange overlapping buffers to minimise memory use
      "p" = plus, "m" = minus
      "n" = normal order, "r" = reciprocal order
   */
   v1en = v1_buf0;
   v1on = v1_buf1;
   v1pn = v1_buf2;
   v1mn = v1_buf0;
   v2en = v2_buf0;
   v2on = v2_buf1;
   v2pn = v2_buf2;
   v2mn = v2_buf0;
   v3pn = v1_buf1;
   v3mn = v1_buf2;
   v3en = v1_buf0;
   v3on = v1_buf1;

   v1er = v1_buf2;
   v1or = v1_buf3;
   v1pr = v1_buf4;
   v1mr = v1_buf2;
   v2er = v2_buf2;
   v2or = v2_buf3;
   v2pr = v2_buf4;
   v2mr = v2_buf2;
   v3pr = v1_buf3;
   v3mr = v1_buf4;
   v3er = v1_buf2;
   v3or = v1_buf3;

   z = TMP_ALLOC(sizeof(mp_limb_t) * 2*w*(n3e + 1));
   zn = z;
   zr = z + w*(n3e + 1);

   /* -------------------------------------------------------------------------
          "normal" evaluation points
   */

   if (!sqr)
   {
      /* multiplication version */

      /*
         evaluate f1e(B^2) and B * f1o(B^2)
         We need max(2 * b*n1e, 2 * b*n1o + b) bits for this packing step,
         which is safe since (n1 + 1) * b + 1 >= max(2 * b*n1e, 2 * b*n1o + b).
         Ditto for f2 below.
      */
      _nmod_poly_KS2_pack(v1en, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1on, op1 + 1, n1o, 2, 2 * b, b, k1);

      /*
         compute  f1(B)  =  f1e(B^2) + B * f1o(B^2)
            and |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      */
      mpn_add_n (v1pn, v1en, v1on, k1);
      v3m_neg = signed_mpn_sub_n(v1mn, v1en, v1on, k1);

      /* evaluate f2e(B^2) and B * f2o(B^2) */
      _nmod_poly_KS2_pack(v2en, op2, n2e, 2, 2 * b, 0, k2);
      _nmod_poly_KS2_pack(v2on, op2 + 1, n2o, 2, 2 * b, b, k2);

      /*
         compute  f2(B)  =  f2e(B^2) + B * f2o(B^2)
            and |f2(-B)| = |f2e(B^2) - B * f2o(B^2)|
      */
      mpn_add_n(v2pn, v2en, v2on, k2);
      v3m_neg ^= signed_mpn_sub_n(v2mn, v2en, v2on, k2);

      /*
         compute  h(B)  =  f1(B)   *  f2(B)
            and |h(-B)| = |f1(-B)| * |f2(-B)|
         hn_neg is set if h(-B) is negative
      */
      flint_mpn_mul(v3pn, v1pn, k1, v2pn, k2);
      flint_mpn_mul(v3mn, v1mn, k1, v2mn, k2);
   }
   else
   {
      /* squaring version */

      /* evaluate f1e(B^2) and B * f1o(B^2) */
      _nmod_poly_KS2_pack(v1en, op1, n1e, 2, 2 * b, 0, k1);
      _nmod_poly_KS2_pack(v1on, op1 + 1, n1o, 2, 2 * b, b, k1);

      /*
         compute  f1(B)  =  f1e(B^2) + B * f1o(B^2)
            and |f1(-B)| = |f1e(B^2) - B * f1o(B^2)|
      */
      mpn_add_n (v1pn, v1en, v1on, k1);
      signed_mpn_sub_n(v1mn, v1en, v1on, k1);

      /*
         compute h(B) =  f1(B)^2
            and h(-B) = |f1(-B)|^2
         hn_neg is cleared since h(-B) is never negative
      */
      flint_mpn_sqr(v3pn, v1pn, k1);
      flint_mpn_sqr(v3mn, v1mn, k1);
      v3m_neg = 0;
   }

   /*
      Each coefficient of h(B) is up to 4b bits long, so h(B) needs at most
      ((n1 + n2 + 2) * b + 1) bits. (The extra +1 is to accommodate carries
      generated by overlapping coefficients.)  The buffer has at least
      ((n1 + n2 + 2) * b + 2) bits. Therefore we can safely store 2*h(B) etc.
   */

   /*
      compute     2 * he(B^2) = h(B) + h(-B)
      and     B * 2 * ho(B^2) = h(B) - h(-B)
   */
   if (v3m_neg)
   {
      mpn_sub_n(v3en, v3pn, v3mn, k3);
      mpn_add_n (v3on, v3pn, v3mn, k3);
   }
   else
   {
      mpn_add_n (v3en, v3pn, v3mn, k3);
      mpn_sub_n (v3on, v3pn, v3mn, k3);
   }

   /* -------------------------------------------------------------------------
          "reciprocal" evaluation points
   */

   /*
      correction factors to take into account that if a polynomial has even
      length, its even and odd coefficients are swapped when the polynomial
      is reversed
   */
   a1 = (n1 & 1) ? 0 : b;
   a2 = (n2 & 1) ? 0 : b;
   a3 = (n3 & 1) ? 0 : b;

   if (!sqr)
   {
      /* multiplication version */

      /* evaluate B^(n1-1) * f1e(1/B^2) and B^(n1-2) * f1o(1/B^2) */
      _nmod_poly_KS2_pack(v1er, op1 + 2*(n1e - 1), n1e, -2, 2 * b, a1, k1);
      _nmod_poly_KS2_pack(v1or, op1 + 1 + 2*(n1o - 1), n1o, -2, 2 * b, b - a1, k1);

      /*
         compute  B^(n1-1) * f1(1/B) =
                      B^(n1-1) * f1e(1/B^2) + B^(n1-2) * f1o(1/B^2)
            and  |B^(n1-1) * f1(-1/B)| =
                     |B^(n1-1) * f1e(1/B^2) - B^(n1-2) * f1o(1/B^2)|
      */
      mpn_add_n(v1pr, v1er, v1or, k1);
      v3m_neg = signed_mpn_sub_n(v1mr, v1er, v1or, k1);

      /* evaluate B^(n2-1) * f2e(1/B^2) and B^(n2-2) * f2o(1/B^2) */
      _nmod_poly_KS2_pack(v2er, op2 + 2*(n2e - 1), n2e, -2, 2 * b, a2, k2);
      _nmod_poly_KS2_pack(v2or, op2 + 1 + 2*(n2o - 1), n2o, -2, 2 * b, b - a2, k2);

      /*
         compute  B^(n2-1) * f2(1/B) =
                      B^(n2-1) * f2e(1/B^2) + B^(n2-2) * f2o(1/B^2)
            and  |B^(n1-1) * f2(-1/B)| =
                     |B^(n2-1) * f2e(1/B^2) - B^(n2-2) * f2o(1/B^2)|
      */
      mpn_add_n (v2pr, v2er, v2or, k2);
      v3m_neg ^= signed_mpn_sub_n(v2mr, v2er, v2or, k2);

      /*
         compute B^(n3-1) * h(1/B) =
                         (B^(n1-1) * f1(1/B)) * (B^(n2-1) * f2(1/B))
             and |B^(n3-1) * h(-1/B)| =
                         |B^(n1-1) * f1(-1/B)| * |B^(n2-1) * f2(-1/B)|
         hr_neg is set if h(-1/B) is negative
      */
      flint_mpn_mul(v3pr, v1pr, k1, v2pr, k2);
      flint_mpn_mul(v3mr, v1mr, k1, v2mr, k2);
   }
   else
   {
      /* squaring version */

      /* evaluate B^(n1-1) * f1e(1/B^2) and B^(n1-2) * f1o(1/B^2) */
      _nmod_poly_KS2_pack(v1er, op1 + 2*(n1e - 1), n1e, -2, 2 * b, a1, k1);
      _nmod_poly_KS2_pack(v1or, op1 + 1 + 2*(n1o - 1), n1o, -2, 2 * b, b - a1, k1);

      /*
         compute  B^(n1-1) * f1(1/B) =
                      B^(n1-1) * f1e(1/B^2) + B^(n1-2) * f1o(1/B^2)
            and  |B^(n1-1) * f1(-1/B)| =
                     |B^(n1-1) * f1e(1/B^2) - B^(n1-2) * f1o(1/B^2)|
      */
      mpn_add_n(v1pr, v1er, v1or, k1);
      signed_mpn_sub_n(v1mr, v1er, v1or, k1);

      /*
         compute B^(n3-1) * h(1/B)  = (B^(n1-1) * f1(1/B))^2
             and B^(n3-1) * h(-1/B) = |B^(n1-1) * f1(-1/B)|^2
         hr_neg is cleared since h(-1/B) is never negative
      */
      flint_mpn_sqr(v3pr, v1pr, k1);
      flint_mpn_sqr(v3mr, v1mr, k1);
      v3m_neg = 0;
   }

   /*
      compute 2 * B^(n3-1) * he(1/B^2)
                     = B^(n3-1) * h(1/B) + B^(n3-1) * h(-1/B)
         and  2 * B^(n3-2) * ho(1/B^2)
                     = B^(n3-1) * h(1/B) - B^(n3-1) * h(-1/B)
   */
   if (v3m_neg)
   {
      mpn_sub_n(v3er, v3pr, v3mr, k3);
      mpn_add_n(v3or, v3pr, v3mr, k3);
   }
   else
   {
      mpn_add_n (v3er, v3pr, v3mr, k3);
      mpn_sub_n (v3or, v3pr, v3mr, k3);
   }

   /* -------------------------------------------------------------------------
          combine "normal" and "reciprocal" information
   */

   /* decompose he(B^2) and B^(2*(n3e-1)) * he(1/B^2) into base-B^2 digits */
   _nmod_poly_KS2_unpack(zn, v3en, n3e + 1, 2 * b, 1);
   _nmod_poly_KS2_unpack(zr, v3er, n3e + 1, 2 * b, a3 + 1);

   /* combine he(B^2) and he(1/B^2) information to get even coefficients of h */
   _nmod_poly_KS2_recover_reduce(res, 2, zn, zr, n3e, 2 * b, mod);

   /* decompose ho(B^2) and B^(2*(n3o-1)) * ho(1/B^2) into base-B^2 digits */
   _nmod_poly_KS2_unpack(zn, v3on, n3o + 1, 2 * b, b + 1);
   _nmod_poly_KS2_unpack(zr, v3or, n3o + 1, 2 * b, b - a3 + 1);

   /* combine ho(B^2) and ho(1/B^2) information to get odd coefficients of h */
   _nmod_poly_KS2_recover_reduce(res + 1, 2, zn, zr, n3o, 2 * b, mod);

   TMP_END;
}

void
nmod_poly_mul_KS4(nmod_poly_t res,
                 const nmod_poly_t poly1, const nmod_poly_t poly2)
{
    slong len_out;

    if ((poly1->length == 0) || (poly2->length == 0))
    {
        nmod_poly_zero(res);
        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        nmod_poly_t temp;
        nmod_poly_init2_preinv(temp, poly1->mod.n, poly1->mod.ninv, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS4(temp->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS4(temp->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
        nmod_poly_swap(res, temp);
        nmod_poly_clear(temp);
    }
    else
    {
        nmod_poly_fit_length(res, len_out);
        if (poly1->length >= poly2->length)
            _nmod_poly_mul_KS4(res->coeffs, poly1->coeffs, poly1->length,
                              poly2->coeffs, poly2->length,
                              poly1->mod);
        else
            _nmod_poly_mul_KS4(res->coeffs, poly2->coeffs, poly2->length,
                              poly1->coeffs, poly1->length,
                              poly1->mod);
    }

    res->length = len_out;
    _nmod_poly_normalise(res);
}
