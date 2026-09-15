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
    /* input-major accumulation with a contiguous, vectorizable inner
       loop over any window (all values fit one small fmpz) */
    slong i, j, jlo, jhi, c;

    if (poly1 == poly2 && len1 == len2 && nlo == 0)
    {
        /* the original _fmpz_poly_sqr_tiny1 / _fmpz_poly_sqrlow_tiny1
           loops, verbatim: the loop shape matters for vectorization */
        if (nhi == 2 * len1 - 1)
        {
            _fmpz_vec_zero(res, 2 * len1 - 1);
            for (i = 0; i < len1; i++)
            {
                c = poly1[i];
                if (c != 0)
                {
                    res[2 * i] += c * c;
                    c *= 2;
                    for (j = i + 1; j < len1; j++)
                        res[i + j] += poly1[j] * c;
                }
            }
        }
        else
        {
            slong n = nhi;
            _fmpz_vec_zero(res, n);
            for (i = 0; i < len1; i++)
            {
                c = poly1[i];
                if (c != 0)
                {
                    if (2 * i < n)
                        res[2 * i] += c * c;
                    c *= 2;
                    for (j = i + 1; j < FLINT_MIN(len1, n - i); j++)
                        res[i + j] += poly1[j] * c;
                }
            }
        }
        return;
    }

    _fmpz_vec_zero(res, nhi - nlo);

    if (nlo == 0)
    {
        /* the original _fmpz_poly_mullow_tiny1 loop, verbatim */
        for (i = 0; i < len1; i++)
        {
            c = poly1[i];
            if (c != 0)
            {
                for (j = 0; j < FLINT_MIN(len2, nhi - i); j++)
                    res[i + j] += c * poly2[j];
            }
        }
        return;
    }

    for (i = 0; i < len1; i++)
    {
        c = poly1[i];

        if (c != 0)
        {
            jlo = FLINT_MAX(0, nlo - i);
            jhi = FLINT_MIN(len2, nhi - i);

            for (j = jlo; j < jhi; j++)
                res[i + j - nlo] += c * poly2[j];
        }
    }
}

static void
_fmpz_poly_mulmid_tiny2(fmpz * res, const fmpz * poly1, slong len1,
                                        const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    if (nlo == 0 && nhi == len1 + len2 - 1 && poly1 != poly2)
    {
        /* full product: input-major two-word accumulation, as in the
           original _fmpz_poly_mul_tiny2 */
        slong i, j, k, c, d;
        ulong hi, lo;
        nn_ptr tmp;
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(2 * (len1 + len2 - 1) * sizeof(ulong));
        flint_mpn_zero(tmp, 2 * (len1 + len2 - 1));

        for (i = 0; i < len1; i++)
        {
            c = poly1[i];

            if (c != 0)
            {
                for (j = 0; j < len2; j++)
                {
                    k = i + j;
                    d = poly2[j];

                    if (d != 0)
                    {
                        smul_ppmm(hi, lo, c, d);
                        add_ssaaaa(tmp[2 * k + 1], tmp[2 * k],
                                   tmp[2 * k + 1], tmp[2 * k], hi, lo);
                    }
                }
            }
        }

        for (i = 0; i < len1 + len2 - 1; i++)
        {
            lo = tmp[2 * i];
            hi = tmp[2 * i + 1];

            if (((slong) hi) >= 0)
            {
                fmpz_set_uiui(res + i, hi, lo);
            }
            else
            {
                sub_ddmmss(hi, lo, 0, 0, hi, lo);
                fmpz_neg_uiui(res + i, hi, lo);
            }
        }
        TMP_END;
        return;
    }

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
    slong len, bits1, bits2, sbits1, sbits2, rbits;

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

    sbits1 = _fmpz_vec_max_bits(poly1, len1);
    sbits2 = (poly1 == poly2 && len1 == len2) ? sbits1 : _fmpz_vec_max_bits(poly2, len2);
    bits1 = FLINT_ABS(sbits1);
    bits2 = FLINT_ABS(sbits2);

#if FLINT_HAVE_FFT_SMALL
    if (((nlo == 0 && nhi < len1 + len2 - 1) ?
            (len2 >= 100 && (bits1 + bits2 <= 40 || bits1 + bits2 >= 128 || len2 >= 130)) :
            ((len2 >= 80 && (bits1 + bits2 <= 40 || bits1 + bits2 >= 128 || len2 >= 100)) ||
             /* tiny coefficients on unbalanced full products beyond the
                tiny1 admission (len2 >= 40 + bsum/2 and len1 >= 70 +
                bsum/2): fft_small/KS rather than the packed mpn path */
             (bits1 + bits2 <= 40 && len2 >= 60 && len1 >= 70 + (bits1 + bits2) / 2))) &&
        !(FLINT_MAX(bits1, bits2) >= 8 * FLINT_MIN(bits1, bits2) &&
          FLINT_MIN(bits1, bits2) <= 60 &&
          FLINT_MAX(bits1, bits2) >= 384 && FLINT_MAX(bits1, bits2) <= 4200))
        if (_fmpz_poly_mul_mid_default_mpn_ctx(res, nlo, nhi, poly1, len1, poly2, len2))
            return;
#endif

    if (bits1 <= SMALL_FMPZ_BITCOUNT_MAX && bits2 <= SMALL_FMPZ_BITCOUNT_MAX &&
        ((nlo == 0 && nhi < len1 + len2 - 1) ?
            (len2 < 50 || (4 * len2 >= 3 * nhi && nhi < 150 + bits1 + bits2)) :
            (len2 < 40 + (bits1 + bits2) / 2 || len1 < 70 + (bits1 + bits2) / 2)))
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

    /* Very short operands with moderate coefficients: classical
       multiplication is the best algorithm throughout (up to 3 terms
       from 64 bits, up to 5 terms from 300 bits, measured on balanced
       and unbalanced shapes) and the remaining selection logic would
       cost a visible fraction of the product. */
    if (len1 <= 5 && FLINT_MAX(bits1, bits2) <= 900 &&
        (len1 <= 3 || FLINT_MAX(bits1, bits2) >= 300))
    {
        _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
        return;
    }

    /* Algorithm selection for all shapes (this is the single dispatcher
       behind _fmpz_poly_mul, _fmpz_poly_mullow, _fmpz_poly_mulhigh,
       _fmpz_poly_sqr and _fmpz_poly_sqrlow). bsum = bits1 + bits2 bounds
       the product coefficients. The boundaries were measured for
       uniformly sized random coefficients on balanced products, low
       and high products (n x n -> n) and middle products (2n x n -> n).
       Relative to a full product, the classical algorithm goes further
       on partial windows (its cost follows the number of terms in the
       window, while the subquadratic methods roughly pay for the full
       product), the packed-limb Karatsuba (with the Karatsuba middle
       product) is strongest on middle products, the Toom-Cook
       evaluation/interpolation gains less on low/high windows and
       nothing on low squarings, and the fft_small-based classical
       multiplication wins low/high windows with large coefficients. */
    {
        slong bsum = bits1 + bits2;
        slong bmin = FLINT_MIN(bits1, bits2);
        slong bmax = FLINT_MAX(bits1, bits2);
        int sparse = 0;

        /* Very unbalanced coefficient sizes between the operands, e.g.
           one operand with word-size coefficients: the evaluation/
           interpolation and packing based algorithms pad the small
           coefficients to the product size, while the direct products
           only pay for the actual sizes. */
        if ((bmin <= 320 || len1 > 2 * len2) && bmax >= 8 * bmin &&
            (bmax >= 512 || (bmin <= 64 && bmax >= 384)))
        {
            slong imb_mpn_max;

            if (len2 <= ((bmin <= 64) ? 30 : (bmin <= 128) ? 16 : 4) && (bmin <= 128 || len2 > 2 || len1 > 2 * len2))
            {
                _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
                return;
            }

            /* with word-size coefficients in the short operand the mpn
               path retains single-limb internal packing (leaving room
               for signs and recursive Karatsuba rounds) and stays ahead
               of the packing based methods up to much larger lengths */
            if (bmin <= 60 && bmax >= 700 && bmax <= 4200)
                imb_mpn_max = 208;
            else if (bmin <= 60 && bmax >= 384 && bmax < 700)
                imb_mpn_max = 144;
            else
                imb_mpn_max = 80;

            if (len2 < imb_mpn_max)
            {
                _fmpz_poly_mulmid_mpn_bits(res, poly1, len1, sbits1, poly2, len2, sbits2, nlo, nhi);
                return;
            }
        }
        else if (bsum >= 800 && len2 <= 32)
        {
            /* Coefficients of varying sizes within the operands hurt the
               evaluation/interpolation and packing based methods the same
               way (they pay for the largest coefficient), while the direct
               products pay for the actual sizes: discount the bit bounds
               to the average sizes for the choices below, and let the
               classical algorithm go further. */
            slong w1 = _fmpz_vec_weight_bits(poly1, len1);
            slong w2 = (poly1 == poly2 && len1 == len2) ? w1 : _fmpz_vec_weight_bits(poly2, len2);
            slong avg = w1 / len1 + w2 / len2;

            if (2 * avg < bsum)
            {
                bsum = FLINT_MAX(avg, bsum / 4);
                sparse = 1;
            }
        }
        slong nfull = len1 + len2 - 1;
        int full = (nlo == 0 && nhi == nfull);
        int edge = !full && (nlo == 0 || nhi == nfull);   /* low or high product */
        int squaring = (poly1 == poly2 && len1 == len2);
        slong cl_max, mpn_max, toom_max;

        if (bsum >= 3000)
        {
            if (squaring && !full)
                toom_max = (bsum <= 12000) ? 15 : 0;
            else if (full)
                toom_max = (bsum <= 4200) ? 16 : (bsum < 8000) ? 12 : (bsum < 16000) ? 10 : (bsum < 32000) ? 8 : 6;
            else if (edge)
                toom_max = (bsum <= 4200) ? 18 : (bsum <= 12000) ? 11 : (bsum < 20000) ? 10 : 4;
            else
                toom_max = (bsum <= 4200) ? 16 : (bsum < 16000) ? 10 : 8;

            /* for squaring, Kronecker substitution beats Toom in a band
               of sizes thanks to GMP's integer squaring (reflecting GMP's
               thresholds; may need re-measuring) */
            if (squaring && bsum >= 60000 && bsum < 200000 && len2 >= 4)
                toom_max = 0;

            /* (toom beats classical from 2 x 2 up at these sizes, except
               on low squarings) */
            if (squaring && !full)
                cl_max = (bsum <= 8000) ? 15 : 4;
            else if (!full)
                cl_max = (bsum <= 8000) ? 8 : 1;
            else
                cl_max = 1;

#if FLINT_HAVE_FFT_SMALL
            /* at very large coefficient sizes, classical multiplication
               with fft_small beats the direct and Toom paths on
               truncated products even at the smallest lengths */
            if (!full && bsum >= 40000 && len2 <= 4)
                if (_fmpz_poly_mulmid_classical_fft_small(res, poly1, len1, poly2, len2, nlo, nhi))
                    return;
#endif

            if (len2 <= cl_max)
            {
                _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
                return;
            }

            /* for full products with very large coefficients and strongly
               unbalanced lengths, a single Kronecker substitution
               multiplication amortizes better than blockwise Toom
               (measured: KS wins by 1.4-1.7x for len2 ~ 12 at ~5000 bits
               per coefficient once len1 > 2 len2, as in Hensel lifting) */
            if (len2 <= toom_max && (full ? (len1 <= 2 * len2 || bsum < 8000) : 1))
            {
                _fmpz_poly_mulmid_toom_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi);
                return;
            }

#if FLINT_HAVE_FFT_SMALL
            /* classical multiplication with fft_small: wins over the
               packing based methods on truncated products with large
               coefficients across a wide length range (it only computes
               the window), and on full products at very large sizes for
               moderate lengths just past the Toom range */
            if (!full && bsum >= 16000)
            {
                slong gmax;

                if (squaring)
                    gmax = (bsum < 80000) ? 61 : 87;
                else
                    gmax = (bsum < 40000) ? 21 : (bsum < 80000) ? 43 : 61;

                if (len2 <= gmax)
                    if (_fmpz_poly_mulmid_classical_fft_small(res, poly1, len1, poly2, len2, nlo, nhi))
                        return;
            }

            if (full && bsum >= 200000 && len2 <= 15)
                if (_fmpz_poly_mulmid_classical_fft_small(res, poly1, len1, poly2, len2, nlo, nhi))
                    return;
#endif
        }
        else
        {
            if (bsum >= 1600)
            {
                /* 800-1500 bits per coefficient */
                if (full && !squaring && len2 <= 4 && len1 <= 2 * len2)
                {
                    _fmpz_poly_mulmid_toom_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi);
                    return;
                }
                if (full && squaring && len2 <= 2)
                {
                    _fmpz_poly_mulmid_toom_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi);
                    return;
                }
                if (!full && !edge && !squaring && len2 <= 8 && len2 > 6)
                {
                    _fmpz_poly_mulmid_toom_karatsuba(res, poly1, len1, poly2, len2, nlo, nhi);
                    return;
                }
                cl_max = full ? 4 : 13;
                mpn_max = full ? (squaring ? 32 : 20) : edge ? 40 : 48;
            }
            else
            {
                cl_max = (bsum <= 200) ? (full ? 2 : 3) : (bsum <= 300) ? (full ? 3 : 4) : (bsum <= 1200) ? (full ? 11 : 15) : (full ? 6 : 13);
                mpn_max = (bsum <= 700) ? 80 : (bsum <= 1200) ? 64 : (full ? (squaring ? 32 : 20) : edge ? 20 : 56);
            }

            if (sparse)
                cl_max *= 2;

            if (len2 <= cl_max && (len1 <= 2 * len2 || bsum > 200 || !full))
            {
                _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
                return;
            }

            if (len2 < mpn_max)
            {
                _fmpz_poly_mulmid_mpn_bits(res, poly1, len1, sbits1, poly2, len2, sbits2, nlo, nhi);
                return;
            }
        }

#if FLINT_HAVE_FFT_SMALL
        /* The fft_small-based KS is so efficient that SS currently
           only wins in a specific medium-size region and for
           huge products when using many threads. */
        if ((squaring && len2 >= 8 && len2 <= 75 && bsum >= 800 && bsum <= 4000) ||
            (len1 + len2 >= 5000 && bsum >= 5000 + (len1 + len2) / 10 && flint_get_num_threads() >= 4) ||
            (!squaring && bsum >= 80000 && (double) (len1 + len2) * bsum > 4e9))
            _fmpz_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
        else
            _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
#else
        {
            slong limbs1 = (bits1 + FLINT_BITS - 1) / FLINT_BITS;
            slong limbs2 = (bits2 + FLINT_BITS - 1) / FLINT_BITS;

            if (limbs1 + limbs2 <= 8)
                _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
            else if ((limbs1 + limbs2) / 2048 > len1 + len2)
                _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
            else if ((limbs1 + limbs2) * FLINT_BITS * 4 < len1 + len2)
                _fmpz_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
            else
                _fmpz_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
        }
#endif
    }
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

