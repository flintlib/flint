/*
    Copyright (C) 2008-2011 William Hart
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint-mparam.h"
#include "mpn_extras.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fft.h"

void _fmpz_poly_mulmid_SS(fmpz * res, const fmpz * poly1, slong len1,
               const fmpz * poly2, slong len2, slong nlo, slong nhi)
{
    slong len_out, loglen, loglen2, n;
    slong res_bits, limbs, size, i;
    ulong * ptr, ** t1, ** t2, ** tt, ** s1, ** ii, ** jj;
    slong bits1, bits2;
    ulong size1, size2;
    int sign = 0;
    int N;
    TMP_INIT;
    TMP_START;

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

    /* The algorithm requires trunc >= 3 */
    if (nhi <= 2)
    {
        if (res == poly1 || res == poly2)
        {
            fmpz * t = _fmpz_vec_init(nhi - nlo);
            _fmpz_poly_mulmid_classical(t, poly1, len1, poly2, len2, nlo, nhi);
            _fmpz_vec_swap(res, t, nhi - nlo);
            _fmpz_vec_clear(t, nhi - nlo);
        }
        else
        {
            _fmpz_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
        }
        return;
    }

    if (len1 < len2)
    {
        FLINT_SWAP(const fmpz *, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    slong trunc = nhi;

    len1 = FLINT_MIN(len1, trunc);
    len2 = FLINT_MIN(len2, trunc);

    len_out = len1 + len2 - 1;
    loglen  = FLINT_CLOG2(len_out);
    loglen2 = FLINT_CLOG2(len2);
    n = (WORD(1) << (loglen - 2));

    bits1 = _fmpz_vec_max_bits(poly1, len1);

    if (poly1 == poly2 && len1 == len2)
        bits2 = bits1;
    else
        bits2 = _fmpz_vec_max_bits(poly2, len2);

    size1 = (FLINT_ABS(bits1) + FLINT_BITS - 1) / FLINT_BITS;
    size2 = (FLINT_ABS(bits2) + FLINT_BITS - 1) / FLINT_BITS;

    /* Start with an upper bound on the number of bits needed */
    res_bits = FLINT_BITS * (size1 + size2) + loglen2 + 1;

    /* round up for sqrt2 trick */
    res_bits = (((res_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (res_bits - 1) / FLINT_BITS + 1; /* initial size of FFT coeffs */
    if (limbs > FFT_MULMOD_2EXPP1_CUTOFF) /* can't be worse than next power of 2 limbs */
        limbs = (WORD(1) << FLINT_CLOG2(limbs));
    size = limbs + 1;

    /* allocate space for ffts */

    N = flint_get_num_threads();
    ii = flint_malloc((4*(n + n*size) + 5*size*N)*sizeof(ulong));
    for (i = 0, ptr = (ulong *) ii + 4*n; i < 4*n; i++, ptr += size)
        ii[i] = ptr;
   t1 = TMP_ALLOC(N*sizeof(ulong *));
   t2 = TMP_ALLOC(N*sizeof(ulong *));
   s1 = TMP_ALLOC(N*sizeof(ulong *));
   tt = TMP_ALLOC(N*sizeof(ulong *));

   t1[0] = ptr;
   t2[0] = t1[0] + size*N;
   s1[0] = t2[0] + size*N;
   tt[0] = s1[0] + size*N;

   for (i = 1; i < N; i++)
   {
      t1[i] = t1[i - 1] + size;
      t2[i] = t2[i - 1] + size;
      s1[i] = s1[i - 1] + size;
      tt[i] = tt[i - 1] + 2*size;
   }

    if (poly1 != poly2)
    {
        jj = flint_malloc(4*(n + n*size)*sizeof(ulong));
        for (i = 0, ptr = (ulong *) jj + 4*n; i < 4*n; i++, ptr += size)
            jj[i] = ptr;
    } else jj = ii;

    /* put coefficients into FFT vecs */
    _fmpz_vec_get_fft(ii, poly1, limbs, len1);
    for (i = len1; i < 4*n; i++)
        flint_mpn_zero(ii[i], limbs + 1);

    if (poly1 != poly2 || len1 != len2)
    {
        _fmpz_vec_get_fft(jj, poly2, limbs, len2);
        for (i = len2; i < 4*n; i++)
            flint_mpn_zero(jj[i], limbs + 1);
    }

    if (bits1 < WORD(0) || bits2 < WORD(0))
    {
        sign = 1;
        bits1 = FLINT_ABS(bits1);
        bits2 = FLINT_ABS(bits2);
    }

    /* Recompute the number of bits/limbs now that we know how large everything is */
    res_bits = bits1 + bits2 + loglen2 + sign;

    /* round up res bits for sqrt2 */
    res_bits = (((res_bits - 1) >> (loglen - 2)) + 1) << (loglen - 2);

    limbs = (res_bits - 1) / FLINT_BITS + 1;
    limbs = fft_adjust_limbs(limbs); /* round up limbs for Nussbaumer */

    fft_convolution(ii, jj, loglen - 2, limbs, len_out, t1, t2, s1, tt);

    _fmpz_vec_set_fft(res, trunc - nlo, ii + nlo, limbs, sign); /* write res */

    flint_free(ii);
    if (poly1 != poly2 || len1 != len2)
        flint_free(jj);

    TMP_END;
}

void
fmpz_poly_mulmid_SS(fmpz_poly_t res, const fmpz_poly_t poly1,
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
    fmpz_poly_fit_length(res, len);
    _fmpz_poly_mulmid_SS(res->coeffs, poly1->coeffs, poly1->length,
                                poly2->coeffs, poly2->length, nlo, nhi);
    _fmpz_poly_set_length(res, len);
    _fmpz_poly_normalise(res);
}
