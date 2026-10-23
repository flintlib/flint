/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

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
#include "fmpzi.h"

#define GUARD_LIMBS 3

/*
    Middle products of Gaussian integer polynomials by Kronecker
    substitution. The real and imaginary parts of each operand are packed
    into integers at a common field size, giving Gaussian integers
    A = Ar + i Ai and B = Br + i Bi; the real and imaginary parts of
    A B are then the packed real and imaginary parts of the polynomial
    product. One flint_mpn_mul_complex (or flint_mpn_sqr_complex) call
    does the work, so above its transform cutoff the four parts are
    transformed once and combined pointwise (4 forward and 2 inverse
    transforms; 2 and 2 for squaring), and below it the Karatsuba or
    schoolbook schemes apply as they would to any complex integer
    product -- including when a part is zero or much shorter than its
    partner.

    Windows are computed as two's-complement windows of the packed real
    and imaginary parts, from which signed fields unpack directly with
    no sign to track: flint_mpn_mullow_complex (flint_mpn_sqrlow_complex)
    when the window starts at the bottom and ends below the top (exact
    residues), and flint_mpn_mulmid_complex (flint_mpn_sqrmid_complex)
    when it starts above the bottom, with guard limbs below the first
    field. Above the transform cutoff the latter reconstructs only the
    CRT coefficients reaching the window and uses a cyclic transform
    when the window allows it, which is where middle products gain
    transform length.
*/

/* pack (poly, len) with the given field size into a zeroed array of
   at least 1 limb, returning the limb count and setting the sign
   (0 or 1) of the packed value */
static slong
_pack_part(nn_ptr * arr, int * sgn, const fmpz * poly, slong len, slong bits)
{
    slong limbs;

    while (len > 0 && fmpz_is_zero(poly + len - 1))
        len--;

    if (len == 0)
    {
        *arr = flint_calloc(1, sizeof(ulong));
        *sgn = 0;
        return 1;
    }

    *sgn = (fmpz_sgn(poly + len - 1) < 0);
    limbs = (bits * len - 1) / FLINT_BITS + 1;
    *arr = flint_calloc(limbs, sizeof(ulong));
    _fmpz_poly_bit_pack(*arr, poly, len, bits, *sgn ? -1 : 0);

    return limbs;
}

void
_fmpzi_poly_mulmid_KS(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    slong bits1, bits2, bits, i, n, zn, ln, zlo;
    slong arn, ain, brn, bin, zr_len, zi_len;
    int ar_sgn, ai_sgn, br_sgn, bi_sgn, squaring, low;
    nn_ptr ar, ai, br, bi, zr, zi;
    fmpz * t, * re, * im;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

    n = nhi - nlo;

    /* low truncation of the inputs */
    len1 = FLINT_MIN(len1, nhi);
    len2 = FLINT_MIN(len2, nhi);

    /* high truncation of the inputs */
    {
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
    }

    squaring = (poly1 == poly2 && len1 == len2);

    /* shallow copies of the parts, read only by the packing */
    t = flint_malloc(2 * (len1 + (squaring ? 0 : len2)) * sizeof(fmpz));
    for (i = 0; i < len1; i++)
    {
        t[i] = *fmpzi_realref(poly1 + i);
        t[len1 + i] = *fmpzi_imagref(poly1 + i);
    }
    if (!squaring)
    {
        for (i = 0; i < len2; i++)
        {
            t[2 * len1 + i] = *fmpzi_realref(poly2 + i);
            t[2 * len1 + len2 + i] = *fmpzi_imagref(poly2 + i);
        }
    }

    bits1 = _fmpz_vec_max_bits(t, 2 * len1);
    bits1 = FLINT_ABS(bits1);
    if (squaring)
        bits2 = bits1;
    else
    {
        bits2 = _fmpz_vec_max_bits(t + 2 * len1, 2 * len2);
        bits2 = FLINT_ABS(bits2);
    }

    if (bits1 == 0 || bits2 == 0)
    {
        flint_free(t);
        for (i = 0; i < n; i++)
            fmpzi_zero(res + i);
        return;
    }

    /* each part of a coefficient of the product is a sum of at most
       2 min(len1, len2) terms below 2^(bits1 + bits2) in absolute value;
       one more bit for the sign of the field */
    bits = bits1 + bits2 + FLINT_BIT_COUNT(FLINT_MIN(len1, len2)) + 2;

    /* a real or imaginary operand makes every part a single product,
       of min(len1, len2) terms */
    if (_fmpz_vec_is_zero(t, len1) || _fmpz_vec_is_zero(t + len1, len1) ||
        (!squaring && (_fmpz_vec_is_zero(t + 2 * len1, len2) ||
                       _fmpz_vec_is_zero(t + 2 * len1 + len2, len2))))
        bits--;

    arn = _pack_part(&ar, &ar_sgn, t, len1, bits);
    ain = _pack_part(&ai, &ai_sgn, t + len1, len1, bits);
    if (!squaring)
    {
        brn = _pack_part(&br, &br_sgn, t + 2 * len1, len2, bits);
        bin = _pack_part(&bi, &bi_sgn, t + 2 * len1 + len2, len2, bits);
    }
    else
    {
        br = bi = NULL;
        brn = bin = 0;
        br_sgn = bi_sgn = 0;
    }

    flint_free(t);

    low = (nhi < len1 + len2 - 1);

    /* limbs holding the fields below nhi */
    ln = (nhi * bits + FLINT_BITS - 1) / FLINT_BITS;

    /* For nlo > 0, only the limbs from just below bit nlo bits - 1 (the
       borrow bit the signed unpacking reads first) are needed. The
       windowed products are within e units of the two's-complement
       window at limb zlo, |e| < 4 (m + 1) B for operands of m limbs; the
       unpacking of fields from nlo up is unaffected as long as
       e B^zlo stays below the distance from the low part
       sum_{k < nlo} r_k 2^(k bits) to +-2^(nlo bits - 1), which the
       bound |r_k| <= 2 (2^loglen - 1) 2^(bits1 + bits2) puts at
       2^(nlo bits - loglen - 1) or more. Three guard limbs leave
       192 - loglen - 2 bits for e, as in _fmpz_poly_mulmid_KS. */
    zlo = 0;
    if (nlo > 0)
        zlo = FLINT_MAX(0, (nlo * bits - 1) / FLINT_BITS - GUARD_LIMBS);

    /* room for the product (or the residues) and for the unpacking,
       which reads up to one limb past the field containing bit
       nhi * bits */
    if (low || zlo > 0)
        zn = ln + 2;
    else
        zn = FLINT_MAX(arn, ain) + (squaring ? FLINT_MAX(arn, ain) : FLINT_MAX(brn, bin)) + 1;
    zn = FLINT_MAX(zn, (nhi * bits) / FLINT_BITS + 2);
    zr = flint_malloc(2 * zn * sizeof(ulong));
    zi = zr + zn;

    if (zlo > 0)
    {
        /* two's-complement windows [zlo, ln), placed at their offset
           (the unpacking reads nothing below limb zlo + GUARD_LIMBS) */
        if (squaring)
            flint_mpn_sqrmid_complex(zr + zlo, zi + zlo, zlo, ln,
                                     ar, arn, ar_sgn, ai, ain, ai_sgn);
        else
            flint_mpn_mulmid_complex(zr + zlo, zi + zlo, zlo, ln,
                                     ar, arn, ar_sgn, ai, ain, ai_sgn,
                                     br, brn, br_sgn, bi, bin, bi_sgn);

        zr_len = zi_len = ln;
    }
    else if (low)
    {
        if (squaring)
            flint_mpn_sqrlow_complex(zr, zi, ln, ar, arn, ar_sgn, ai, ain, ai_sgn);
        else
            flint_mpn_mullow_complex(zr, zi, ln, ar, arn, ar_sgn, ai, ain, ai_sgn,
                                     br, brn, br_sgn, bi, bin, bi_sgn);

        /* two's-complement residues: unpacked as they stand, the limbs
           above them being outside every field read */
        zr_len = zi_len = ln;
    }
    else if (squaring)
        flint_mpn_sqr_complex(zr, &zr_len, zi, &zi_len,
                              ar, arn, ar_sgn, ai, ain, ai_sgn);
    else
        flint_mpn_mul_complex(zr, &zr_len, zi, &zi_len,
                              ar, arn, ar_sgn, ai, ain, ai_sgn,
                              br, brn, br_sgn, bi, bin, bi_sgn);

    flint_free(ar);
    flint_free(ai);
    if (!squaring)
    {
        flint_free(br);
        flint_free(bi);
    }

    /* limbs above the significant length are unspecified */
    flint_mpn_zero(zr + FLINT_ABS(zr_len), zn - FLINT_ABS(zr_len));
    flint_mpn_zero(zi + FLINT_ABS(zi_len), zn - FLINT_ABS(zi_len));

    /* the inputs have been consumed, so res may alias them */
    re = _fmpz_vec_init(2 * n);
    im = re + n;

    _fmpz_poly_bit_unpack(re, nlo, nhi, zr, bits, zr_len < 0 ? -1 : 0);
    _fmpz_poly_bit_unpack(im, nlo, nhi, zi, bits, zi_len < 0 ? -1 : 0);

    for (i = 0; i < n; i++)
    {
        fmpz_swap(fmpzi_realref(res + i), re + i);
        fmpz_swap(fmpzi_imagref(res + i), im + i);
    }

    _fmpz_vec_clear(re, 2 * n);
    flint_free(zr);
}
