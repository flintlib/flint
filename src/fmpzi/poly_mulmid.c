/*
    Copyright (C) 2026 Fredrik Johansson

    Developed using Claude Opus 5.5

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpzi.h"
#include "fft.h"
#include "gr.h"
#include "gr_poly.h"

/*
    Algorithm selection for Gaussian integer polynomial windows.

    The boundaries were fitted to timings of every candidate on
    uniformly random coefficients (x86-64, one thread) for the shapes
    n x n -> 2n - 1, squares, n x n -> [0, n), n x n -> [n - 1, 2n - 1),
    the Newton middle product (2n - 1) x n -> [n - 1, 2n - 1) and
    4n x n -> 5n - 1, with real and imaginary parts of equal sizes, or
    a quarter, or zero. The variables are b, the average over the two
    operands of the larger part size in bits, and q, the number of
    coefficient products in the window per operand coefficient, which
    puts the crossovers of all the shapes at similar values (q is n / 2
    for n x n products, n / 3 for the Newton window, n / 4 for low and
    high halves and 4n / 5 for 4n x n).

      - At most 12 products below 6144 bits: the classical algorithm.
      - A vanishing or much smaller part: four real products, which skip
        the vanishing ones (a real times a complex polynomial is two).
      - From 8192 bits for windows, 16384 bits for full products and
        about 200000 bits for squares, up to a bound on q depending on b
        and on whether the
        window is the full product: the classical algorithm on fft_small
        transforms, which pays for the products themselves and not for
        the length of a Kronecker product.
      - Below a bound on q depending on b, on the shape and on squaring:
        real products, three (Karatsuba) from about 256 bits, or for
        squares from 6144 bits.
      - Otherwise Kronecker substitution (one complex integer product of
        the packed parts), except for Schoenhage-Strassen in the bands
        where it measured faster: squares at 384-3072 bits, short full
        products at 768-1536 bits, and middle products whose cyclic
        transform saves at least a quarter of the length, short ones at
        768-3072 bits and long ones from 24576 bits; and for huge
        products (sum of the part sizes from 32768 bits, more than
        2.5e8 times the sum of the lengths), for memory.
*/

/* size classes, upper bounds in bits */
static const slong _bits_class[8] = { 96, 192, 384, 768, 1536, 3072, 6144, 12288 };

/* q from which the transformed algorithms (KS or SS) beat the real
   products, per size class */
static const double _min_q_full[9] = { 30, 16, 16, 13, 5, 3.5, 3.5, 3, 2 };
static const double _min_q_window[9] = { 21, 21, 8, 8, 3.5, 2.7, 2.5, 2, 2 };
static const double _min_q_sqr[9] = { 1e300, 256, 22, 11, 5.5, 4, 5.5, 4, 2 };

static int
_size_class(slong b)
{
    int c = 0;

    while (c < 8 && b >= _bits_class[c])
        c++;

    return c;
}

/* products a_i b_j with nlo <= i + j < nhi (len1 >= len2) */
static double
_window_terms(slong len1, slong len2, slong nlo, slong nhi)
{
    double t = 0;
    slong k;

    for (k = nlo; k < nhi; k++)
        t += FLINT_MIN(k, len1 - 1) - FLINT_MAX(0, k - len2 + 1) + 1;

    return t;
}

/* q up to which the classical algorithm on fft_small transforms beats
   the transformed algorithms (squares count each pair once) */
static double
_classical_fft_small_q(slong b, int full, int squaring)
{
    if (squaring)
        return (b >= 196608) ? 4.0 : 0.0;

    if (b < 8192)
        return 0.0;

    if (full)
        return (b < 16384) ? 0.0 : (b < 65536) ? 3.0 : (b < 131072) ? 3.5 : (b < 262144) ? 6.5 : 8.0;
    else
        return (b < 16384) ? 3.0 : (b < 32768) ? 4.5 : (b < 65536) ? 5.0 : (b < 131072) ? 8.0 : 8.5;
}

/* does the cyclic Schoenhage-Strassen transform save a quarter or more
   of the truncated length? (lengths after the truncation of the inputs) */
static int
_SS_cyclic_saves(slong len1, slong len2, slong nlo, slong nhi)
{
    slong len_out = len1 + len2 - 1, trunc, depth;

    if (len_out < 4)
        return 0;

    depth = fft_convolution_choose_mid(&trunc, len_out, nlo, nhi);

    return trunc == (WORD(4) << depth) && trunc != len_out &&
           4 * trunc <= 3 * len_out;
}

static void
_fmpzi_poly_mulmid_reorder(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi, int karatsuba)
{
    gr_ctx_t ctx, zctx;

    gr_ctx_init_fmpzi(ctx);
    gr_ctx_init_fmpz(zctx);
    GR_MUST_SUCCEED(_gr_poly_mulmid_complex_reorder(res, poly1, len1, poly2, len2,
                                                    nlo, nhi, karatsuba, ctx, zctx));
    gr_ctx_clear(zctx);
    gr_ctx_clear(ctx);
}

static void
_fmpzi_poly_mulmid_classical(fmpzi_struct * res, const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2, slong nlo, slong nhi)
{
    gr_ctx_t ctx;

    gr_ctx_init_fmpzi(ctx);
    GR_MUST_SUCCEED(_gr_poly_mulmid_classical(res, poly1, len1, poly2, len2,
                                              nlo, nhi, ctx));
    gr_ctx_clear(ctx);
}

void
_fmpzi_poly_mulmid(fmpzi_struct * res,
    const fmpzi_struct * poly1, slong len1,
    const fmpzi_struct * poly2, slong len2,
    slong nlo, slong nhi)
{
    slong i, ra, ia, rb, ib, ma, mb, b;
    int squaring, full, c;
    double terms, q;

    FLINT_ASSERT(len1 != 0);
    FLINT_ASSERT(len2 != 0);
    FLINT_ASSERT(nlo >= 0);
    FLINT_ASSERT(nlo < nhi);
    FLINT_ASSERT(nhi <= len1 + len2 - 1);

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

    if (len1 < len2)
    {
        FLINT_SWAP(const fmpzi_struct *, poly1, poly2);
        FLINT_SWAP(slong, len1, len2);
    }

    squaring = (poly1 == poly2 && len1 == len2);

    if (len2 == 1)
    {
        fmpzi_t c0;
        fmpzi_init(c0);
        fmpzi_set(c0, poly2);
        for (i = nlo; i < nhi; i++)
            fmpzi_mul(res + i - nlo, poly1 + i, c0);
        fmpzi_clear(c0);
        return;
    }

    ra = ia = rb = ib = 0;
    for (i = 0; i < len1; i++)
    {
        slong br = fmpz_bits(fmpzi_realref(poly1 + i));
        slong bi = fmpz_bits(fmpzi_imagref(poly1 + i));
        ra = FLINT_MAX(ra, br);
        ia = FLINT_MAX(ia, bi);
    }
    if (squaring)
    {
        rb = ra;
        ib = ia;
    }
    else
    {
        for (i = 0; i < len2; i++)
        {
            slong br = fmpz_bits(fmpzi_realref(poly2 + i));
            slong bi = fmpz_bits(fmpzi_imagref(poly2 + i));
            rb = FLINT_MAX(rb, br);
            ib = FLINT_MAX(ib, bi);
        }
    }

    ma = FLINT_MAX(ra, ia);
    mb = FLINT_MAX(rb, ib);

    if (ma == 0 || mb == 0)
    {
        for (i = 0; i < nhi - nlo; i++)
            fmpzi_zero(res + i);
        return;
    }

    terms = _window_terms(len1, len2, nlo, nhi);
    b = (ma + mb) / 2;

    if (terms <= 12 && b < 6144)
    {
        _fmpzi_poly_mulmid_classical(res, poly1, len1, poly2, len2, nlo, nhi);
        return;
    }

    /* a vanishing or much smaller part: the real products, which skip
       vanishing parts, against packings sized by the larger parts */
    if (ra == 0 || ia == 0 || rb == 0 || ib == 0 ||
        2 * (ra + ia + rb + ib) < 3 * (ma + mb))
    {
        _fmpzi_poly_mulmid_reorder(res, poly1, len1, poly2, len2, nlo, nhi, 0);
        return;
    }

    c = _size_class(b);
    full = (nlo == 0 && nhi == len1 + len2 - 1);
    q = terms / (double) (len1 + len2);

    if ((squaring ? q / 2 : q) <= _classical_fft_small_q(b, full, squaring) &&
        _fmpzi_poly_mulmid_classical_fft_small(res, poly1, len1, poly2, len2, nlo, nhi))
        return;

    if (q < (squaring ? _min_q_sqr[c] : (full ? _min_q_full[c] : _min_q_window[c])))
    {
        int karatsuba = squaring ? (b >= 6144) : (b >= 384 || (b >= 192 && q >= 4));
        _fmpzi_poly_mulmid_reorder(res, poly1, len1, poly2, len2, nlo, nhi, karatsuba);
        return;
    }

    /* huge products: Schoenhage-Strassen within about 10% of the time
       (from 16384-bit parts) in a fifth of the memory of Kronecker
       substitution, whose fft_small product expands its operands about
       20 times. The bound keeps the peak memory of KS below about
       600 MB: the four sequential real products used before took half
       the memory of one Gaussian KS product of the same size, so a
       higher bound regresses against them */
    if (ma + mb >= 32768 && (double) (len1 + len2) * (ma + mb) > 2.5e8)
    {
        _fmpzi_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
        return;
    }

    if (squaring ? ((c == 3 && q <= 32) || c == 4 || c == 5)
                 : ((full && c == 4 && q <= 32) ||
                    (!full && nlo > 0 && _SS_cyclic_saves(len1, len2, nlo, nhi) &&
                     (((c == 4 || c == 5) && q <= 22) ||
                      (b >= 49152 && q >= 10) || (b >= 24576 && q >= 21)))))
    {
        _fmpzi_poly_mulmid_SS(res, poly1, len1, poly2, len2, nlo, nhi);
        return;
    }

    _fmpzi_poly_mulmid_KS(res, poly1, len1, poly2, len2, nlo, nhi);
}
