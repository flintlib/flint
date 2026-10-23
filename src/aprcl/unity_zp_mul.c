/*
    Copyright (C) 2015 Vladimir Glazachev
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mod.h"
#include "fmpz_mod_poly.h"
#include "fmpz_poly.h"
#include "aprcl.h"

/*
    Sets f to g h (or g^2 if h == g) reduced by the cyclotomic polynomial
    Phi_{p^k} and modulo n: the product on the integers is computed by
    _fmpz_poly_mul or _fmpz_poly_sqr directly from the coefficients of
    g and h into the temporary space t (which must hold
    glen + hlen - 1 values), the cyclotomic reduction is done there with
    additions and subtractions, using x^{p^k} = 1 and then
    Phi_{p^k}(x) = Phi_p(x^{p^{k-1}}), i.e.
    x^d = -(1 + x^q + ... + x^{(p-2) q}) with q = p^{k-1}, d = (p-1) q,
    and finally each of the d coefficients is reduced modulo n once.
    The inputs need not be reduced by Phi_{p^k}.
*/
void
_unity_zp_mul_reduce(unity_zp f, const unity_zp g, const unity_zp h, fmpz * t)
{
    slong glen = g->poly->length, hlen = h->poly->length, len, i, j;
    ulong p = f->p, q = n_pow(p, f->exp - 1), d = (p - 1) * q;

    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx), fmpz_mod_ctx_modulus(g->ctx)));
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx), fmpz_mod_ctx_modulus(h->ctx)));

    if (glen == 0 || hlen == 0)
    {
        fmpz_mod_poly_zero(f->poly, f->ctx);
        return;
    }

    if (g == h)
        _fmpz_poly_sqr(t, g->poly->coeffs, glen);
    else if (glen >= hlen)
        _fmpz_poly_mul(t, g->poly->coeffs, glen, h->poly->coeffs, hlen);
    else
        _fmpz_poly_mul(t, h->poly->coeffs, hlen, g->poly->coeffs, glen);

    len = glen + hlen - 1;

    /* first x^{p^k} = 1 (Phi_{p^k} divides x^{p^k} - 1), which leaves at
       most q coefficients for the (p - 1)-term reduction */
    for (i = len - 1; i >= (slong) (p * q); i--)
        fmpz_add(t + i - p * q, t + i - p * q, t + i);
    len = FLINT_MIN(len, (slong) (p * q));

    for (i = len - 1; i >= (slong) d; i--)
    {
        if (!fmpz_is_zero(t + i))
            for (j = 0; j < (slong) p - 1; j++)
                fmpz_sub(t + i - d + j * q, t + i - d + j * q, t + i);
    }

    len = FLINT_MIN(len, (slong) d);
    fmpz_mod_poly_fit_length(f->poly, len, f->ctx);
    for (i = 0; i < len; i++)
        fmpz_mod_set_fmpz(f->poly->coeffs + i, t + i, f->ctx);
    _fmpz_mod_poly_set_length(f->poly, len);
    _fmpz_mod_poly_normalise(f->poly);
}


/*
    Straight-line multiplication and squaring programs for the small
    prime powers used most heavily by APRCL, adapted from the original
    hand-optimized routines by Vladimir Glazachev: minimal numbers of
    fmpz multiplications via Karatsuba-style schemes exploiting the ring
    structure. Unlike the originals these read the operands' coefficients
    directly (no fmpz_mod_poly_get_coeff_fmpz copies; the operands must
    have full length deg = phi(p^k)), keep each unreduced output in the
    scratch slot where it is produced (or a dead slot, by an O(1)
    fmpz_swap), and reduce all outputs into f at the very end, so that f
    may alias g or h. The scratch requirement is that of the originals.
    They win over the generic _unity_zp_mul_reduce for moduli up to a
    few thousand bits, where the multiplication count dominates; above
    that the generic path gains from Toom-Cook.
*/

/* reduce the d outputs held in the scratch slots slot[i] into f
   (not inlined: keeps d opaque, which avoids a spurious GCC
   -Wstringop-overread on the unrolled stores) */
FLINT_STATIC_NOINLINE void
_zp_store(unity_zp f, const fmpz * s, const unsigned char * slot, slong d)
{
    slong i;

    fmpz_mod_poly_fit_length(f->poly, d, f->ctx);
    for (i = 0; i < d; i++)
        fmpz_mod_set_fmpz(f->poly->coeffs + i, s + slot[i], f->ctx);
    _fmpz_mod_poly_set_length(f->poly, d);
    _fmpz_mod_poly_normalise(f->poly);
}

static void
_zp_ar1(fmpz * c, const fmpz * a, const fmpz * b, fmpz * s)
{
    fmpz_mul(c + 0, a + 0, b + 0);   /* c0 = a0 * b0 */
    fmpz_mul(s + 13, a + 1, b + 1);   /* d1 = a1 * b1 */
    fmpz_mul(c + 4, a + 2, b + 2);   /* c4 = a2 * b2 */
    fmpz_add(s + 11, a + 0, a + 1);   /* m1 = a0 + a1 */
    fmpz_add(s + 12, b + 0, b + 1);   /* m2 = b0 + b1 */
    fmpz_mul(s + 15, s + 11, s + 12);   /* d3 = m1 * m2 */
    fmpz_add(s + 11, a + 0, a + 2);   /* m1 = a0 + a2 */
    fmpz_add(s + 12, b + 0, b + 2);   /* m2 = b0 + b2 */
    fmpz_mul(s + 16, s + 11, s + 12);   /* d4 = m1 * m2 */
    fmpz_add(s + 11, a + 1, a + 2);   /* m1 = a1 + a2 */
    fmpz_add(s + 12, b + 1, b + 2);   /* m2 = b1 + b2 */
    fmpz_mul(s + 17, s + 11, s + 12);   /* d5 = m1 * m2 */
    fmpz_add(s + 14, c + 0, s + 13);   /* d2 = c0 + d1 */
    fmpz_sub(c + 1, s + 15, s + 14);   /* c1 = d3 - d2 */
    fmpz_add(s + 14, s + 16, s + 13);   /* d2 = d4 + d1 */
    fmpz_add(s + 16, c + 0, c + 4);   /* d4 = c0 + c4 */
    fmpz_sub(c + 2, s + 14, s + 16);   /* c2 = d2 - d4 */
    fmpz_add(s + 14, s + 13, c + 4);   /* d2 = d1 + c4 */
    fmpz_sub(c + 3, s + 17, s + 14);   /* c3 = d5 - d2 */
}

static void
_zp_ar2(fmpz * c, const fmpz * a, const fmpz * b, fmpz * s)
{
    fmpz_mul(c + 0, a + 0, b + 0);   /* c0 = a0 * b0 */
    fmpz_mul(s + 20, a + 1, b + 1);   /* d1 = a1 * b1 */
    fmpz_mul(s + 21, a + 2, b + 2);   /* d2 = a2 * b2 */
    fmpz_mul(c + 6, a + 3, b + 3);   /* c6 = a3 * b3 */
    fmpz_add(s + 15, a + 0, a + 1);   /* m1 = a0 + a1 */
    fmpz_add(s + 16, b + 0, b + 1);   /* m2 = b0 + b1 */
    fmpz_mul(s + 22, s + 15, s + 16);   /* d3 = m1 * m2 */
    fmpz_add(s + 15, a + 0, a + 2);   /* m1 = a0 + a2 */
    fmpz_add(s + 16, b + 0, b + 2);   /* m2 = b0 + b2 */
    fmpz_mul(s + 23, s + 15, s + 16);   /* d4 = m1 * m2 */
    fmpz_add(s + 17, a + 2, a + 3);   /* m3 = a2 + a3 */
    fmpz_add(s + 18, b + 2, b + 3);   /* m4 = b2 + b3 */
    fmpz_mul(s + 24, s + 17, s + 18);   /* d5 = m3 * m4 */
    fmpz_add(s + 17, a + 1, a + 3);   /* m3 = a1 + a3 */
    fmpz_add(s + 18, b + 1, b + 3);   /* m4 = b1 + b3 */
    fmpz_mul(s + 25, s + 17, s + 18);   /* d6 = m3 * m4 */
    fmpz_add(s + 26, c + 0, s + 20);   /* d7 = c0 + d1 */
    fmpz_sub(c + 1, s + 22, s + 26);   /* c1 = d3 - d7 */
    fmpz_add(s + 26, c + 0, s + 21);   /* d7 = c0 + d2 */
    fmpz_add(s + 27, s + 20, s + 23);   /* d8 = d1 + d4 */
    fmpz_sub(c + 2, s + 27, s + 26);   /* c2 = d8 - d7 */
    fmpz_add(s + 19, s + 15, s + 17);   /* m5 = m1 + m3 */
    fmpz_add(s + 17, s + 16, s + 18);   /* m3 = m2 + m4 */
    fmpz_add(s + 26, s + 21, c + 6);   /* d7 = d2 + c6 */
    fmpz_sub(c + 5, s + 24, s + 26);   /* c5 = d5 - d7 */
    fmpz_mul(s + 26, s + 17, s + 19);   /* d7 = m3 * m5 */
    fmpz_add(s + 27, c + 1, c + 5);   /* d8 = c1 + c5 */
    fmpz_add(s + 28, s + 27, s + 25);   /* d9 = d8 + d6 */
    fmpz_add(s + 27, s + 28, s + 23);   /* d8 = d9 + d4 */
    fmpz_sub(c + 3, s + 26, s + 27);   /* c3 = d7 - d8 */
    fmpz_add(s + 26, s + 25, s + 21);   /* d7 = d6 + d2 */
    fmpz_add(s + 27, s + 20, c + 6);   /* d8 = d1 + c6 */
    fmpz_sub(c + 4, s + 26, s + 27);   /* c4 = d7 - d8 */
}

static void
_zp_mul3(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    fmpz_mul(s + 6, x + 0, y + 0);   /*  d1 = x0 * y0    */
    fmpz_mul(s + 7, x + 1, y + 1);   /*  d2 = x1 * y1    */
    fmpz_sub(s + 4, x + 0, x + 1);   /*  m1 = x0 - x1    */
    fmpz_sub(s + 5, y + 1, y + 0);   /*  m2 = y1 - y0    */
    fmpz_mul(s + 8, s + 4, s + 5);   /*  d3 = m1 * m2    */
    fmpz_add(s + 8, s + 8, s + 6);   /*  d3 = d3 + d1    */
    /* z1 = d3 mod n: output 1 stays in s + 8 */
    fmpz_sub(s + 0, s + 6, s + 7);   /*  x0 = d1 - d2    */
    /* z0 = x0 mod n: output 0 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[2]) {0, 8}, 2);
}

static void
_zp_sqr3(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 2, x + 0, x + 1);   /*  m1 = x0 - x1    */
    fmpz_add(s + 3, x + 0, x + 1);   /*  m2 = x0 + x1    */
    fmpz_mul(s + 4, s + 2, s + 3);   /*  d1 = m1 * m2    */
    fmpz_add(s + 3, s + 2, x + 0);   /*  m2 = m1 + m0    */
    fmpz_swap(s + 2, s + 4);   /*  y0 = d1 mod n   */
    fmpz_mul(s + 4, x + 1, s + 3);   /*  d1 = x1 * m2    */
    /* y1 = d1 mod n: output 1 stays in s + 4 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[2]) {2, 4}, 2);
}

static void
_zp_mul4(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    fmpz_add(s + 4, x + 0, x + 1);   /*  m1 = x0 + x1    */
    fmpz_add(s + 5, y + 0, y + 1);   /*  m2 = y0 + y1    */
    fmpz_sub(s + 6, y + 1, y + 0);   /*  m3 = y1 - y0    */
    fmpz_mul(s + 7, s + 4, y + 0);   /*  d1 = m1 * y0    */
    fmpz_mul(s + 8, s + 5, x + 1);   /*  d2 = m2 * x1    */
    fmpz_mul(s + 9, s + 6, x + 0);   /*  d3 = m3 * x0    */
    fmpz_sub(s + 0, s + 7, s + 8);   /*  t[0] = d1 - d2  */
    fmpz_swap(s + 4, s + 0);   /*  z0 = t[0] mod n */
    fmpz_add(s + 0, s + 7, s + 9);   /*  t[0] = d1 + d3  */
    /* z1 = t[0] mod n: output 1 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[2]) {4, 0}, 2);
}

static void
_zp_sqr4(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 2, x + 0, x + 1);   /*  m1 = x0 - x1    */
    fmpz_add(s + 3, x + 0, x + 1);   /*  m2 = x0 + x1    */
    fmpz_mul(s + 4, s + 2, s + 3);   /*  d1 = m1 * m2    */
    fmpz_add(s + 2, x + 0, x + 0);   /*  m1 = x0 + x0    */
    fmpz_swap(s + 3, s + 4);   /*  y0 = d1 mod n   */
    fmpz_mul(s + 4, s + 2, x + 1);   /*  d1 = m1 * x1    */
    /* y1 = d1 mod n: output 1 stays in s + 4 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[2]) {3, 4}, 2);
}

static void
_zp_mul5(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    fmpz_sub(s + 8, x + 1, x + 3);   /* m1 = x1 - x3  */
    fmpz_sub(s + 9, y + 1, y + 3);   /* m2 = y1 - y3  */
    fmpz_sub(s + 10, x + 2, x + 3);   /* m3 = x2 - x3  */
    fmpz_sub(s + 11, y + 3, y + 2);   /* m4 = y3 - y2  */
    fmpz_sub(s + 12, x + 0, x + 1);   /* m5 = x0 - x1  */
    fmpz_sub(s + 13, y + 1, y + 0);   /* m6 = y1 - y0  */
    fmpz_sub(s + 14, x + 0, x + 2);   /* m7 = x0 - x2  */
    fmpz_sub(s + 15, y + 2, y + 0);   /* m8 = y2 - y0  */
    fmpz_mul(s + 16, x + 0, y + 0);   /* d0 = x0 * y0  */
    fmpz_mul(s + 18, s + 8, s + 9);   /* d2 = m1 * m2  */
    fmpz_add(s + 17, s + 16, s + 18);   /* d1 = d0 + d2  */
    fmpz_mul(s + 18, s + 10, s + 11);   /* d2 = m3 * m4  */
    fmpz_mul(s + 19, s + 12, s + 13);   /* d3 = m5 * m6  */
    fmpz_mul(s + 20, s + 14, s + 15);   /* d4 = m7 * m8  */
    fmpz_mul(s + 21, x + 1, y + 1);   /* d5 = x1 * y1  */
    fmpz_mul(s + 22, x + 2, y + 2);   /* d6 = x2 * y2  */
    fmpz_mul(s + 23, x + 3, y + 3);   /* d7 = x3 * y3  */
    fmpz_add(s + 24, s + 17, s + 18);   /* d8 = d1 + d2  */
    fmpz_sub(s + 0, s + 24, s + 21);   /* x0 = d8 - d5  */
    fmpz_swap(s + 12, s + 0);   /* z0 = x0 mod n */
    fmpz_add(s + 24, s + 17, s + 19);   /* d8 = d1 + d3  */
    fmpz_sub(s + 0, s + 24, s + 22);   /* x0 = d8 - d6  */
    fmpz_swap(s + 13, s + 0);   /* z1 = x0 mod n */
    fmpz_add(s + 24, s + 17, s + 20);   /* d8 = d1 + d4  */
    fmpz_sub(s + 0, s + 24, s + 23);   /* x0 = d8 - d7  */
    fmpz_swap(s + 21, s + 0);   /* z2 = x0 mod n */
    fmpz_sub(s + 10, s + 8, s + 14);   /* m3 = m1 - m7  */
    fmpz_add(s + 11, s + 9, s + 15);   /* m4 = m2 + m8  */
    fmpz_mul(s + 17, s + 10, s + 11);   /* d1 = m3 * m4  */
    fmpz_add(s + 24, s + 16, s + 17);   /* d8 = d0 + d1  */
    fmpz_add(s + 23, s + 24, s + 18);   /* d7 = d8 + d2  */
    fmpz_add(s + 24, s + 23, s + 19);   /* d8 = d7 + d3  */
    fmpz_add(s + 0, s + 24, s + 20);   /* d7 = d8 + d4  */
    /* z3 = d7 mod n: output 3 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[4]) {12, 13, 21, 0}, 4);
}

static void
_zp_sqr5(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 4, x + 0, x + 2);   /*  m1 = x0 - x2    */
    fmpz_add(s + 5, x + 0, x + 2);   /*  m2 = x0 + x2    */
    fmpz_sub(s + 6, x + 2, x + 1);   /*  m3 = x2 - x1    */
    fmpz_sub(s + 7, x + 0, x + 3);   /*  m4 = x0 - x3    */
    fmpz_sub(s + 8, x + 1, x + 0);   /*  m5 = x1 - x0    */
    fmpz_sub(s + 9, x + 2, x + 3);   /*  m6 = x2 - x3    */
    fmpz_sub(s + 10, x + 1, x + 3);   /*  m7 = x1 - x3    */
    fmpz_add(s + 11, x + 3, x + 3);   /*  m8 = x3 + x3    */
    fmpz_mul(s + 12, s + 4, s + 5);   /*  d1 = m1 * m2    */
    fmpz_mul(s + 13, s + 6, s + 11);   /*  d2 = m3 * m8    */
    fmpz_add(s + 14, s + 12, s + 13);   /*  d3 = d1 + d2    */
    fmpz_swap(s + 16, s + 14);   /*  y0 = d3 mod n   */
    fmpz_add(s + 11, s + 8, s + 10);   /*  m8 = m5 + m7    */
    fmpz_mul(s + 13, s + 7, s + 11);   /*  d2 = m4 * m8    */
    fmpz_add(s + 15, s + 12, s + 13);   /*  d4 = d1 + d2    */
    /* y1 = d4 mod n: output 1 stays in s + 15 */
    fmpz_add(s + 6, s + 4, x + 0);   /*  m3 = m1 + x0    */
    fmpz_mul(s + 12, x + 2, s + 6);   /*  d1 = x2 * m3    */
    fmpz_sub(s + 5, s + 10, x + 3);   /*  m2 = m7 - x3    */
    fmpz_mul(s + 13, s + 5, x + 1);   /*  d2 = m2 * x1    */
    fmpz_add(s + 14, s + 12, s + 13);   /*  d3 = d1 + d2    */
    fmpz_swap(s + 4, s + 14);   /*  y2 = d3 mod n   */
    fmpz_add(s + 10, s + 9, s + 9);   /*  m7 = m6 + m6    */
    fmpz_mul(s + 13, s + 10, s + 8);   /*  d2 = m7 * m5    */
    fmpz_add(s + 14, s + 12, s + 13);   /*  d3 = d1 + d2    */
    /* y3 = d3 mod n: output 3 stays in s + 14 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[4]) {16, 15, 4, 14}, 4);
}

static void
_zp_mul7(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    _zp_ar1(s + 50, x + 0, y + 0, s);
    _zp_ar1(s + 56, x + 3, y + 3, s);
    fmpz_sub(s + 0, x + 0, x + 3);   /*  a0 = x0 - x3    */
    fmpz_sub(s + 1, x + 1, x + 4);   /*  a1 = x1 - x4    */
    fmpz_sub(s + 2, x + 2, x + 5);   /*  a2 = x2 - x5    */
    fmpz_sub(s + 3, y + 3, y + 0);   /*  b0 = y3 - y0    */
    fmpz_sub(s + 4, y + 4, y + 1);   /*  b1 = y4 - y1    */
    fmpz_sub(s + 5, y + 5, y + 2);   /*  b2 = y5 - y2    */
    _zp_ar1(s + 61, s + 0, s + 3, s);
    fmpz_add(s + 68, s + 56, s + 64);   /*  d18 = d6 + d14  */
    fmpz_add(s + 66, s + 68, s + 53);   /*  d16 = d18 + d3  */
    fmpz_add(s + 68, s + 57, s + 65);   /*  d18 = d7 + d15  */
    fmpz_add(s + 67, s + 68, s + 54);   /*  d17 = d18 + d4  */
    fmpz_add(s + 68, s + 53, s + 61);   /*  d18 = d3 + d11  */
    fmpz_add(s + 53, s + 68, s + 50);   /*  d3 = d18 + d0   */
    fmpz_add(s + 68, s + 54, s + 62);   /*  d18 = d4 + d12  */
    fmpz_add(s + 54, s + 68, s + 51);   /*  d4 = d18 + d1   */
    fmpz_add(s + 55, s + 52, s + 63);   /*  d5 = d2 + d13   */
    fmpz_add(s + 63, s + 53, s + 56);   /*  d13 = d3 + d6   */
    fmpz_add(s + 64, s + 54, s + 57);   /*  d14 = d4 + d7   */
    fmpz_add(s + 65, s + 55, s + 58);   /*  d15 = d5 + d8   */
    fmpz_add(s + 56, s + 66, s + 59);   /*  d6 = d16 + d9   */
    fmpz_add(s + 57, s + 67, s + 60);   /*  d7 = d17 + d10  */
    fmpz_add(s + 68, s + 50, s + 57);   /*  d18 = d10 + d7  */
    fmpz_sub(s + 0, s + 68, s + 56);   /*  a0 = d18 - d6   */
    fmpz_swap(s + 1, s + 0);   /*  z0 = a0 mod n   */
    fmpz_add(s + 68, s + 51, s + 58);   /*  d18 = d1 + d8   */
    fmpz_sub(s + 0, s + 68, s + 56);   /*  a0 = d18 - d6   */
    fmpz_swap(s + 2, s + 0);   /*  z1 = a0 mod n   */
    fmpz_add(s + 68, s + 52, s + 59);   /*  d18 = d2 + d9   */
    fmpz_sub(s + 0, s + 68, s + 56);   /*  a0 = d18 - d6   */
    fmpz_swap(s + 3, s + 0);   /*  z2 = a0 mod n   */
    fmpz_add(s + 68, s + 63, s + 60);   /*  d18 = d13 + d10 */
    fmpz_sub(s + 0, s + 68, s + 56);   /*  a0 = d18 - d6   */
    fmpz_swap(s + 4, s + 0);   /*  z3 = a0 mod n   */
    fmpz_sub(s + 0, s + 64, s + 56);   /*  a0 = d14 - d6   */
    fmpz_swap(s + 5, s + 0);   /*  z4 = a0 mod n   */
    fmpz_sub(s + 0, s + 65, s + 56);   /*  a0 = d15 - d6   */
    /* z5 = a0 mod n: output 5 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[6]) {1, 2, 3, 4, 5, 0}, 6);
}

static void
_zp_sqr7(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 6, x + 0, x + 1);   /*  m1 = x0 - x1    */
    fmpz_sub(s + 7, x + 1, x + 2);   /*  m2 = x1 - x2    */
    fmpz_sub(s + 8, x + 2, x + 3);   /*  m3 = x2 - x3    */
    fmpz_sub(s + 9, x + 3, x + 4);   /*  m4 = x3 - x4    */
    fmpz_sub(s + 10, x + 5, x + 4);   /*  m5 = x5 - x4    */
    fmpz_add(s + 11, s + 6, s + 7);   /*  m6 = x6 + x7    */
    fmpz_add(s + 12, s + 7, s + 8);   /*  m7 = m2 + m3    */
    fmpz_add(s + 13, s + 8, s + 9);   /*  m8 = m3 + m4    */
    fmpz_sub(s + 14, x + 3, x + 5);   /*  m9 = x3 - x5    */
    fmpz_add(s + 15, s + 8, s + 11);   /*  m10 = m3 + m6   */
    fmpz_add(s + 16, s + 9, s + 12);   /*  m11 = m4 + m7   */
    fmpz_add(s + 18, s + 11, s + 13);   /*  m13 = m6 + m8   */
    fmpz_add(s + 19, s + 12, s + 14);   /*  m14 = m7 + m9   */
    fmpz_add(s + 21, x + 0, x + 1);   /*  m16 = x0 + x1   */
    fmpz_add(s + 22, x + 0, s + 15);   /*  m17 = m8 + m2   */
    fmpz_mul(s + 24, x + 3, s + 22);   /*  d1 = x3 * m17   */
    fmpz_sub(s + 22, s + 19, x + 4);   /*  m17 = m14 - m9  */
    fmpz_add(s + 23, s + 19, x + 4);   /*  m18 = m14 + m9  */
    fmpz_mul(s + 25, s + 22, s + 23);   /*  d2 = m17 * m18  */
    fmpz_sub(s + 22, s + 13, s + 7);   /*  m17 = m8 - m2   */
    fmpz_mul(s + 26, s + 16, s + 22);   /*  d3 = m11 * m17  */
    fmpz_add(s + 22, s + 19, s + 14);   /*  m17 = m14 + m9  */
    fmpz_mul(s + 27, s + 22, s + 12);   /*  d4 = m17 * m7   */
    fmpz_add(s + 22, x + 1, x + 1);   /*  m17 = x1 + x1   */
    fmpz_mul(s + 28, s + 22, s + 11);   /*  d5 = m17 * m6   */
    fmpz_mul(s + 29, s + 6, s + 21);   /*  d6 = m1 * m16   */
    fmpz_add(s + 22, s + 8, s + 8);   /*  m17 = m3 + m3   */
    fmpz_add(s + 7, x + 0, s + 18);   /*  m2 = x0 + m13   */
    fmpz_mul(s + 30, s + 22, s + 10);   /*  d7 = m17 * m5   */
    fmpz_add(s + 31, s + 24, s + 25);   /*  d8 = d1 + d2    */
    fmpz_add(s + 24, s + 31, s + 26);   /*  d1 = d8 + d5    */
    fmpz_swap(s + 11, s + 24);   /*  y3 = d1 mod n   */
    fmpz_add(s + 31, s + 26, s + 27);   /*  d8 = d3 + d4    */
    fmpz_add(s + 24, s + 31, s + 28);   /*  d1 = d8 + d5    */
    fmpz_swap(s + 17, s + 24);   /*  y1 = d1 mod n   */
    fmpz_add(s + 31, s + 27, s + 29);   /*  d8 = d4 + d6    */
    fmpz_add(s + 24, s + 31, s + 30);   /*  d1 = d8 + d7    */
    fmpz_swap(s + 18, s + 24);   /*  y0 = d1 mod n   */
    fmpz_add(s + 22, s + 12, s + 19);   /*  m17 = m7 + m14  */
    fmpz_mul(s + 24, s + 14, s + 22);   /*  d1 = m9 * m17   */
    fmpz_sub(s + 22, s + 13, x + 5);   /*  m17 = m8 - x5   */
    fmpz_add(s + 23, x + 2, s + 10);   /*  m18 = x2 + m5   */
    fmpz_mul(s + 25, s + 22, s + 23);   /*  d2 = m17 * m18  */
    fmpz_mul(s + 26, s + 7, x + 4);   /*  d3 = m2 * x4    */
    fmpz_add(s + 22, s + 8, s + 13);   /*  m17 = m3 + m8   */
    fmpz_mul(s + 27, s + 22, s + 9);   /*  d4 = m17 * m4   */
    fmpz_add(s + 22, s + 6, s + 6);   /*  m17 = m1 + m1   */
    fmpz_mul(s + 28, s + 22, s + 10);   /*  d5 = m17 * m5   */
    fmpz_sub(s + 22, s + 19, s + 10);   /*  m17 = m14 - m5  */
    fmpz_mul(s + 29, s + 22, s + 16);   /*  d6 = m17 * m11  */
    fmpz_add(s + 22, x + 2, x + 2);   /*  m17 = x2 + x2   */
    fmpz_mul(s + 30, s + 22, s + 15);   /*  d7 = m17 * m10  */
    fmpz_add(s + 31, s + 24, s + 25);   /*  d8 = d1 + d2    */
    fmpz_add(s + 24, s + 31, s + 26);   /*  d1 = d8 + d3    */
    fmpz_swap(s + 6, s + 24);   /*  y4 = d1 mod n   */
    fmpz_add(s + 31, s + 26, s + 27);   /*  d8 = d3 + d4    */
    fmpz_add(s + 24, s + 31, s + 28);   /*  d1 = d8 + d5    */
    fmpz_swap(s + 7, s + 24);   /*  y5 = d1 mod n   */
    fmpz_add(s + 31, s + 27, s + 29);   /*  d8 = d4 + d6    */
    fmpz_add(s + 24, s + 31, s + 30);   /*  d1 = d8 + d7    */
    /* y2 = d1 mod n: output 2 stays in s + 24 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[6]) {18, 17, 24, 11, 6, 7}, 6);
}

static void
_zp_mul8(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    fmpz_add(s + 8, x + 1, x + 3);   /*  m1 = x1 + x3    */
    fmpz_add(s + 9, y + 1, y + 3);   /*  m2 = y1 + y3    */
    fmpz_add(s + 10, x + 2, x + 3);   /*  m3 = x2 + x3    */
    fmpz_add(s + 11, y + 2, y + 3);   /*  m4 = y2 + y3    */
    fmpz_add(s + 12, x + 0, x + 1);   /*  m5 = x0 + x1    */
    fmpz_add(s + 13, y + 0, y + 1);   /*  m6 = y0 + y1    */
    fmpz_add(s + 14, x + 0, x + 2);   /*  m7 = x0 + x2    */
    fmpz_add(s + 15, y + 0, y + 2);   /*  m8 = y0 + y2    */
    fmpz_mul(s + 16, x + 0, y + 0);   /*  d0 = x0 * y0    */
    fmpz_mul(s + 17, x + 1, y + 1);   /*  d1 = x1 * y1    */
    fmpz_mul(s + 18, x + 2, y + 2);   /*  d2 = x2 * y2    */
    fmpz_mul(s + 19, x + 3, y + 3);   /*  d3 = x3 * y3    */
    fmpz_mul(s + 22, s + 12, s + 13);   /*  d6 = m5 * m6    */
    fmpz_mul(s + 23, s + 14, s + 15);   /*  d7 = m7 * m8    */
    fmpz_mul(s + 24, s + 8, s + 9);   /*  d8 = m1 * m2    */
    fmpz_mul(s + 25, s + 10, s + 11);   /*  d9 = m3 * m4    */
    fmpz_add(s + 10, s + 8, s + 14);   /*  m3 = m1 + m7    */
    fmpz_add(s + 11, s + 9, s + 15);   /*  m4 = m2 + m8    */
    fmpz_mul(s + 20, s + 10, s + 11);   /*  d4 = m3 * m4    */
    fmpz_add(s + 26, s + 16, s + 17);   /*  d10 = d0 + d1   */
    fmpz_add(s + 27, s + 18, s + 19);   /*  d11 = d2 + d3   */
    fmpz_add(s + 28, s + 26, s + 19);   /*  d12 = d10 + d3  */
    fmpz_add(s + 21, s + 24, s + 18);   /*  d5 = d8 + d2    */
    fmpz_sub(s + 0, s + 28, s + 21);   /*  t[0] = d12 - d5 */
    fmpz_swap(s + 8, s + 0);   /*  z0 = t[0] mod n */
    fmpz_add(s + 28, s + 22, s + 27);   /*  d12 = d6 + d11  */
    fmpz_add(s + 21, s + 26, s + 25);   /*  d5 = d10 + d9   */
    fmpz_sub(s + 0, s + 28, s + 21);   /*  t[0] = d12 - d5 */
    fmpz_swap(s + 9, s + 0);   /*  z1 = t[0] mod n */
    fmpz_add(s + 28, s + 17, s + 23);   /*  d12 = d1 + d7   */
    fmpz_add(s + 21, s + 16, s + 27);   /*  d5 = d0 + d11   */
    fmpz_sub(s + 0, s + 28, s + 21);   /*  t[0] = d12 - d5 */
    fmpz_swap(s + 10, s + 0);   /*  z2 = t[0] mod n */
    fmpz_add(s + 28, s + 23, s + 22);   /*  d12 = d7 + d6   */
    fmpz_add(s + 21, s + 28, s + 24);   /*  d5 = d12 + d8   */
    fmpz_add(s + 28, s + 21, s + 25);   /*  d12 = d5 + d8   */
    fmpz_add(s + 19, s + 26, s + 27);   /*  d3 = d10 + d11  */
    fmpz_add(s + 21, s + 19, s + 20);   /*  d5 = d3 + d4    */
    fmpz_sub(s + 0, s + 21, s + 28);   /*  t[0] = d5 - d12 */
    /* z3 = t[0] mod n: output 3 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[4]) {8, 9, 10, 0}, 4);
}

static void
_zp_sqr8(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 4, x + 0, x + 2);   /*  m1 = x0 - x2    */
    fmpz_add(s + 5, x + 0, x + 2);   /*  m2 = x0 + x2    */
    fmpz_sub(s + 6, x + 1, x + 3);   /*  m3 = x1 - x3    */
    fmpz_add(s + 7, x + 1, x + 3);   /*  m4 = x1 + x3    */
    fmpz_add(s + 8, x + 0, x + 0);   /*  m5 = x0 + x0    */
    fmpz_add(s + 9, x + 1, x + 1);   /*  m6 = x1 + x1    */
    fmpz_add(s + 10, s + 4, s + 6);   /*  m7 = m1 + m3    */
    fmpz_add(s + 11, s + 5, s + 7);   /*  m8 = m2 + m4    */
    fmpz_mul(s + 12, s + 4, s + 5);   /*  d1 = m1 * m2    */
    fmpz_mul(s + 13, s + 6, s + 7);   /*  d2 = m3 * m4    */
    fmpz_mul(s + 14, s + 9, x + 3);   /*  d3 = m6 * x3    */
    fmpz_mul(s + 15, s + 8, x + 2);   /*  d4 = m5 * x2    */
    fmpz_add(s + 5, x + 2, x + 3);   /*  m2 = x2 + x3    */
    fmpz_sub(s + 16, s + 12, s + 14);   /*  d5 = d1 - d3    */
    fmpz_swap(s + 6, s + 16);   /*  y0 = d5 mod n   */
    fmpz_add(s + 17, s + 13, s + 15);   /*  d6 = d2 + d4    */
    fmpz_swap(s + 7, s + 17);   /*  y2 = d6 mod n   */
    fmpz_mul(s + 16, s + 10, s + 11);   /*  d5 = m7 * m8    */
    fmpz_add(s + 17, s + 12, s + 13);   /*  d6 = d1 + d2    */
    fmpz_sub(s + 13, s + 16, s + 17);   /*  d2 = d5 - d6    */
    fmpz_swap(s + 10, s + 13);   /*  y1 = d2 mod n   */
    fmpz_add(s + 4, s + 8, s + 9);   /*  m1 = m5 + m6    */
    fmpz_mul(s + 12, s + 4, s + 5);   /*  d1 = m1 * m2    */
    fmpz_add(s + 17, s + 14, s + 15);   /*  d6 = d3 + d4    */
    fmpz_sub(s + 13, s + 12, s + 17);   /*  d2 = d1 - d6    */
    /* y3 = d2 mod n: output 3 stays in s + 13 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[4]) {6, 10, 7, 13}, 4);
}

static void
_zp_mul9(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    _zp_ar1(s + 32, x + 0, y + 0, s);
    _zp_ar1(s + 38, x + 3, y + 3, s);
    fmpz_sub(s + 0, x + 0, x + 3);   /*  a0 = x0 - x3    */
    fmpz_sub(s + 1, x + 1, x + 4);   /*  a1 = x1 - x4    */
    fmpz_sub(s + 2, x + 2, x + 5);   /*  a2 = x2 - x5    */
    fmpz_sub(s + 3, y + 3, y + 0);   /*  b0 = y3 - y0    */
    fmpz_sub(s + 4, y + 4, y + 1);   /*  b1 = y4 - y1    */
    fmpz_sub(s + 5, y + 5, y + 2);   /*  b2 = y5 - y2    */
    _zp_ar1(s + 43, s + 0, s + 3, s);
    fmpz_add(s + 50, s + 38, s + 46);   /*  d18 = d6 + d14  */
    fmpz_add(s + 48, s + 50, s + 35);   /*  d16 = d18 + d3  */
    fmpz_add(s + 50, s + 39, s + 47);   /*  d18 = d7 + d15  */
    fmpz_add(s + 49, s + 50, s + 36);   /*  d17 = d18 + d4  */
    fmpz_add(s + 50, s + 35, s + 43);   /*  d18 = d3 + d11  */
    fmpz_add(s + 35, s + 50, s + 32);   /*  d3 = d18 + d0   */
    fmpz_add(s + 50, s + 36, s + 44);   /*  d18 = d4 + d12  */
    fmpz_add(s + 36, s + 50, s + 33);   /*  d4 = d18 + d1   */
    fmpz_add(s + 37, s + 34, s + 45);   /*  d5 = d2 + d13   */
    fmpz_sub(s + 0, s + 32, s + 48);   /*  a0 = d0 - d16   */
    fmpz_swap(s + 1, s + 0);   /*  z0 = a0 mod n   */
    fmpz_sub(s + 0, s + 33, s + 49);   /*  a0 = d1 - d17   */
    fmpz_swap(s + 2, s + 0);   /*  z1 = a0 mod n   */
    fmpz_sub(s + 0, s + 34, s + 40);   /*  a0 = d2 - d8    */
    fmpz_swap(s + 3, s + 0);   /*  z2 = a0 mod n   */
    /* z5 = d5 mod n: output 5 stays in s + 37 */
    fmpz_add(s + 50, s + 35, s + 38);   /*  d18 = d3 + d6   */
    fmpz_add(s + 51, s + 48, s + 41);   /*  d19 = d16 + d9  */
    fmpz_sub(s + 0, s + 50, s + 51);   /*  a0 = d18 - d19  */
    fmpz_swap(s + 4, s + 0);   /*  z3 = a0 mod n   */
    fmpz_add(s + 50, s + 36, s + 39);   /*  d18 = d4 + d7   */
    fmpz_add(s + 51, s + 42, s + 49);   /*  d19 = d10 + d17 */
    fmpz_sub(s + 0, s + 50, s + 51);   /*  a0 = d18 - d19  */
    /* z4 = a0 mod n: output 4 stays in s + 0 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[6]) {1, 2, 3, 4, 0, 37}, 6);
}

static void
_zp_sqr9(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_sub(s + 0, x + 0, x + 3);   /*  a0 = x0 - x3    */
    fmpz_sub(s + 1, x + 1, x + 4);   /*  a1 = x1 - x4    */
    fmpz_sub(s + 2, x + 2, x + 5);   /*  a2 = x2 - x5    */
    fmpz_add(s + 3, x + 0, x + 3);   /*  b0 = x0 + x3    */
    fmpz_add(s + 4, x + 1, x + 4);   /*  b1 = x1 + x4    */
    fmpz_add(s + 5, x + 2, x + 5);   /*  b2 = x2 + x5    */
    _zp_ar1(s + 26, s + 0, s + 3, s);
    fmpz_add(s + 3, x + 0, s + 0);   /*  b0 = x0 + a0    */
    fmpz_add(s + 4, x + 1, s + 1);   /*  b1 = x1 + a1    */
    fmpz_add(s + 5, x + 2, s + 2);   /*  b2 = x2 + a2    */
    _zp_ar1(s + 6, x + 3, s + 3, s);
    fmpz_sub(s + 0, s + 26, s + 9);   /*  a0 = d0 - c3    */
    fmpz_swap(s + 3, s + 0);   /*  y0 = a0 mod n   */
    fmpz_sub(s + 0, s + 27, s + 10);   /*  a0 = d1 - c4    */
    fmpz_swap(s + 4, s + 0);   /*  y1 = a0 mod n   */
    /* y2 = d2 mod n: output 2 stays in s + 28 */
    fmpz_add(s + 0, s + 29, s + 6);   /*  a0 = d3 + c0    */
    fmpz_sub(s + 1, s + 0, s + 9);   /*  a1 = a0 - c3    */
    fmpz_swap(s + 5, s + 1);   /*  y3 = a1 mod n   */
    fmpz_add(s + 0, s + 30, s + 7);   /*  a0 = d4 + c1    */
    fmpz_sub(s + 1, s + 0, s + 10);   /*  a1 = a0 - c4    */
    /* y4 = a1 mod n: output 4 stays in s + 1 */
    /* y5 = c2 mod n: output 5 stays in s + 8 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[6]) {3, 4, 28, 5, 1, 8}, 6);
}

static void
_zp_mul16(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    const fmpz * y = h->poly->coeffs;
    fmpz_add(s + 0, x + 0, x + 4);   /*  a0 = x0 + x4    */
    fmpz_add(s + 1, x + 1, x + 5);   /*  a1 = x1 + x5    */
    fmpz_add(s + 2, x + 2, x + 6);   /*  a2 = x2 + x6    */
    fmpz_add(s + 3, x + 3, x + 7);   /*  a3 = x3 + x7    */
    _zp_ar2(s + 50, s + 0, y + 0, s);
    fmpz_add(s + 0, y + 0, y + 4);   /*  a0 = y0 + y4    */
    fmpz_add(s + 1, y + 1, y + 5);   /*  a1 = y1 + y5    */
    fmpz_add(s + 2, y + 2, y + 6);   /*  a2 = y2 + y6    */
    fmpz_add(s + 3, y + 3, y + 7);   /*  a3 = y3 + y7    */
    _zp_ar2(s + 57, s + 0, x + 4, s);
    fmpz_sub(s + 0, y + 4, y + 0);   /*  a0 = y4 - y0    */
    fmpz_sub(s + 1, y + 5, y + 1);   /*  a1 = y5 - y1    */
    fmpz_sub(s + 2, y + 6, y + 2);   /*  a2 = y6 - y2    */
    fmpz_sub(s + 3, y + 7, y + 3);   /*  a3 = y7 - y3    */
    _zp_ar2(s + 8, s + 0, x + 0, s);
    fmpz_add(s + 1, s + 54, s + 57);   /*  a1 = d4 + d7    */
    fmpz_add(s + 2, s + 1, s + 12);   /*  a2 = a1 + c3    */
    fmpz_sub(s + 0, s + 50, s + 2);   /*  a0 = d0 - a2    */
    fmpz_swap(s + 3, s + 0);   /*  z0 = a0 mod n   */
    fmpz_add(s + 1, s + 55, s + 58);   /*  a1 = d5 + d8    */
    fmpz_add(s + 2, s + 1, s + 13);   /*  a2 = a1 + c4    */
    fmpz_sub(s + 0, s + 51, s + 2);   /*  a0 = d1 - a2    */
    fmpz_swap(s + 12, s + 0);   /*  z1 = a0 mod n   */
    fmpz_add(s + 1, s + 56, s + 59);   /*  a1 = d6 + d9    */
    fmpz_add(s + 2, s + 1, s + 14);   /*  a2 = a1 + c5    */
    fmpz_sub(s + 0, s + 52, s + 2);   /*  a0 = d2 - a2    */
    fmpz_swap(s + 13, s + 0);   /*  z2 = a0 mod n   */
    fmpz_sub(s + 0, s + 53, s + 60);   /*  a0 = d3 - d10   */
    fmpz_swap(s + 14, s + 0);   /*  z3 = a0 mod n   */
    fmpz_add(s + 1, s + 54, s + 50);   /*  a1 = d4 + d0    */
    fmpz_add(s + 2, s + 1, s + 8);   /*  a2 = a1 + c0    */
    fmpz_sub(s + 0, s + 2, s + 61);   /*  a0 = a2 - d11   */
    fmpz_swap(s + 8, s + 0);   /*  z4 = a0 mod n   */
    fmpz_add(s + 1, s + 55, s + 51);   /*  a1 = d5 + d1    */
    fmpz_add(s + 2, s + 1, s + 9);   /*  a2 = a1 + c0    */
    fmpz_sub(s + 0, s + 2, s + 62);   /*  a0 = a2 - d12   */
    fmpz_swap(s + 9, s + 0);   /*  z5 = a0 mod n   */
    fmpz_add(s + 1, s + 56, s + 52);   /*  a1 = d6 + d2    */
    fmpz_add(s + 2, s + 1, s + 10);   /*  a2 = a1 + c1    */
    fmpz_sub(s + 0, s + 2, s + 63);   /*  a0 = a2 - a13   */
    /* z6 = a0 mod n: output 6 stays in s + 0 */
    fmpz_add(s + 1, s + 53, s + 11);   /*  a1 = d3 + c2    */
    /* z7 = a1 mod n: output 7 stays in s + 1 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[8]) {3, 12, 13, 14, 8, 9, 0, 1}, 8);
}

static void
_zp_sqr16(unity_zp f, const unity_zp g, fmpz * s)
{
    const fmpz * x = g->poly->coeffs;
    fmpz_add(s + 0, x + 0, x + 4);   /*  a0 = x0 + x4    */
    fmpz_add(s + 1, x + 1, x + 5);   /*  a1 = x1 + x5    */
    fmpz_add(s + 2, x + 2, x + 6);   /*  a2 = x2 + x6    */
    fmpz_add(s + 3, x + 3, x + 7);   /*  a3 = x3 + x7    */
    fmpz_sub(s + 4, x + 0, x + 4);   /*  b0 = x0 - x4    */
    fmpz_sub(s + 5, x + 1, x + 5);   /*  b1 = x1 - x5    */
    fmpz_sub(s + 6, x + 2, x + 6);   /*  b2 = x2 - x6    */
    fmpz_sub(s + 7, x + 3, x + 7);   /*  b3 = x3 - x7    */
    _zp_ar2(s + 38, s + 0, s + 4, s);
    fmpz_add(s + 0, x + 0, x + 0);   /*  a0 = x0 + x0    */
    fmpz_add(s + 1, x + 1, x + 1);   /*  a1 = x1 + x1    */
    fmpz_add(s + 2, x + 2, x + 2);   /*  a2 = x2 + x2    */
    fmpz_add(s + 3, x + 3, x + 3);   /*  a3 = x3 + x3    */
    _zp_ar2(s + 8, s + 0, x + 4, s);
    fmpz_sub(s + 16, s + 38, s + 12);   /*  d7 = d0 - c4    */
    fmpz_swap(s + 0, s + 16);   /*  y0 = d7 mod n   */
    fmpz_sub(s + 16, s + 39, s + 13);   /*  d7 = d1 - c5    */
    fmpz_swap(s + 1, s + 16);   /*  y1 = d7 mod n   */
    fmpz_sub(s + 16, s + 40, s + 14);   /*  d7 = d2 - c6    */
    fmpz_swap(s + 2, s + 16);   /*  y2 = d7 mod n   */
    /* y3 = d3 mod n: output 3 stays in s + 41 */
    fmpz_add(s + 16, s + 42, s + 8);   /*  d7 = d4 + c0    */
    fmpz_swap(s + 3, s + 16);   /*  y4 = d7 mod n   */
    fmpz_add(s + 16, s + 43, s + 9);   /*  d7 = d5 + c1    */
    fmpz_swap(s + 8, s + 16);   /*  y5 = d7 mod n   */
    fmpz_add(s + 16, s + 44, s + 10);   /*  d7 = d6 + c2    */
    /* y6 = d7 mod n: output 6 stays in s + 16 */
    /* y7 = c3 mod n: output 7 stays in s + 11 */
    /* reduce the outputs (unreduced, in the scratch slots listed) into
       f, after all reads of the inputs so that f may alias g or h */
    _zp_store(f, s, (const unsigned char[8]) {0, 1, 2, 41, 3, 8, 16, 11}, 8);
}

/* per-p^k modulus bit thresholds up to which the straight-line programs
   are used (0 = none available); measured against the generic path */
static slong
_unity_zp_special_max_bits(ulong p, ulong exp)
{
    ulong pk = n_pow(p, exp);

    switch (pk)
    {
        case 3: return 6000;
        case 4: return 6000;
        case 5: return 6000;
        case 7: return 6000;
        case 8: return 6000;
        case 9: return 6000;
        case 16: return 6000;
        default: return 0;
    }
}

int
_unity_zp_mul_special(unity_zp f, const unity_zp g, const unity_zp h, fmpz * s)
{
    slong d = (f->p - 1) * n_pow(f->p, f->exp - 1);

    if (g->poly->length != d || h->poly->length != d)
        return 0;

    if (fmpz_bits(fmpz_mod_ctx_modulus(f->ctx)) > (ulong) _unity_zp_special_max_bits(f->p, f->exp))
        return 0;

    switch (n_pow(f->p, f->exp))
    {
        case 3: _zp_mul3(f, g, h, s); break;
        case 4: _zp_mul4(f, g, h, s); break;
        case 5: _zp_mul5(f, g, h, s); break;
        case 7: _zp_mul7(f, g, h, s); break;
        case 8: _zp_mul8(f, g, h, s); break;
        case 9: _zp_mul9(f, g, h, s); break;
        case 16: _zp_mul16(f, g, h, s); break;
        default: return 0;
    }

    return 1;
}

int
_unity_zp_sqr_special(unity_zp f, const unity_zp g, fmpz * s)
{
    slong d = (f->p - 1) * n_pow(f->p, f->exp - 1);

    if (g->poly->length != d)
        return 0;

    if (fmpz_bits(fmpz_mod_ctx_modulus(f->ctx)) > (ulong) _unity_zp_special_max_bits(f->p, f->exp))
        return 0;

    switch (n_pow(f->p, f->exp))
    {
        case 3: _zp_sqr3(f, g, s); break;
        case 4: _zp_sqr4(f, g, s); break;
        case 5: _zp_sqr5(f, g, s); break;
        case 7: _zp_sqr7(f, g, s); break;
        case 8: _zp_sqr8(f, g, s); break;
        case 9: _zp_sqr9(f, g, s); break;
        case 16: _zp_sqr16(f, g, s); break;
        default: return 0;
    }

    return 1;
}

void
unity_zp_mul(unity_zp f, const unity_zp g, const unity_zp h)
{
    slong i;

    if (g->poly->length == 0 || h->poly->length == 0)
    {
        fmpz_mod_poly_zero(f->poly, f->ctx);
        return;
    }

    {
        UNITY_ZP_MUL_BEGIN(t, g->poly->length, h->poly->length)
        _unity_zp_mul_reduce(f, g, h, t);
        UNITY_ZP_MUL_END(t)
    }
}

void
unity_zp_mul_inplace(unity_zp f, const unity_zp g, const unity_zp h, fmpz_t * t)
{
    /* the preallocated space holds SQUARING_SPACE values; larger products
       (of unreduced operands or for larger p^k) get their own space */
    if (g->poly->length + h->poly->length - 1 <= SQUARING_SPACE)
    {
        if (!_unity_zp_mul_special(f, g, h, (fmpz *) t))
            _unity_zp_mul_reduce(f, g, h, (fmpz *) t);
    }
    else
        unity_zp_mul(f, g, h);
}

void
unity_zp_mul_scalar_ui(unity_zp f, const unity_zp g, ulong s)
{
    FLINT_ASSERT(fmpz_equal(fmpz_mod_ctx_modulus(f->ctx),
                            fmpz_mod_ctx_modulus(g->ctx)));

    fmpz_mod_poly_scalar_mul_ui(f->poly, g->poly, s, f->ctx);
}
