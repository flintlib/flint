/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "longlong.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_interpolation_weights(fmpz * w, fmpz_t wden, const fmpq * xs, slong len)
{
    slong i, j;
    fmpz_t t;
    fmpz_init(t);

    /* Weights: b[i]^(n-1) / w[i] = prod_(j!=i) (a[i]*b[j] - a[j]*b[j]) with ys[i] = c[i]/d[i]*/
    for (i = 0; i < len; i++)
    {
        fmpz_one(w + i);
        for (j = 0; j < len; j++)
        {
            if (i == j) continue;
            fmpz_mul(t, fmpq_numref(xs + i), fmpq_denref(xs + j));
            fmpz_submul(t, fmpq_numref(xs + j), fmpq_denref(xs + i));
            fmpz_mul(w + i, w + i, t);
        }
    }

    /* Same denominator reduction */
    _fmpz_vec_lcm(wden, w, len);
    for (i = 0; i < len; i++) {
        fmpz_pow_ui(t, fmpq_denref(xs + i), len - 1);
        fmpz_mul(t, t, wden);
        fmpz_divexact(w + i, t, w + i);
    }

    fmpz_clear(t);
}


void
_fmpq_poly_interpolate_fast_precomp(fmpz * poly, fmpz_t den,
    const fmpq * ys, fmpz * const * tree, const fmpz * weights, slong len)
{
    fmpz * t, * u,  * pa, * pb;
    fmpz_t yden;
    slong i, pow, left;
    ulong j;

    t = _fmpz_vec_init(2 * len);
    u = t + len;

    fmpz_init(yden);
    fmpz_one(yden);
    for (i = 0; i < len; i++)
        fmpz_lcm(yden, yden, fmpq_denref(ys + i));
    fmpz_mul(den, den, yden);

    /* Integer coefficients are l[i] = (yden / d[i]) * c[i] * w[i] */
    for (i = 0; i < len; i++) {
        fmpz_divexact(poly + i, yden, fmpq_denref(ys + i));
        fmpz_mul(poly + i, poly + i, fmpq_numref(ys + i));
        fmpz_mul(poly + i, poly + i, weights + i );
    }

    /* Compute numerator of sum_j l[j] / (b[j]*x - a[j]) */
    for (j = 0; j < FLINT_CLOG2(len); j++)
    {
        pow = (WORD(1) << j);
        pa = tree[j];
        pb = poly;
        left = len;

        while (left >= 2 * pow)
        {
            _fmpz_poly_mul(t, pa, pow + 1, pb + pow, pow);
            _fmpz_poly_mul(u, pa + pow + 1, pow + 1, pb, pow);
            _fmpz_vec_add(pb, t, u, 2 * pow);

            left -= 2 * pow;
            pa = pa + 2 * pow + 2;
            pb = pb + 2 * pow;
        }

        if (left > pow)
        {
            _fmpz_poly_mul(t, pa, pow + 1, pb + pow, left - pow);
            _fmpz_poly_mul(u, pb, pow, pa + pow + 1, left - pow + 1);
            _fmpz_vec_add(pb, t, u, left);
        }
    }

    _fmpz_vec_clear(t, 2 * len);
    fmpz_clear(yden);
}

/* Perform Lagrange interpolation with basis (bi*x-ai) */
void
_fmpq_poly_interpolate_fast(fmpz * poly, fmpz_t den,
    const fmpq * xs, const fmpq * ys, slong len)
{
    fmpz ** tree;
    fmpz * w;

    // Build tree of subproducts of (bi*x-ai) where xs[i]=ai/bi
    tree = _fmpz_poly_tree_alloc(len);
    _fmpz_poly_tree_build_fmpq_vec(tree, xs, len);

    // Build naively the interpolation weights (cheap)
    w = _fmpz_vec_init(len);
    _fmpq_poly_interpolation_weights(w, den, xs, len);

    // Perform fast polynomial operations using the precomputed tree and weights (avoiding dens)
    _fmpq_poly_interpolate_fast_precomp(poly, den, ys, tree, w, len);

    _fmpz_vec_clear(w, len);
    _fmpz_poly_tree_free(tree, len);
}

void
fmpq_poly_interpolate_fast(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)
{
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fast(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}
