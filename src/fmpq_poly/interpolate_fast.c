/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "profiler.h"
#include "longlong.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_vec.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"


void
_fmpq_poly_interpolate_fast_precomp(fmpz * poly, fmpz_t den, const fmpq * xs,
    const fmpq * ys, fmpz * const * tree, const fmpz * weights, slong len)
{
    fmpz * t, * u,  * pa, * pb;
    fmpz_t t1;
    slong i, pow, left;
    timeit_t t0;

    if (len == 0)
        return;

    t = _fmpz_vec_init(2 * len);
    u = t + len;
    fmpz_init(t1);

    timeit_start(t0);
    _fmpz_vec_lcm(den, weights, len);
    /* Integer coefficients are (den / w[i]) * c[i] * b[i]^(n-1) */
    for (i = 0; i < len; i++) {
        fmpz_divexact(poly + i, den, weights + i);
        fmpz_mul(poly + i, poly + i, fmpq_numref(ys + i));
        fmpz_pow_ui(t1, fmpq_denref(xs + i), len - 1);
        fmpz_mul(poly + i, poly + i, t1);
    }
    timeit_stop(t0);
    flint_printf("Int coeffs: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);
    timeit_start(t0);
    for (i = 0; i < FLINT_CLOG2(len); i++)
    {
        pow = (WORD(1) << i);
        pa = tree[i];
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
    timeit_stop(t0);
    flint_printf("Poly: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);

    _fmpz_vec_clear(t, 2 * len);
    fmpz_clear(t1);
}

void
_fmpq_poly_interpolate_fast(fmpz * poly, fmpz_t den,
    const fmpq * xs, const fmpq * ys, slong len)
{
    fmpz ** tree;
    fmpz * w;
    fmpz_t t;
    slong i, j;
    timeit_t t0;

    timeit_start(t0);
    // Build tree of subproducts of (bi*x-ai) where xs[i]=ai/bi
    tree = _fmpz_poly_tree_alloc(len);
    _fmpz_poly_tree_build_fmpq_vec(tree, xs, len);
    timeit_stop(t0);
    flint_printf("Tree building: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);

    timeit_start(t0);
    // Build efficiently (in practice) the interpolation weights
    fmpz_init(t);
    w = _fmpz_vec_init(len);
    /* Weights: w[i] = d[i] * prod_(j!=i) (a[i]*b[j] - a[j]*b[j]) with ys[i] = c[i]/d[i]*/
    for (i = 0; i < len; i++)
    {
        fmpz_set(w + i, fmpq_denref(ys + i));
        for (j = 0; j < len; j++)
        {
            if (i != j)
            {
                fmpz_mul(t, fmpq_numref(xs + i), fmpq_denref(xs + j));
                fmpz_submul(t, fmpq_numref(xs + j), fmpq_denref(xs + i));
                fmpz_mul(w + i, w + i, t);
            }
        }
    }
    timeit_stop(t0);
    flint_printf("Weights: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);

    timeit_start(t0);
    // Perform fast polynomial operations using the precomputed tree and weights (avoiding dens)
    _fmpq_poly_interpolate_fast_precomp(poly, den, xs, ys, tree, w, len);
    timeit_stop(t0);
    flint_printf("Interp: cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);


    _fmpz_vec_clear(w, len);
    fmpz_clear(t);
    _fmpz_poly_tree_free(tree, len);
}

void
fmpq_poly_interpolate_fast(fmpq_poly_t poly, const fmpq * xs, const fmpq * ys, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(poly);
    }
    else if (n == 1)
    {
        fmpq_poly_set_fmpq(poly, ys);
    }
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fast(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}

