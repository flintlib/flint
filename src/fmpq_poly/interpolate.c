/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_interpolate_fmpz_vec(fmpz * poly, fmpz_t den,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz *P, *Q, *w;
    fmpz_t t;

    slong i, j;

    /* Constant */
    if (n == 1)
    {
        fmpz_set(poly, ys);
        fmpz_one(den);
        return;
    }

    /* Linear */
    if (n == 2)
    {
        fmpz_sub(den, xs, xs + 1);
        fmpz_sub(poly + 1, ys, ys + 1);
        fmpz_mul(poly, xs, ys + 1);
        fmpz_submul(poly, xs + 1, ys);
        return;
    }

    fmpz_init(t);

    P = _fmpz_vec_init(n + 1);
    Q = _fmpz_vec_init(n);
    w = _fmpz_vec_init(n);

    /* P = (x-x[0])*(x-x[1])*...*(x-x[n-1]) */
    _fmpz_poly_product_roots_fmpz_vec(P, xs, n);

    /* Weights */
    for (i = 0; i < n; i++)
    {
        fmpz_one(w + i);
        for (j = 0; j < n; j++)
        {
            if (i != j)
            {
                fmpz_sub(t, xs + i, xs + j);
                fmpz_mul(w + i, w + i, t);
            }
        }
    }

    _fmpz_vec_zero(poly, n);
    _fmpz_vec_lcm(den, w, n);

    for (i = 0; i < n; i++)
    {
        /* Q = P / (x - x[i]) */
        _fmpz_poly_div_root(Q, P, n + 1, xs + i);

        /* result += Q * weight(i) */
        fmpz_divexact(t, den, w + i);
        fmpz_mul(t, t, ys + i);
        _fmpz_vec_scalar_addmul_fmpz(poly, Q, n, t);
    }

    _fmpz_vec_clear(P, n + 1);
    _fmpz_vec_clear(Q, n);
    _fmpz_vec_clear(w, n);
    fmpz_clear(t);
}

void
fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(poly);
    }
    else if (n == 1)
    {
        fmpq_poly_set_fmpz(poly, ys);
    }
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fmpz_vec(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}

void
fmpq_poly_interpolate_fmpz_fmpq_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpq * ys, slong n)
{
    fmpz *yi;
    fmpz_t Y_den, tmp;

    slong i;

    if (n == 0)
    {
        fmpq_poly_zero(poly);
        return;
    }

    if (n == 1)
    {
        fmpq_poly_set_fmpq(poly, ys);
        return;
    }

    yi = _fmpz_vec_init(n);
    fmpz_init(Y_den);
    fmpz_init(tmp);

    // LCM of denominators of the ys -> Y_den
    fmpz_set_ui(Y_den, 1);
    for (i = 0; i < n; i++)
        fmpz_lcm(Y_den, Y_den, fmpq_denref(ys + i));

    // Convert to integer samples
    for (i = 0; i < n; i++)
    {
        // yi = ci * (Y_den / di) where ys = ci/di
        fmpz_divexact(tmp, Y_den, fmpq_denref(ys + i));
        fmpz_mul(yi + i, fmpq_numref(ys + i), tmp);
    }

    // Allocation and interpolation
    fmpq_poly_fit_length(poly, n);
    _fmpq_poly_interpolate_fmpz_vec(poly->coeffs, poly->den, xs, yi, n);

    // Ajust denominators: P / Y_den
    fmpz_mul(poly->den, poly->den, Y_den);
    _fmpq_poly_set_length(poly, n);
    fmpq_poly_canonicalise(poly);

    _fmpz_vec_clear(yi, n);
    fmpz_clear(Y_den);
    fmpz_clear(tmp);
}

void
fmpq_poly_interpolate_fmpq_vec(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    fmpz *xi, *yi;
    fmpz_t X_den, Y_den, tmp;

    slong i;

    if (n == 0)
    {
        fmpq_poly_zero(poly);
        return;
    }

    if (n == 1)
    {
        fmpq_poly_set_fmpq(poly, ys);
        return;
    }

    xi = _fmpz_vec_init(n);
    yi = _fmpz_vec_init(n);
    fmpz_init(X_den);
    fmpz_init(Y_den);
    fmpz_init(tmp);

    // LCM of denominators of the xs -> X_den
    fmpz_set_ui(X_den, 1);
    for (i = 0; i < n; i++)
        fmpz_lcm(X_den, X_den, fmpq_denref(xs + i));

    // LCM of denominators of the ys -> Y_den
    fmpz_set_ui(Y_den, 1);
    for (i = 0; i < n; i++)
        fmpz_lcm(Y_den, Y_den, fmpq_denref(ys + i));

    // Convert to integer samples
    for (i = 0; i < n; i++)
    {
        // xi = ai * (X_den / bi) where xs = ai/bi
        fmpz_divexact(tmp, X_den, fmpq_denref(xs + i));
        fmpz_mul(xi + i, fmpq_numref(xs + i), tmp);

        // yi = ci * (Y_den / di) where ys = ci/di
        fmpz_divexact(tmp, Y_den, fmpq_denref(ys + i));
        fmpz_mul(yi + i, fmpq_numref(ys + i), tmp);
    }

    // Allocation and interpolation
    fmpq_poly_fit_length(poly, n);
    _fmpq_poly_interpolate_fmpz_vec(poly->coeffs, poly->den, xi, yi, n);

    // Adjust variable: P(x/X_den) => multiply x^k coeff by X_den^k
    fmpz_set(tmp, X_den);
    for (i = 1; i < n; i++)
    {
        fmpz_mul(poly->coeffs + i, poly->coeffs + i, tmp);
        fmpz_mul(tmp, tmp, X_den); // tmp = X_den^i
    }
    // Ajust denominators: P / Y_den
    fmpz_mul(poly->den, poly->den, Y_den);
    _fmpq_poly_set_length(poly, n);
    fmpq_poly_canonicalise(poly);

    _fmpz_vec_clear(xi, n);
    _fmpz_vec_clear(yi, n);
    fmpz_clear(X_den);
    fmpz_clear(Y_den);
    fmpz_clear(tmp);
}
