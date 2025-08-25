/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_interpolate_barycentric(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    slong i, j;
    fmpz *P, *Q, *w;
    fmpz_t t, t1;

    fmpz_init(t);
    fmpz_init(t1);

    P = _fmpz_vec_init(n + 1);
    Q = _fmpz_vec_init(n);
    w = _fmpz_vec_init(n);

    /* P = (b[0]*x-a[0])*(b[1]*x-a[1])*...*(b[n-1]*x-a[n-1]) with xs[i]=a[i]/b[i]*/
    _fmpz_poly_product_roots_fmpq_vec(P, xs, n);

    /* Weights: w[i] = d[i] * prod_(j!=i) (a[i]*b[j] - a[j]*b[j]) with ys[i] = c[i]/d[i]*/
    for (i = 0; i < n; i++)
    {
        fmpz_set(w + i, fmpq_denref(ys + i));
        for (j = 0; j < n; j++)
        {
            if (i == j) continue;
            fmpz_mul(t, fmpq_numref(xs + i), fmpq_denref(xs + j));
            fmpz_submul(t, fmpq_numref(xs + j), fmpq_denref(xs + i));
            fmpz_mul(w + i, w + i, t);
        }
    }

    _fmpz_vec_zero(poly, n);
    _fmpz_vec_lcm(den, w, n);

    for (i = 0; i < n; i++)
    {
        /* Q = P / (b[i]*x - a[i]) with xs[i]=a[i]/b[i]*/
        _fmpz_poly_divexact_root_fmpq(Q, P, n + 1, xs + i);

        /* poly += (den / w[i]) * c[i] * b[i]^(n-1) * Q(x) */
        fmpz_divexact(t, den, w + i);
        fmpz_mul(t, t, fmpq_numref(ys + i));
        fmpz_pow_ui(t1, fmpq_denref(xs + i), n - 1);
        fmpz_mul(t, t, t1);

        _fmpz_vec_scalar_addmul_fmpz(poly, Q, n, t);
    }

    _fmpz_vec_clear(P, n + 1);
    _fmpz_vec_clear(Q, n);
    _fmpz_vec_clear(w, n);
    fmpz_clear(t);
    fmpz_clear(t1);
}

void
fmpq_poly_interpolate_barycentric(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_barycentric(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}
