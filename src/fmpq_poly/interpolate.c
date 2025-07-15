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
#include "fmpq_vec.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"


// computes A / (b*x-a) where x = a / b is a root of A = prod_i(bi*x-ai)
void
_fmpz_poly_div_root_fmpq(fmpz * Q, const fmpz * A, slong len, const fmpq_t x)
{
    if (len < 2)
        return;

    fmpz_t r, t;

    fmpz_init(r);
    fmpz_init(t);

    // r = A[len - 1] / b
    fmpz_divexact(r, A + len - 1, fmpq_denref(x));

    for (slong i = len - 2; i > 0; i--)
    {
        // Q[i] = (A[i] + Q[i+1] * a) / b
        fmpz_set(t, A + i);
        fmpz_addmul(t, r, fmpq_numref(x));
        fmpz_divexact(t, t, fmpq_denref(x));

        fmpz_swap(Q + i, r);
        fmpz_swap(r, t);
    }
    // Q[0] = (A[0] + Q[1] * a) / b
    fmpz_swap(Q, r);

    fmpz_clear(r);
    fmpz_clear(t);
}


void
_fmpq_poly_interpolate_fmpq_vec(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    slong i, j;

    /* Constant */
    if (n == 1)
    {
        fmpz_set(poly, fmpq_numref(ys));
        fmpz_set(den, fmpq_denref(ys));
        return;
    }

    /* Linear */
    if (n == 2)
    {
        fmpz_t Dx, Dy;
        fmpz *xi, *yi;

        fmpz_init(Dx);
        fmpz_init(Dy);
        xi = _fmpz_vec_init(2);
        yi = _fmpz_vec_init(2);

        fmpz_lcm(Dx, fmpq_denref(xs), fmpq_denref(xs + 1));
        fmpz_lcm(Dy, fmpq_denref(ys), fmpq_denref(ys + 1));

        for (i = 0; i < 2; i++) {
            // xi = Dx * xs
            fmpz_divexact(xi + i, Dx, fmpq_denref(xs + i));
            fmpz_mul(xi + i, xi + i, fmpq_numref(xs + i));
            // yi = Dy * ys
            fmpz_divexact(yi + i, Dy, fmpq_denref(ys + i));
            fmpz_mul(yi + i, yi + i, fmpq_numref(ys + i));
        }
        fmpz_sub(den, xi, xi + 1);
        fmpz_mul(den, den, Dy);

        fmpz_sub(poly + 1, yi, yi + 1);
        fmpz_mul(poly, xi, yi + 1);
        fmpz_submul(poly, xi + 1, yi);
        return;
    }

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
            if (i != j)
            {
                fmpz_mul(t, fmpq_numref(xs + i), fmpq_denref(xs + j));
                fmpz_submul(t, fmpq_numref(xs + j), fmpq_denref(xs + i));
                fmpz_mul(w + i, w + i, t);
            }
        }
    }

    _fmpz_vec_zero(poly, n);
    _fmpz_vec_lcm(den, w, n);

    for (i = 0; i < n; i++)
    {
        /* Q = P / (b[i]*x - a[i]) with xs[i]=a[i]/b[i]*/
        _fmpz_poly_div_root_fmpq(Q, P, n + 1, xs + i);

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
fmpq_poly_interpolate_fmpq_vec(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
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
        _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}

void
fmpq_poly_interpolate_fmpz_fmpq_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpq * ys, slong n)
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
        fmpq *xs1 = _fmpq_vec_init(n);
        _fmpq_vec_set_fmpz_vec(xs1, xs, n);

        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs1, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
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
        fmpq *xs1 = _fmpq_vec_init(n);
        fmpq *ys1 = _fmpq_vec_init(n);
        _fmpq_vec_set_fmpz_vec(xs1, xs, n);
        _fmpq_vec_set_fmpz_vec(ys1, ys, n);

        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs1, ys1, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}

