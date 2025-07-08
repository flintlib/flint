/*
    Copyright (C) 2011 Fredrik Johansson

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
_fmpq_poly_interpolate_fmpz_fmpq_vec(fmpz * poly, fmpz_t den,
                                    const fmpz * xs, const fmpq * ys, slong n)
{
    fmpz *P, *Q, *w;
    fmpz_t t;

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
        fmpz_t Dy;
        fmpz *yi;

        fmpz_init(Dy);
        yi = _fmpz_vec_init(2);

        fmpz_lcm(Dy, fmpq_denref(ys), fmpq_denref(ys + 1));

        for (i=0; i < 2; i++) {
            fmpz_divexact(yi + i, Dy, fmpq_denref(ys + i));
            fmpz_mul(yi + i, yi + i, fmpq_numref(ys + i));
        }
        fmpz_sub(den, xs, xs + 1);
        fmpz_mul(den, den, Dy);

        fmpz_sub(poly + 1, yi, yi + 1);
        fmpz_mul(poly, xs, yi + 1);
        fmpz_submul(poly, xs + 1, yi);
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
        fmpz_set(w + i, fmpq_denref(ys + i));
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
        fmpz_mul(t, t, fmpq_numref(ys + i));
        _fmpz_vec_scalar_addmul_fmpz(poly, Q, n, t);
    }

    _fmpz_vec_clear(P, n + 1);
    _fmpz_vec_clear(Q, n);
    _fmpz_vec_clear(w, n);
    fmpz_clear(t);
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
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_fmpz_fmpq_vec(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
}

void
_fmpz_poly_div_root_fmpq(fmpz *Q, const fmpz *A, slong len, const fmpq_t x)
{
    if (len < 2)
        return;

    fmpz_t r, t;

    fmpz_init(r);
    fmpz_init(t);

    // r = A[len - 1]
    fmpz_divexact(r, A + len - 1, fmpq_denref(x));

    for (slong i = len - 2; i > 0; i--)
    {
        //x = a / b
        // t = r * a
        fmpz_mul(t, r, fmpq_numref(x));
        // t = t + b * A[i]
        fmpz_add(t, t, A + i);
        fmpz_divexact(t, t, fmpq_denref(x));
        // Q[i] = r
        fmpz_swap(Q + i, r);
        // r = t
        fmpz_swap(r, t);
    }

    // Q[0] = r
    fmpz_swap(Q, r);

    fmpz_clear(r);
    fmpz_clear(t);
}

void
_fmpq_poly_product_roots_fmpq_vec(fmpz * poly, const fmpq * xs, slong n)
{
    if (n == 0)
    {
        fmpz_one(poly);
    }
    else if (n < 20)
    {
        slong i, j;

        fmpz_set(poly + n, fmpq_denref(xs));
        fmpz_neg(poly + n - 1, fmpq_numref(xs));

        for (i = 1; i < n; i++)
        {
            fmpz_mul(poly + n - i - 1, poly + n - i, fmpq_numref(xs + i));
            fmpz_neg(poly + n - i - 1, poly + n - i - 1);
            for (j = 0; j < i; j++) {
                fmpz_mul(poly + n - i + j, poly + n - i + j, fmpq_denref(xs + i));
                fmpz_submul(poly + n - i + j, poly + n - i + j + 1, fmpq_numref(xs + i));
            }
            fmpz_mul(poly + n, poly + n, fmpq_denref(xs + i));
        }

    }
    else
    {
        slong m;
        fmpz *t_poly;

        m = (n + 1) / 2;

        t_poly = _fmpz_vec_init(n + 2);

        _fmpq_poly_product_roots_fmpq_vec(t_poly, xs, m);
        _fmpq_poly_product_roots_fmpq_vec(t_poly + m + 1, xs + m, n - m);
        _fmpz_poly_mul(poly, t_poly, m + 1, t_poly + m + 1, n - m + 1);

        _fmpz_vec_clear(t_poly, n + 2);
    }
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

        for (i=0; i < 2; i++) {
            fmpz_divexact(xi + i, Dx, fmpq_denref(xs + i));
            fmpz_mul(xi + i, xi + i, fmpq_numref(xs + i));
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

    fmpz *P, *Q;
    fmpz *w;
    fmpz_t t, t1;

    fmpz_init(t);
    fmpz_init(t1);

    P = _fmpz_vec_init(n + 1);
    Q = _fmpz_vec_init(n);
    w = _fmpz_vec_init(n);

    /* P = (x-x[0])*(x-x[1])*...*(x-x[n-1]) */
    _fmpq_poly_product_roots_fmpq_vec(P, xs, n);

    /* Weights */
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
    fmpz_one(den);
    for (i=0; i < n; i++)
        fmpz_lcm(den, den, w + i);

    for (i = 0; i < n; i++)
    {
        /* Q = P / (x - x[i]) */ //This has to be adapted
        _fmpz_poly_div_root_fmpq(Q, P, n + 1, xs + i);

        /* result += Q * weight(i) */

        fmpz_divexact(t, den, w + i);

        fmpz_pow_ui(t1, fmpq_denref(xs + i), n - 1);
        fmpz_mul(t, t, t1);
        fmpz_mul(t, t, fmpq_numref(ys + i));

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