/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/
#include "math.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

static int
_fmpq_vec_has_unique_entries(const fmpq * x, slong n)
{
    fmpq * t;
    slong i;
    int ok = 1;

    t = _fmpq_vec_init(n);
    for (i = 0; i < n; i++)
        fmpq_set(t + i, x + i);
    _fmpq_vec_sort(t, n);
    for (i = 1; i < n && ok; i++)
        if (fmpq_equal(t + i - 1, t + i))
            ok = 0;

    _fmpq_vec_clear(t, n);

    return ok;
}


int
_fmpq_poly_interpolate_fmpq_vec(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    /* Constant */
    if (n == 1)
    {
        fmpz_set(poly, fmpq_numref(ys));
        fmpz_set(den, fmpq_denref(ys));
        return 1;
    }

    if (!_fmpq_vec_has_unique_entries(xs, n))
        return 0;

    /* Linear */
    if (n == 2)
    {
        slong i;
        fmpz_t Dx, Dy;
        fmpz *xi, *yi;

        fmpz_init(Dx);
        fmpz_init(Dy);
        xi = _fmpz_vec_init(2);
        yi = _fmpz_vec_init(2);

        fmpz_lcm(Dx, fmpq_denref(xs), fmpq_denref(xs + 1));
        fmpz_lcm(Dy, fmpq_denref(ys), fmpq_denref(ys + 1));

        for (i = 0; i < 2; i++) {
            // xs = xi / Dx
            fmpz_divexact(xi + i, Dx, fmpq_denref(xs + i));
            fmpz_mul(xi + i, xi + i, fmpq_numref(xs + i));
            // ys = yi / Dy
            fmpz_divexact(yi + i, Dy, fmpq_denref(ys + i));
            fmpz_mul(yi + i, yi + i, fmpq_numref(ys + i));
        }
        // den = Dy * (x0 - x1)
        fmpz_sub(den, xi, xi + 1);
        fmpz_mul(den, den, Dy);
        // deg1: Dx * (y0 - y1)
        fmpz_sub(poly + 1, yi, yi + 1);
        fmpz_mul(poly + 1, poly + 1, Dx);
        // deg0: y1 * x0 - x1 * y0
        fmpz_mul(poly, yi + 1, xi);
        fmpz_submul(poly, xi + 1, yi);

        fmpz_clear(Dx);
        fmpz_clear(Dy);
        _fmpz_vec_clear(xi, 2);
        _fmpz_vec_clear(yi, 2);
        return 1;
    }

    /*  Estimate on f = (1/D) * sum_{k=1..n-1} fk * X^k,
        whith xi = a/b and yi = c/d:

        ci / di = sum_k fk * ai^k * bi^(n-1-k) / D * bi^(n-1)
        => ht(yi) ~= ht(f) - (n - 1) * ht(xi)
    */
    flint_bitcnt_t xbits, ybits, fbits;
    xbits = _fmpq_vec_max_height_bits(xs, n);
    ybits = _fmpq_vec_max_height_bits(ys, n);
    fbits = FLINT_MAX( ybits - ( n - 1) *  FLINT_MAX(xbits - 1,0), 0);

    if ((fbits < FLINT_BITS -  2 && n > 10) || n > 2 * (sqrt(fbits) + 1))
        _fmpq_poly_interpolate_multi_mod(poly, den, xs, ys, n);
    else if (n < 80)
        _fmpq_poly_interpolate_barycentric(poly, den, xs, ys, n);
    else
        _fmpq_poly_interpolate_fast(poly, den, xs, ys, n);

    return 1;
}

int
fmpq_poly_interpolate_fmpq_vec(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    int ok = 1;
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq_poly_fit_length(poly, n);
        ok = _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);
    }
    return ok;
}

int
fmpq_poly_interpolate_fmpz_fmpq_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpq * ys, slong n)
{
    int ok = 1;
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq *xs1 = _fmpq_vec_init(n);
        _fmpq_vec_set_fmpz_vec(xs1, xs, n);

        fmpq_poly_fit_length(poly, n);
        ok = _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs1, ys, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);

        _fmpq_vec_clear(xs1, n);
    }
    return ok;
}

int
fmpq_poly_interpolate_fmpz_vec(fmpq_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    int ok = 1;
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq *xs1 = _fmpq_vec_init(n);
        fmpq *ys1 = _fmpq_vec_init(n);
        _fmpq_vec_set_fmpz_vec(xs1, xs, n);
        _fmpq_vec_set_fmpz_vec(ys1, ys, n);

        fmpq_poly_fit_length(poly, n);
        ok = _fmpq_poly_interpolate_fmpq_vec(poly->coeffs, poly->den, xs1, ys1, n);
        _fmpq_poly_set_length(poly, n);
        fmpq_poly_canonicalise(poly);

        _fmpq_vec_clear(xs1, n);
        _fmpq_vec_clear(ys1, n);
    }
    return ok;
}
