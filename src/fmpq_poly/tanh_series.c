/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpq_poly.h"

void
_fmpq_poly_tanh_series(fmpz * g, fmpz_t gden,
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    fmpz * t;
    fmpz * u;
    fmpz_t tden;
    fmpz_t uden;

    hlen = FLINT_MIN(hlen, n);

    t = _fmpz_vec_init(n);
    u = _fmpz_vec_init(n);
    fmpz_init(tden);
    fmpz_init(uden);

    /* tanh(x) = (exp(2x)-1)/(exp(2x)+1) */
    _fmpq_poly_scalar_mul_ui(t, tden, h, hden, hlen, UWORD(2));
    _fmpq_poly_exp_series(u, uden, t, tden, hlen, n);
    _fmpz_vec_set(t, u, n);
    fmpz_set(tden, uden);
    fmpz_zero(t);               /* t[0] = 0 */
    fmpz_mul_ui(u, uden, UWORD(2));  /* u[0] = 2 */
    _fmpq_poly_div_series(g, gden, t, tden, n, u, uden, n, n);
    _fmpq_poly_canonicalise(g, gden, n);

    _fmpz_vec_clear(t, n);
    _fmpz_vec_clear(u, n);
    fmpz_clear(tden);
    fmpz_clear(uden);
}

void fmpq_poly_tanh_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length && !fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_tanh_series). Constant term != 0.\n");
    }

    if (poly->length == 0 || n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_tanh_series(res->coeffs, res->den,
            poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_tanh_series(t->coeffs, t->den,
            poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
