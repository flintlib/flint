/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"


#define FLINT_REVERSE_NEWTON_CUTOFF 4

void
_fmpq_poly_revert_series_newton(fmpz * Qinv, fmpz_t den,
        const fmpz * Q, const fmpz_t Qden, slong Qlen, slong n)
{
    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen <= 2)
    {
        fmpz_zero(Qinv);

        if (Qlen == 2)
        {
            fmpz_set(Qinv + 1, Qden);
            fmpz_set(den, Q + 1);
            _fmpq_poly_canonicalise(Qinv, den, 2);
        }

        _fmpz_vec_zero(Qinv + 2, n - 2);
    }
    else
    {
        fmpz *T, *U, *V;
        fmpz_t Tden, Uden, Vden;

        T = _fmpz_vec_init(n);
        U = _fmpz_vec_init(n);
        V = _fmpz_vec_init(n);
        fmpz_init(Tden);
        fmpz_init(Uden);
        fmpz_init(Vden);

        FLINT_NEWTON_INIT(FLINT_REVERSE_NEWTON_CUTOFF, n)

        FLINT_NEWTON_BASECASE(k)
        _fmpq_poly_revert_series_lagrange(Qinv, den, Q, Qden, Qlen, k);
        _fmpz_vec_zero(Qinv + k, n - k);
        FLINT_NEWTON_END_BASECASE

        FLINT_NEWTON_LOOP(k0, k)
        _fmpq_poly_compose_series(T, Tden, Q, Qden, FLINT_MIN(Qlen, k), Qinv, den, k0, k);
        _fmpq_poly_derivative(U, Uden, T, Tden, k); fmpz_zero(U + k - 1);
        fmpz_zero(T + 1);
        _fmpq_poly_div_series(V, Vden, T, Tden, k, U, Uden, k, k);
        _fmpq_poly_canonicalise(V, Vden, k);
        _fmpq_poly_derivative(T, Tden, Qinv, den, k);
        _fmpq_poly_mullow(U, Uden, V, Vden, k, T, Tden, k, k);
        _fmpq_poly_sub(Qinv, den, Qinv, den, k, U, Uden, k);

        FLINT_NEWTON_END_LOOP
        FLINT_NEWTON_END

        _fmpq_poly_canonicalise(Qinv, den, n);

        _fmpz_vec_clear(T, n);
        _fmpz_vec_clear(U, n);
        _fmpz_vec_clear(V, n);
        fmpz_clear(Tden);
        fmpz_clear(Uden);
        fmpz_clear(Vden);
    }
}

void
fmpq_poly_revert_series_newton(fmpq_poly_t res,
            const fmpq_poly_t poly, slong n)
{
    if (poly->length < 2 || !fmpz_is_zero(poly->coeffs)
                         || fmpz_is_zero(poly->coeffs + 1))
    {
        flint_printf("Exception (fmpq_poly_revert_series_newton). Input must have \n"
               "zero constant term and nonzero coefficient of x^1.\n");
        flint_abort();
    }

    if (n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (res != poly)
    {
        fmpq_poly_fit_length(res, n);
        _fmpq_poly_revert_series_newton(res->coeffs,
                res->den, poly->coeffs, poly->den, poly->length, n);
    }
    else
    {
        fmpq_poly_t t;
        fmpq_poly_init2(t, n);
        _fmpq_poly_revert_series_newton(t->coeffs,
                t->den, poly->coeffs, poly->den, poly->length, n);
        fmpq_poly_swap(res, t);
        fmpq_poly_clear(t);
    }

    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}
