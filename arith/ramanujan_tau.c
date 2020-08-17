/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "arith.h"

void arith_ramanujan_tau_series(fmpz_poly_t res, slong n)
{
    slong j, k, jv, kv;
    fmpz_t tmp;
    fmpz_poly_fit_length(res, n);
    _fmpz_vec_zero(res->coeffs, n);
    _fmpz_poly_set_length(res, n);
    fmpz_init(tmp);
    for (j = jv = 0; jv < n; jv += ++j)
    {
        fmpz_set_ui(tmp, 2*j+1);
        for (k = kv = 0; jv + kv < n; kv += ++k)
        {
            if ((j+k) & 1)
                fmpz_submul_ui(res->coeffs + jv+kv, tmp, 2*k+1);
            else
                fmpz_addmul_ui(res->coeffs + jv+kv, tmp, 2*k+1);
        }
    }
    fmpz_poly_sqrlow(res, res, n-1);
    fmpz_poly_sqrlow(res, res, n-1);
    fmpz_poly_shift_left(res, res, 1);
    fmpz_clear(tmp);
}

void _arith_ramanujan_tau(fmpz_t res, fmpz_factor_t factors)
{
    fmpz_poly_t poly;
    fmpz_t tau_p, p_11, next, this, prev;
    slong k, r;
    ulong max_prime;

    max_prime = UWORD(1);
    for (k = 0; k < factors->num; k++)
    {
        /* TODO: handle overflow properly */
        max_prime = FLINT_MAX(max_prime, fmpz_get_ui(factors->p + k));
    }

    fmpz_poly_init(poly);
    arith_ramanujan_tau_series(poly, max_prime + 1);

    fmpz_one(res);
    fmpz_init(tau_p);
    fmpz_init(p_11);
    fmpz_init(next);
    fmpz_init(this);
    fmpz_init(prev);

    for (k = 0; k < factors->num; k++)
    {
        ulong p = fmpz_get_ui(factors->p + k);

        fmpz_set(tau_p, poly->coeffs + p);
        fmpz_set_ui(p_11, p);
        fmpz_pow_ui(p_11, p_11, 11);
        fmpz_one(prev);
        fmpz_set(this, tau_p);

        for (r = 1; r < factors->exp[k]; r++)
        {
            fmpz_mul(next, tau_p, this);
            fmpz_submul(next, p_11, prev);
            fmpz_set(prev, this);
            fmpz_set(this, next);
        }
        fmpz_mul(res, res, this);
    }

    fmpz_clear(tau_p);
    fmpz_clear(p_11);
    fmpz_clear(next);
    fmpz_clear(this);
    fmpz_clear(prev);
    fmpz_poly_clear(poly);
}

void arith_ramanujan_tau(fmpz_t res, const fmpz_t n)
{
    fmpz_factor_t factors;

    if (fmpz_sgn(n) <= 0)
    {
        fmpz_zero(res);
        return;
    }

    fmpz_factor_init(factors);
    fmpz_factor(factors, n);
    _arith_ramanujan_tau(res, factors);
    fmpz_factor_clear(factors);
}
