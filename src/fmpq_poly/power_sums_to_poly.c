/*
    Copyright (C) 2016 Vincent Delecroix

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_power_sums_to_poly(fmpz * res, const fmpz * poly, const fmpz_t den,
                              slong len)
{
    slong i, k;
    fmpz_t f;
    slong d;
    ulong a;

    fmpz_init(f);
    fmpz_divexact(f, poly, den);
    d = fmpz_get_ui(f);

    fmpz_one(f);
    for (k = 1; k <= d; k++)
    {
        if (k < len)
        {
			fmpz_mul(res + d - k, poly + k, f);
            _fmpz_vec_dot_general(res + d - k, res + d - k, 0, res + d - k + 1, poly + 1, 0, k - 1);

        }
		else
        {
            _fmpz_vec_dot_general(res + d - k, NULL, 0, res + d - k + 1, poly + 1, 0, len - 1);
        }

        a = n_gcd(FLINT_ABS(fmpz_fdiv_ui(res + d - k, k)), k);
        fmpz_divexact_ui(res + d - k, res + d - k, a);
        if (a != (ulong) k)
        {
            a = k / a;
            for (i = d - k + 1; i < d; i++)
                fmpz_mul_ui(res + i, res + i, a);
            fmpz_mul_ui(f, f, a);
        }
        fmpz_neg(res + d - k, res + d - k);
        fmpz_mul(f, f, den);

        for (i = d - k + 1; i < d; i++)
            fmpz_mul(res + i, res + i, den);
    }
    fmpz_swap(res + d, f);
    fmpz_clear(f);
}

void
fmpq_poly_power_sums_to_fmpz_poly(fmpz_poly_t res, const fmpq_poly_t Q)
{
    if (Q->length == 0)
        fmpz_poly_one(res);
    else
    {
        fmpz_t fd;
        slong d;
        fmpz_init(fd);
        fmpz_divexact(fd, Q->coeffs, Q->den);
        d = fmpz_get_ui(fd);
        fmpz_clear(fd);
        fmpz_poly_fit_length(res, d + 1);
        _fmpq_poly_power_sums_to_poly(res->coeffs, Q->coeffs, Q->den,
                                      Q->length);
        _fmpz_poly_set_length(res, d + 1);
        _fmpz_poly_normalise(res);
        _fmpz_poly_primitive_part(res->coeffs, res->coeffs, d + 1);
    }
}

void
fmpq_poly_power_sums_to_poly(fmpq_poly_t res, const fmpq_poly_t Q)
{
    if (Q->length == 0)
        fmpq_poly_one(res);
    else
    {
        fmpz_t fd;
        slong d;
        fmpz_init(fd);
        fmpz_divexact(fd, Q->coeffs, Q->den);
        d = fmpz_get_ui(fd);
        fmpz_clear(fd);

        if (res == Q)
        {
            fmpq_poly_t t;
            fmpq_poly_init(t);
            fmpq_poly_fit_length(res, d + 1);
            _fmpq_poly_power_sums_to_poly(t->coeffs, Q->coeffs, Q->den,
                                          Q->length);
            fmpq_poly_swap(res, t);
            fmpq_poly_clear(t);
        }
        else
        {
            fmpq_poly_fit_length(res, d + 1);
            _fmpq_poly_power_sums_to_poly(res->coeffs, Q->coeffs, Q->den,
                                          Q->length);
        }
        _fmpq_poly_set_length(res, d + 1);
        _fmpq_poly_normalise(res);
        _fmpq_poly_make_monic(res->coeffs, res->den, res->coeffs, res->den,
                              res->length);
    }
}
