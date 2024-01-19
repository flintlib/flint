/*
    Copyright (C) 2008, 2009, 2016 William Hart
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2021 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

void
fmpz_poly_scalar_abs(fmpz_poly_t res, const fmpz_poly_t poly)
{
    fmpz_poly_fit_length(res, poly->length);

    _fmpz_vec_scalar_abs(res->coeffs, poly->coeffs, poly->length);

    _fmpz_poly_set_length(res, poly->length);
}

void
fmpz_poly_scalar_addmul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             const fmpz_t x)
{
    if (!fmpz_is_zero(x) && !fmpz_poly_is_zero(poly2))
    {
        fmpz_poly_fit_length(poly1, poly2->length);

        if (poly2->length > poly1->length)
            _fmpz_vec_zero(poly1->coeffs + poly1->length,
                    poly2->length - poly1->length);

        _fmpz_vec_scalar_addmul_fmpz(poly1->coeffs,
                poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, FLINT_MAX(poly1->length, poly2->length));
        _fmpz_poly_normalise(poly1);
    }
}

void
fmpz_poly_scalar_addmul_si(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             slong x)
{
    if (x != 0 && !fmpz_poly_is_zero(poly2))
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        if (poly2->length > poly1->length)
            _fmpz_vec_zero(poly1->coeffs + poly1->length,
			   poly2->length - poly1->length);
        _fmpz_vec_scalar_addmul_si(poly1->coeffs,
                                    poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, FLINT_MAX(poly1->length, poly2->length));
        _fmpz_poly_normalise(poly1);
    }
}

void
fmpz_poly_scalar_addmul_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             ulong x)
{
    if (x != 0 && !fmpz_poly_is_zero(poly2))
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        if (poly2->length > poly1->length)
            _fmpz_vec_zero(poly1->coeffs + poly1->length,
			   poly2->length - poly1->length);
        _fmpz_vec_scalar_addmul_ui(poly1->coeffs,
                                    poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, FLINT_MAX(poly1->length, poly2->length));
        _fmpz_poly_normalise(poly1);
    }
}

void
fmpz_poly_scalar_divexact_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                               const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_divexact_fmpz). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_divexact_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_divexact_si(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             slong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_divexact_si). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_divexact_si(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_divexact_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             ulong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_divexact_ui). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_divexact_ui(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_fdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           ulong exp)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_fdiv_q_2exp(poly1->coeffs, poly2->coeffs,
        poly2->length, exp);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_fdiv_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_fdiv_fmpz). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_fdiv_q_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_fdiv_si(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                         slong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_fdiv_si). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_fdiv_q_si(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_fdiv_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                         ulong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_fdiv_ui). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_fdiv_q_ui(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_mod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
    }
    else
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_scalar_mod_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, poly2->length);
        _fmpz_poly_normalise(poly1);
    }
}

void
fmpz_poly_scalar_smod_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2, const fmpz_t x)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
    }
    else
    {
        fmpz_poly_fit_length(poly1, poly2->length);
        _fmpz_vec_scalar_smod_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, poly2->length);
        _fmpz_poly_normalise(poly1);

    }
}

void
fmpz_poly_scalar_mul_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                        ulong exp)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_mul_2exp(poly1->coeffs, poly2->coeffs, poly2->length, exp);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_mul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                          const fmpz_t x)
{
    /* Either scalar or input poly is zero */
    if ((*x == 0) || (poly2->length == 0))
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_mul_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_mul_si(fmpz_poly_t poly1, const fmpz_poly_t poly2, slong x)
{
    slong i;

    /* Either scalar or input poly is zero */
    if ((x == WORD(0)) || (poly2->length == 0))
    {
        fmpz_poly_zero(poly1);
        return;
    }

    /* Special case, multiply by 1 */
    if (x == WORD(1))
    {
        fmpz_poly_set(poly1, poly2);
        return;
    }

    /* Special case, multiply by -1 */
    if (x == WORD(-1))
    {
        fmpz_poly_neg(poly1, poly2);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);

    for (i = 0; i < poly2->length; i++)
        fmpz_mul_si(poly1->coeffs + i, poly2->coeffs + i, x);

    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_mul_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2, ulong x)
{
    slong i;

    /* Either scalar or input poly is zero */
    if ((x == 0) || (poly2->length == 0))
    {
        fmpz_poly_zero(poly1);
        return;
    }

    /* Special case, multiply by 1 */
    if (x == 1)
    {
        fmpz_poly_set(poly1, poly2);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);

    for (i = 0; i < poly2->length; i++)
        fmpz_mul_ui(poly1->coeffs + i, poly2->coeffs + i, x);

    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_submul_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                             const fmpz_t x)
{
    if (!fmpz_is_zero(x) && !fmpz_poly_is_zero(poly2))
    {
        fmpz_poly_fit_length(poly1, poly2->length);

        if (poly2->length > poly1->length)
            _fmpz_vec_zero(poly1->coeffs + poly1->length,
                    poly2->length - poly1->length);

        _fmpz_vec_scalar_submul_fmpz(poly1->coeffs,
                poly2->coeffs, poly2->length, x);
        _fmpz_poly_set_length(poly1, FLINT_MAX(poly1->length, poly2->length));
        _fmpz_poly_normalise(poly1);
    }
}

void
fmpz_poly_scalar_tdiv_2exp(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           ulong exp)
{
    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_tdiv_q_2exp(poly1->coeffs, poly2->coeffs,
        poly2->length, exp);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_tdiv_fmpz(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                           const fmpz_t x)
{
    if (fmpz_is_zero(x))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_tdiv_fmpz). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_tdiv_q_fmpz(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_tdiv_si(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                         slong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_tdiv_si). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_tdiv_q_si(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}

void
fmpz_poly_scalar_tdiv_ui(fmpz_poly_t poly1, const fmpz_poly_t poly2,
                         ulong x)
{
    if (x == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (fmpz_poly_scalar_tdiv_ui). Division by zero.\n");
    }

    if (poly2->length == 0)
    {
        fmpz_poly_zero(poly1);
        return;
    }

    fmpz_poly_fit_length(poly1, poly2->length);
    _fmpz_vec_scalar_tdiv_q_ui(poly1->coeffs, poly2->coeffs, poly2->length, x);
    _fmpz_poly_set_length(poly1, poly2->length);
}
