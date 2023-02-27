/*
    Copyright (C) 2008, 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/* Assumes len > 0. */
void _fmpz_poly_sqr_classical(fmpz *rop, const fmpz *op, slong len)
{
    if (len == 1)  /* Special case */
    {
        fmpz_mul(rop, op, op);
    }
    else   /* Ordinary case */
    {
        slong i;

        _fmpz_vec_scalar_mul_fmpz(rop, op, len, op);

        _fmpz_vec_scalar_mul_fmpz(rop + len, op + 1, len - 1, op + len - 1);

        for (i = 1; i < len - 1; i++)
            _fmpz_vec_scalar_addmul_fmpz(rop + i + 1, op + 1, i - 1, op + i);

        for (i = 1; i < 2 * len - 2; i++)
            fmpz_mul_ui(rop + i, rop + i, 2);

        for (i = 1; i < len - 1; i++)
            fmpz_addmul(rop + 2 * i, op + i, op + i);
    }
}

void fmpz_poly_sqr_classical(fmpz_poly_t rop, const fmpz_poly_t op)
{
    slong len;

    if (op->length == 0)
    {
        fmpz_poly_zero(rop);
        return;
    }

    len = 2 * op->length - 1;

    if (rop == op)
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, len);
        _fmpz_poly_sqr_classical(t->coeffs, op->coeffs, op->length);
        fmpz_poly_swap(rop, t);
        fmpz_poly_clear(t);
    }
    else
    {
        fmpz_poly_fit_length(rop, len);
        _fmpz_poly_sqr_classical(rop->coeffs, op->coeffs, op->length);
    }

    _fmpz_poly_set_length(rop, len);
}

