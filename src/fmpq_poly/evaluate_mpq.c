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
#include "fmpq.h"

void 
fmpq_poly_evaluate_mpq(mpq_t res, const fmpq_poly_t poly, const mpq_t a)
{
    fmpq_t r, b;

    fmpq_init(r);
    fmpq_init(b);
    fmpq_set_mpq(b, a);
    fmpq_poly_evaluate_fmpq(r, poly, b);
    fmpq_get_mpq(res, r);
    fmpq_clear(r);
    fmpq_clear(b);
}
