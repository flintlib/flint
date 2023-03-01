/*
    Copyright (C) 2011 Sebastian Pancratz
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpq_poly.h"

void fmpq_poly_get_coeff_fmpz(fmpz_t x, const fmpq_poly_t poly, slong n)
{
    if (n >= poly->length)  /* Coefficient is beyond the end of poly */
        fmpz_zero(x);
    else
        fmpz_set(x, poly->coeffs + n);
}

