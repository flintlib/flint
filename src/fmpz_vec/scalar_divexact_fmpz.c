/*
    Copyright (C) 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"

void
_fmpz_vec_scalar_divexact_fmpz(fmpz * vec1, const fmpz * vec2,
                               slong len2, const fmpz_t x)
{
    fmpz c = *x;

    flint_printf("\n_FMPZ_VEC_SCALAR_DIVEXACT_FMPZ\n");
    flint_printf("len = %wd\n", len2);
    flint_printf("x = "); fmpz_print(x); flint_printf("\n");
    flint_printf("vec = [ "); _fmpz_vec_print(vec2, len2); flint_printf(" ]\n\n");
    fflush(stdout);

    if (!COEFF_IS_MPZ(c))
    {
        if (c == 1)
            _fmpz_vec_set(vec1, vec2, len2);
        else if (c == -1)
            _fmpz_vec_neg(vec1, vec2, len2);
        else
            _fmpz_vec_scalar_divexact_si(vec1, vec2, len2, c);
    }
    else
    {
        slong i;
        for (i = 0; i < len2; i++)
        {
            flint_printf("i = %wd\n", i); fflush(stdout);
            fmpz_divexact(vec1 + i, vec2 + i, x);
        }
    }
}
