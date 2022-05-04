/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "flint.h"
#include "flint-impl.h"
#include "fmpz-conversions.h"

int fmpz_fprint(FILE * file, const fmpz_t x)
{
	if (!COEFF_IS_MPZ(*x))
        return fprintf(file, WORD_FMT "d", *x);
	else 
        return (int) mpz_out_str(file, 10, COEFF_TO_PTR(*x));
}

