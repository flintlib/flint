/*
    Copyright (C) 2010 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <gmp.h>

#include "fmpz.h"

int fmpz_print(const fmpz_t x)
{
	if (!COEFF_IS_MPZ(*x)) 
        return flint_printf("%wd", *x);
	else 
        return (int) mpz_out_str(stdout, 10, COEFF_TO_PTR(*x));
}

