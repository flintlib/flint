/*
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmp.h"
#include "flint.h"
#include "fmpz-conversions.h"
#ifdef LONGSLONG
# define flint_mpz_cmp_ui mpz_cmp_ui
#else
# include "gmpcompat.h"
#endif

int fmpz_equal_ui(const fmpz_t f, ulong g)
{
    fmpz c = *f;

    return !COEFF_IS_MPZ(c) ? ((c >= 0) & (c == g)) : 
                              !flint_mpz_cmp_ui((mpz_ptr) COEFF_TO_PTR(c), g);
}
