/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2010, 2012 Sebastian Pancratz
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "gmpcompat.h"
#include "fmpz.h"

int fmpz_equal(const fmpz_t f, const fmpz_t g)
{
    if (f == g) return 1;  /* aliased inputs */

    if (!COEFF_IS_MPZ(*f)) return (*f == *g);  /* if f is large it can't be equal to g */
    else if (!COEFF_IS_MPZ(*g)) return 0;  /* f is large, so if g isn't... */
    else return (mpz_cmp(COEFF_TO_PTR(*f), COEFF_TO_PTR(*g)) == 0);
}

int fmpz_equal_si(const fmpz_t f, slong g)
{
    fmpz c = *f;

    return !COEFF_IS_MPZ(c) ? (c == g) : !flint_mpz_cmp_si(COEFF_TO_PTR(c), g);
}

int fmpz_equal_ui(const fmpz_t f, ulong g)
{
    fmpz c = *f;

    return !COEFF_IS_MPZ(c) ? ((c >= 0) & (c == g)) :
                              !flint_mpz_cmp_ui(COEFF_TO_PTR(c), g);
}
