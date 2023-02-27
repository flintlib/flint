/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
fmpz_jacobi(const fmpz_t a, const fmpz_t p)
{
    fmpz c = *p;
    fmpz d = *a;
    mpz_t t, u;
    int r;

    if (!COEFF_IS_MPZ(c) && !COEFF_IS_MPZ(d))
        return n_jacobi(d, c);

    if (COEFF_IS_MPZ(c) && COEFF_IS_MPZ(d))
        return mpz_jacobi(COEFF_TO_PTR(d), COEFF_TO_PTR(c));

    if (d == 0)
        return 0; /* a is zero and p is large */

    flint_mpz_init_set_readonly(t, a);
    flint_mpz_init_set_readonly(u, p);

    r = mpz_jacobi(t, u);

    flint_mpz_clear_readonly(t);
    flint_mpz_clear_readonly(u);

    return r;
}
