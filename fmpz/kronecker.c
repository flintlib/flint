/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "fmpz.h"

int fmpz_kronecker(const fmpz_t a, const fmpz_t n)
{
    fmpz A = *a;
    fmpz N = *n;
    mpz_t aa, nn;
    int r;

    if (!COEFF_IS_MPZ(A) && !COEFF_IS_MPZ(N))
        return z_kronecker(A, N);

    if (COEFF_IS_MPZ(A) && COEFF_IS_MPZ(N))
        return mpz_kronecker(COEFF_TO_PTR(A), COEFF_TO_PTR(N));

    flint_mpz_init_set_readonly(aa, a);
    flint_mpz_init_set_readonly(nn, n);

    r = mpz_kronecker(aa, nn);

    flint_mpz_clear_readonly(aa);
    flint_mpz_clear_readonly(nn);

    return r;
}

