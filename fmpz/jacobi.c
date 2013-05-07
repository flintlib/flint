/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2012 William Hart

******************************************************************************/

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

    if (d == 0)
       return 0;

    if (c == 2)
        return 1;

    if (!COEFF_IS_MPZ(c) && !COEFF_IS_MPZ(d))
        return n_jacobi(d, c);

    if (COEFF_IS_MPZ(c) && COEFF_IS_MPZ(d))
        return mpz_jacobi(COEFF_TO_PTR(d), COEFF_TO_PTR(c));

    flint_mpz_init_set_readonly(t, a);
    flint_mpz_init_set_readonly(u, p);

    r = mpz_jacobi(t, u);

    flint_mpz_clear_readonly(t);
    flint_mpz_clear_readonly(u);

    return r;
}
