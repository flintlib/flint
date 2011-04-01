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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

long fmpz_remove(fmpz_t rop, const fmpz_t op, const fmpz_t f)
{
    long ans;
    mpz_t x, y, z;

    if ((fmpz_sgn(f) <= 0) || fmpz_is_one(f))
    {
        printf("Exception:  factor f <= 1 in fmpz_remove\n");
        abort();
    }

    mpz_init(x);
    mpz_init(y);
    mpz_init(z);

    fmpz_get_mpz(y, op);
    fmpz_get_mpz(z, f);
    ans = mpz_remove(x, y, z);
    fmpz_set_mpz(rop, x);

    mpz_clear(x);
    mpz_clear(y);
    mpz_clear(z);

    return ans;
}

