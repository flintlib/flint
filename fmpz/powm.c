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

    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

void
fmpz_powm(fmpz_t f, const fmpz_t g, const fmpz_t e, const fmpz_t m)
{
    if (fmpz_sgn(m) <= 0)
    {
        printf("Exception (fmpz_powm).  Modulus is less than 1.\n");
        abort();
    }

    if (!COEFF_IS_MPZ(*e))
    {
        fmpz_powm_ui(f, g, *e, m);
        return;
    }

    if (fmpz_is_one(m))
    {
        fmpz_zero(f);
        return;
    }

    /* TODO:  Implement this properly! */
    {
        mpz_t f2, g2, e2, m2;

        mpz_init(f2);
        mpz_init(g2);
        mpz_init(e2);
        mpz_init(m2);

        fmpz_get_mpz(f2, f);
        fmpz_get_mpz(g2, g);
        fmpz_get_mpz(e2, e);
        fmpz_get_mpz(m2, m);

        mpz_powm(f2, g2, e2, m2);

        fmpz_set_mpz(f, f2);

        mpz_clear(f2);
        mpz_clear(g2);
        mpz_clear(e2);
        mpz_clear(m2);
    }
}

