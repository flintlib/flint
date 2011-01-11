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
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("content....");
    fflush(stdout);

    flint_randinit(state);

    /* Check that content(a f) = abs(a) content(f) */
    for (i = 0; i < 10000; i++)
    {
        fmpz_t a, b;
        fmpq_poly_t f;
        mpq_t x, y, z;

        fmpz_init(a);
        fmpz_init(b);
        fmpq_poly_init(f);
        fmpq_poly_randtest(f, state, n_randint(100, state), 200);
        fmpz_randtest(a, state, 100);
        fmpz_randtest_not_zero(b, state, 100);
        mpq_init(x);
        mpq_init(y);
        mpq_init(z);

        fmpz_get_mpz(mpq_numref(x), a);
        fmpz_get_mpz(mpq_denref(x), b);
        mpq_canonicalize(x);

        fmpq_poly_content(y, f);
        fmpq_poly_scalar_mul_mpq(f, f, x);
        mpq_abs(x, x);
        mpq_mul(y, x, y);
        fmpq_poly_content(z, f);

        result = (mpq_equal(y, z));
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("y = %Qd\nz = %Qd\n", y, z);
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(b);
        mpq_clear(x);
        mpq_clear(y);
        mpq_clear(z);
        fmpq_poly_clear(f);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
