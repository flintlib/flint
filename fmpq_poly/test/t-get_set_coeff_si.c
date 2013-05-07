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
#include <gmp.h>
#include "flint.h"
#include "fmpq_poly.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, j, result;
    ulong cflags = 0UL;

    long n;
    mpq_t n_mpq;
    flint_rand_t state;
    flint_randinit(state);

    printf("get/set_coeff_si....");
    fflush(stdout);

    mpq_init(n_mpq);

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        fmpq_poly_t a;
        long coeff, len;

        fmpq_poly_init(a);
        len = (long) n_randint(state, 100) + 1;

        for (j = 0; j < 1000; j++)
        {
            n = z_randtest(state);
            coeff = n_randint(state, len);
            fmpq_poly_set_coeff_si(a, coeff, n);
            fmpq_poly_get_coeff_mpq(n_mpq, a, coeff);

            cflags |= fmpq_poly_is_canonical(a) ? 0 : 1;
            result = (mpz_cmp_ui(mpq_denref(n_mpq), 1) == 0 
                   && mpz_cmp_si(mpq_numref(n_mpq), n) == 0
                   && !cflags);
            if (!result)
            {
                printf("FAIL:\n");
                printf("a      = "), fmpq_poly_debug(a), printf("\n");
                printf("len    = %ld\n", len);
                printf("coeff  = %ld\n", coeff);
                printf("cflags = %lu\n", cflags);
                printf("n      = %ld\n", n);
                gmp_printf("n_mpq  = %Qd\n", n_mpq);
                abort();
            }
        }
        fmpq_poly_clear(a);
    }

    flint_randclear(state);
    mpq_clear(n_mpq);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
