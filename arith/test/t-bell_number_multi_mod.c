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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <mpfr.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
    flint_rand_t state;
    long i;

    printf("bell_number_multi_mod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        long n;
        fmpz_t b1, b2;

        fmpz_init(b1);
        fmpz_init(b2);

        n = n_randint(state, 500);

        arith_bell_number_bsplit(b1, n);
        arith_bell_number_multi_mod(b2, n);

        if (!fmpz_equal(b1, b2))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            abort();
        }

        fmpz_clear(b1);
        fmpz_clear(b2);
    }

    flint_randclear(state);
    mpfr_free_cache();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
