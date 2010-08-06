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

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    gmp_randstate_t state;
    gmp_randinit_default(state);

    printf("get/set_mpz....");
    fflush(stdout);

    for (i = 0; i < 100000; i++)
    {
        fmpz_t a;
        mpz_t b, c;
        mp_bitcnt_t bits;

        mpz_init(b);
        mpz_init(c);

        bits = n_randint(200) + 1;

        mpz_rrandomb(b, state, bits);
        if (n_randint(2))
            mpz_neg(b, b);

        fmpz_init(a);

        fmpz_set_mpz(a, b);
        fmpz_get_mpz(c, a);

        result = (mpz_cmp(b, c) == 0);
        if (!result)
        {
            printf("FAIL:\n");
            gmp_printf("b = %Zd, c = %Zd\n", b, c);
            abort();
        }

        fmpz_clear(a);

        mpz_clear(b);
        mpz_clear(c);
    }

    gmp_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
