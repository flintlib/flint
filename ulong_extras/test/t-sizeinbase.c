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
#include <string.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    mp_limb_t n;
    int base, size1, size2;
    len_t rep;
    mpz_t t;
    char * str;

    flint_rand_t state;

    printf("sizeinbase....");
    fflush(stdout);

    flint_randinit(state);
    mpz_init(t);
    str = flint_malloc((FLINT_BITS + 1) * sizeof(char));

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        n = n_randtest(state);
        base = 2 + n_randint(state, 34);

        size1 = n_sizeinbase(n, base);

        mpz_set_ui(t, n);

        mpz_get_str(str, base, t);
        size2 = strlen(str);

        if (size1 != size2)
        {
            printf("FAIL: n = %lu, base = %d\n", n, base);
            printf("n_sizeinbase: %d, strlen: %d\n", size1, size2);
            abort();
        }
    }

    flint_free(str);
    mpz_clear(t);

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
