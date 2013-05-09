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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("max_limbs....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a;
        len_t len, bits;
        mp_size_t limbs, limbs2;

        len = n_randint(state, 100);

        a = _fmpz_vec_init(len);
        bits = n_randint(state, 200);
        limbs = (bits + FLINT_BITS - 1) / FLINT_BITS;
        _fmpz_vec_randtest(a, state, len, bits);

        limbs2 = _fmpz_vec_max_limbs(a, len);

        result = (limbs >= limbs2);
        if (!result)
        {
            printf("FAIL:\n");
            printf("bits   = %ld\n", bits);
            printf("limbs  = %ld\n", limbs);
            printf("a      = {"), _fmpz_vec_print(a, len), printf("}\n");
            printf("limbs2 = %ld\n", limbs2);
            abort();
        }

        _fmpz_vec_clear(a, len);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
