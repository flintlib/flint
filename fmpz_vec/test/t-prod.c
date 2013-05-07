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

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2011 Fredrik Johansson

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

    printf("prod....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t x, y, z;

        long len1 = n_randint(state, 100);
        long len2 = n_randint(state, 100);

        a = _fmpz_vec_init(len1 + len2);
        b = a + len1;

        _fmpz_vec_randtest(a, state, len1 + len2, 200);

        fmpz_init(x);
        fmpz_init(y);
        fmpz_init(z);

        _fmpz_vec_prod(x, a, len1);
        _fmpz_vec_prod(y, b, len2);
        fmpz_mul(x, x, y);
        _fmpz_vec_prod(z, a, len1 + len2);

        result = (fmpz_equal(x, z));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, len1), printf("\n\n");
            _fmpz_vec_print(b, len2), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len1 + len2);

        fmpz_clear(x);
        fmpz_clear(y);
        fmpz_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
