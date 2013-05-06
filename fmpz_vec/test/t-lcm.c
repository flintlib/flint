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

    printf("lcm....");
    fflush(stdout);

    flint_randinit(state);

    /* Check that lcm(a f) = abs(a) lcm(f) */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t a, c, d;
        fmpz *f;
        long len = n_randint(state, 100);

        fmpz_init(a);
        fmpz_init(c);
        fmpz_init(d);
        f = _fmpz_vec_init(len);
        _fmpz_vec_randtest(f, state, len, 200);
        fmpz_randtest(a, state, 100);

        _fmpz_vec_lcm(c, f, len);

        if (len == 0)
        {
            result = fmpz_is_one(c);
        }
        else
        {
            _fmpz_vec_scalar_mul_fmpz(f, f, len, a);
            fmpz_abs(a, a);
            fmpz_mul(c, a, c);
            _fmpz_vec_lcm(d, f, len);
            result = (fmpz_equal(c, d));
        }

        if (!result)
        {
            fmpz_print(c), printf("\n\n");
            fmpz_print(d), printf("\n\n");
            abort();
        }

        fmpz_clear(a);
        fmpz_clear(c);
        fmpz_clear(d);
        _fmpz_vec_clear(f, len);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
