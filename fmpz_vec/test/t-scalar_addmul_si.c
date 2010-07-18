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
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int result;
    printf("scalar_addmul_si....");
    fflush(stdout);

    _fmpz_vec_randinit();

    // Compare with alternative method of computation
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b, *c, *d;
        ulong length = n_randint(100);
        long x;

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        c = _fmpz_vec_init(length);
        d = _fmpz_vec_init(length);

        _fmpz_vec_randtest(a, length, n_randint(200));
        _fmpz_vec_randtest(b, length, n_randint(200));
        _fmpz_vec_copy(c, b, length);

        x = n_randbits(n_randint(FLINT_BITS - 1));
        if (n_randint(2))
            x = -x;

        _fmpz_vec_scalar_addmul_si(b, a, length, x);
        _fmpz_vec_scalar_mul_si(d, a, length, x);
        _fmpz_vec_add(c, c, d, length);

        result = (_fmpz_vec_equal(b, c, length));
        if (!result)
        {
            printf("FAIL:\n");
            printf("x = %ld\n", x);
            _fmpz_vec_print(b, length), printf("\n\n");
            _fmpz_vec_print(c, length), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, length);
        _fmpz_vec_clear(b, length);
        _fmpz_vec_clear(c, length);
        _fmpz_vec_clear(d, length);
    }

    _fmpz_vec_randclear();

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
