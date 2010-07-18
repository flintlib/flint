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
    printf("scalar_addmul_fmpz....");
    fflush(stdout);

    _fmpz_vec_randinit();

    // Compare with fmpz_vec_scalar_addmul_si
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b, *c;
        fmpz_t n1;
        fmpz_init(n1);
        ulong length = n_randint(100);
        long n = (long) n_randbits(FLINT_BITS - 1);
        if (n_randint(2))
            n = -n;
        fmpz_set_si(n1, n);

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        c = _fmpz_vec_init(length);
        _fmpz_vec_randtest(a, length, n_randint(200));
        _fmpz_vec_randtest(b, length, n_randint(200));
        _fmpz_vec_copy(c, b, length);

        _fmpz_vec_scalar_addmul_fmpz(b, a, length, n1);
        _fmpz_vec_scalar_addmul_si(c, a, length, n);

        result = (_fmpz_vec_equal(c, b, length));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(c, length), printf("\n\n");
            _fmpz_vec_print(b, length), printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
        _fmpz_vec_clear(a, length);
        _fmpz_vec_clear(b, length);
        _fmpz_vec_clear(c, length);
    }

    // Compute a different way
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b, *c, *d;
        fmpz_t n1;
        fmpz_init(n1);
        ulong length = n_randint(100);
        fmpz_randtest(n1, 200);

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        c = _fmpz_vec_init(length);
        d = _fmpz_vec_init(length);
        _fmpz_vec_randtest(a, length, n_randint(200));
        _fmpz_vec_randtest(b, length, n_randint(200));
        _fmpz_vec_copy(c, b, length);

        _fmpz_vec_scalar_addmul_fmpz(b, a, length, n1);
        _fmpz_vec_scalar_mul_fmpz(d, a, length, n1);
        _fmpz_vec_add(c, c, d, length);

        result = (_fmpz_vec_equal(c, b, length));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(c, length), printf("\n\n");
            _fmpz_vec_print(b, length), printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
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
