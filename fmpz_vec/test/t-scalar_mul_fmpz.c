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
    Copyright (C) 2010 Sebastian Pancratz

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
    printf("scalar_mul_fmpz....");
    fflush(stdout);

    _fmpz_vec_randinit();

    // Check aliasing of a and b
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b;
        fmpz_t n;
        ulong length = n_randint(100);
        fmpz_init(n);
        fmpz_randtest(n, 100);
        if (n_randint(2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        _fmpz_vec_randtest(a, length, n_randint(200));

        _fmpz_vec_scalar_mul_fmpz(b, a, length, n);
        _fmpz_vec_scalar_mul_fmpz(a, a, length, n);

        result = (_fmpz_vec_equal(a, b, length));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, length), printf("\n\n");
            _fmpz_vec_print(b, length), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, length);
        _fmpz_vec_clear(b, length);
        fmpz_clear(n);
    }

    // Check that n (a + b) == na + nb
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b, *lhs, *rhs;
        fmpz_t n;
        ulong length = n_randint(100);
        fmpz_init(n);
        fmpz_randtest(n, 100);
        if (n_randint(2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        lhs = _fmpz_vec_init(length);
        rhs = _fmpz_vec_init(length);
        _fmpz_vec_randtest(a, length, n_randint(200));
        _fmpz_vec_randtest(b, length, n_randint(200));

        _fmpz_vec_scalar_mul_fmpz(lhs, a, length, n);
        _fmpz_vec_scalar_mul_fmpz(rhs, b, length, n);
        _fmpz_vec_add(rhs, lhs, rhs, length);
        _fmpz_vec_add(lhs, a, b, length);
        _fmpz_vec_scalar_mul_fmpz(lhs, lhs, length, n);

        result = (_fmpz_vec_equal(lhs, rhs, length));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, length), printf("\n\n");
            _fmpz_vec_print(b, length), printf("\n\n");
            _fmpz_vec_print(lhs, length), printf("\n\n");
            _fmpz_vec_print(rhs, length), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, length);
        _fmpz_vec_clear(b, length);
        _fmpz_vec_clear(lhs, length);
        _fmpz_vec_clear(rhs, length);
        fmpz_clear(n);
    }

    // Check that n2 * (n1 a) == (n1 * n2) a
    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz *a, *b;
        fmpz_t n1, n2, n;
        ulong length = n_randint(100);
        fmpz_init(n1);
        fmpz_init(n2);
        fmpz_init(n);
        fmpz_randtest(n1, 100);
        fmpz_randtest(n2, 100);
        if (n_randint(2))
            fmpz_neg(n1, n1);
        if (n_randint(2))
            fmpz_neg(n2, n2);

        a = _fmpz_vec_init(length);
        b = _fmpz_vec_init(length);
        _fmpz_vec_randtest(a, length, n_randint(200));

        _fmpz_vec_scalar_mul_fmpz(b, a, length, n1);
        _fmpz_vec_scalar_mul_fmpz(b, b, length, n2);
        fmpz_mul(n, n1, n2);
        _fmpz_vec_scalar_mul_fmpz(a, a, length, n);

        result = (_fmpz_vec_equal(a, b, length));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, length), printf("\n\n");
            _fmpz_vec_print(b, length), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, length);
        _fmpz_vec_clear(b, length);
        fmpz_clear(n1);
        fmpz_clear(n2);
        fmpz_clear(n);
    }

    _fmpz_vec_randclear();

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
