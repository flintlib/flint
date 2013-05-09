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
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("get_cfrac....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq_t x, r;
        fmpz *c1, *c2;
        len_t n1, n2, bound;

        fmpq_init(x);
        fmpq_init(r);

        fmpq_randtest(x, state, 1 + n_randint(state, 1000));
        bound = fmpq_cfrac_bound(x);

        c1 = _fmpz_vec_init(bound);
        c2 = _fmpz_vec_init(bound);

        n1 = fmpq_get_cfrac(c1, r, x, bound);

        if (!fmpq_is_zero(r))
        {
            printf("FAIL: expected zero remainder\n");
            abort();
        }

        /* Test chaining */
        n2 = 0;
        while (1)
        {
            n2 += fmpq_get_cfrac(c2 + n2, x, x, n_randint(state, 10));
            if (fmpq_is_zero(x))
                break;
            fmpq_inv(x, x);
        }

        if (n1 != n2)
        {
            printf("FAIL: n1 = %ld, n2 = %ld\n", n1, n2);
            abort();
        }

        if (!_fmpz_vec_equal(c1, c2, n1))
        {
            printf("FAIL: vectors not equal\n");
            _fmpz_vec_print(c1, n1); printf("\n");
            _fmpz_vec_print(c2, n2); printf("\n");
            abort();
        }

        _fmpz_vec_clear(c1, bound);
        _fmpz_vec_clear(c2, bound);
        fmpq_clear(x);
        fmpq_clear(r);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
