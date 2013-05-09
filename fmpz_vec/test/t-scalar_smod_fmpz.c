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

    Copyright (C) 2011 Sebastian Pancratz

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

    printf("scalar_smod_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p;
        fmpz *a, *b;
        len_t len = n_randint(state, 100);

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 1);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_smod_fmpz(b, a, len, p);
        _fmpz_vec_scalar_smod_fmpz(a, a, len, p);

        result = (_fmpz_vec_equal(a, a, len));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, len), printf("\n\n");
            _fmpz_vec_print(b, len), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(p);
    }

    /* Check the result is reduced */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz_t p, lo, hi;
        fmpz *a, *b;
        len_t j, len = n_randint(state, 100);

        fmpz_init(p);
        fmpz_init(lo);
        fmpz_init(hi);
        fmpz_randtest_unsigned(p, state, 100);
        fmpz_add_ui(p, p, 1);
        if (fmpz_cmp_ui(p, 2) > 0)
        {
            fmpz_fdiv_q_2exp(hi, p, 1);
            fmpz_neg(lo, hi);
        }
        else if (fmpz_cmp_ui(p, 2) == 0)
        {
            fmpz_zero(lo);
            fmpz_one(hi);
        }
        else
        {
            fmpz_zero(lo);
            fmpz_zero(hi);
        }

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_scalar_smod_fmpz(b, a, len, p);

        result = 1;
        for (j = 0; j < len; j++)
            result &= (fmpz_cmp(lo, b + j) <= 0 && fmpz_cmp(b + j, hi) <= 0);

        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, len), printf("\n\n");
            _fmpz_vec_print(b, len), printf("\n\n");
            fmpz_print(p), printf("\n\n");
            fmpz_print(lo), printf("\n\n");
            fmpz_print(hi), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        fmpz_clear(p);
        fmpz_clear(lo);
        fmpz_clear(hi);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
