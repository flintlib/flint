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

    printf("scalar_submul_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    /* Compare with fmpz_vec_scalar_submul_si */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c;
        len_t len, n;
        fmpz_t n1;
        len = n_randint(state, 100);
        n = (len_t) n_randbits(state, FLINT_BITS - 1);
        if (n_randint(state, 2))
            n = -n;
        fmpz_init(n1);
        fmpz_set_si(n1, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);
        _fmpz_vec_set(c, b, len);

        _fmpz_vec_scalar_submul_fmpz(b, a, len, n1);
        _fmpz_vec_scalar_submul_si(c, a, len, n);

        result = (_fmpz_vec_equal(c, b, len));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(c, len), printf("\n\n");
            _fmpz_vec_print(b, len), printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);
    }

    /* Compute a different way */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c, *d;
        len_t len = n_randint(state, 100);
        fmpz_t n1;
        fmpz_init(n1);
        fmpz_randtest(n1, state, 200);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);
        d = _fmpz_vec_init(len);
        _fmpz_vec_randtest(a, state, len, 200);
        _fmpz_vec_randtest(b, state, len, 200);
        _fmpz_vec_set(c, b, len);

        _fmpz_vec_scalar_submul_fmpz(b, a, len, n1);
        _fmpz_vec_scalar_mul_fmpz(d, a, len, n1);
        _fmpz_vec_sub(c, c, d, len);

        result = (_fmpz_vec_equal(c, b, len));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(c, len), printf("\n\n");
            _fmpz_vec_print(b, len), printf("\n\n");
            abort();
        }

        fmpz_clear(n1);
        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);
        _fmpz_vec_clear(d, len);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
