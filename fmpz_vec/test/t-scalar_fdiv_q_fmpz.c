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

    printf("scalar_fdiv_q_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b, *c;
        fmpz_t n;
        mpz_t d, e, f, m;
        len_t i;
        len_t len = n_randint(state, 100);

        fmpz_init(n);
        fmpz_randtest_not_zero(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _fmpz_vec_init(len);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_set(b, a, len);

        _fmpz_vec_scalar_fdiv_q_fmpz(c, a, len, n);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        for (i = 0; i < len; i++)
        {
            fmpz_get_mpz(m, n);
            fmpz_get_mpz(d, b + i);

            mpz_fdiv_q(e, d, m);

            fmpz_get_mpz(f, c + i);

            result = (mpz_cmp(f, e) == 0);
            if (!result)
            {
                printf("FAIL:\n");
                gmp_printf("d = %Zd, m = %Zd, e = %Zd, f = %Zd\n", d, m, e, f);
                abort();
            }
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _fmpz_vec_clear(c, len);

        fmpz_clear(n);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }    

    /* Test aliasing of a and c */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        fmpz_t n;
        mpz_t d, e, f, m;
        len_t i;
        len_t len = n_randint(state, 100);

        fmpz_init(n);
        fmpz_randtest_not_zero(n, state, 100);
        if (n_randint(state, 2))
            fmpz_neg(n, n);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_set(b, a, len);

        _fmpz_vec_scalar_fdiv_q_fmpz(a, a, len, n);

        mpz_init(d);
        mpz_init(e);
        mpz_init(f);
        mpz_init(m);

        for (i = 0; i < len; i++)
        {
            fmpz_get_mpz(m, n);
            fmpz_get_mpz(d, b + i);

            mpz_fdiv_q(e, d, m);

            fmpz_get_mpz(f, a + i);

            result = (mpz_cmp(f, e) == 0);
            if (!result)
            {
                printf("FAIL:\n");
                gmp_printf("d = %Zd, m = %Zd, e = %Zd, f = %Zd\n", d, m, e, f);
                abort();
            }
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);

        fmpz_clear(n);
        mpz_clear(d);
        mpz_clear(e);
        mpz_clear(f);
        mpz_clear(m);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
