/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of Arb.

    Arb is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "bernoulli.h"
#include "flint/ulong_extras.h"
#include "flint/nmod_poly.h"
#include "flint/nmod_vec.h"

int main()
{
    slong iter;
    flint_rand_t state;
    slong n, bound;
    mp_limb_t p, pinv, m1, m2;
    nmod_poly_t A;

    flint_printf("fmpq_vec....");
    fflush(stdout);
    flint_randinit(state);

    bound = 1000 * FLINT_MIN(1.0, arb_test_multiplier());

    p = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 0);
    pinv = n_preinvert_limb(p);

    nmod_poly_init(A, p);
    nmod_poly_set_coeff_ui(A, 1, 1);
    nmod_poly_exp_series(A, A, bound);
    nmod_poly_shift_right(A, A, 1);
    nmod_poly_inv_series(A, A, bound);

    m1 = 1;
    for (n = 0; n < A->length; n++)
    {
        A->coeffs[n] = n_mulmod2_preinv(A->coeffs[n], m1, p, pinv);
        m1 = n_mulmod2_preinv(m1, n + 1, p, pinv);
    }

    for (iter = 0; iter < 100 * arb_test_multiplier(); iter++)
    {
        slong a, b, num;
        fmpq * res;
        slong i;

        b = n_randint(state, bound);
        a = n_randint(state, b + 1);
        num = b - a;

        flint_set_num_threads(1 + n_randint(state, 4));

        res = _fmpq_vec_init(num);

        for (i = 0; i < num; i++)
            fmpq_randtest(res + i, state, 100);

        bernoulli_fmpq_vec_no_cache(res, a, num);

        for (i = 0; i < num; i++)
        {
            m1 = fmpz_fdiv_ui(fmpq_numref(res + i), p);
            m2 = fmpz_fdiv_ui(fmpq_denref(res + i), p);
            m2 = n_invmod(m2, p);
            m1 = n_mulmod2_preinv(m1, m2, p, pinv);
            m2 = nmod_poly_get_coeff_ui(A, a + i);

            if (m1 != m2)
            {
                flint_printf("FAIL:\n");
                flint_printf("a = %wd, b = %wd, num = %wd, n = %wd\n", a, b, num, a + i);
                flint_printf("m1 = %wu mod %wu\n", m1, p);
                flint_printf("m2 = %wu mod %wu\n", m2, p);
                flint_abort();
            }
        }

        _fmpq_vec_clear(res, num);
    }

    nmod_poly_clear(A);

    flint_randclear(state);
    flint_cleanup_master();
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}

