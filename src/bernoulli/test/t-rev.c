/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz_extras.h"
#include "bernoulli.h"

TEST_FUNCTION_START(bernoulli_rev, state)
{
    slong nmax, n, bound, count;
    mp_limb_t p, pinv, m1, m2;
    nmod_poly_t A;

    bound = 100000 * FLINT_MIN(1.0, 0.1 * flint_test_multiplier());

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

    for (nmax = 0; nmax < bound; nmax = 1.5 * nmax + 2)
    {
        fmpz_t numer, denom;
        bernoulli_rev_t iter;

        fmpz_init(numer);
        fmpz_init(denom);

        nmax += (nmax % 2);

        bernoulli_rev_init(iter, nmax);

        if (nmax < 8000)
            count = 4000;
        else
            count = 100;

        /* flint_printf("nmax = %wd, count = %wd\n", nmax, count); */

        for (n = nmax; n >= 0 && count > 0; n -= 2, count--)
        {
            bernoulli_rev_next(numer, denom, iter);

            m1 = fmpz_fdiv_ui(numer, p);
            m2 = fmpz_fdiv_ui(denom, p);
            m2 = n_invmod(m2, p);
            m1 = n_mulmod2_preinv(m1, m2, p, pinv);
            m2 = nmod_poly_get_coeff_ui(A, n);

            if (m1 != m2)
            {
                flint_printf("FAIL:\n");
                flint_printf("nmax = %wd, n = %wd\n", nmax, n);
                flint_printf("m1 = %wu mod %wu\n", m1, p);
                flint_printf("m2 = %wu mod %wu\n", m2, p);
                flint_abort();
            }
        }

        bernoulli_rev_clear(iter);

        fmpz_clear(numer);
        fmpz_clear(denom);
    }

    nmod_poly_clear(A);

    TEST_FUNCTION_END(state);
}
