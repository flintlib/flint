/*
    Copyright (C) 2012, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"
#include "fmpz_extras.h"
#include "bernoulli.h"
#include "bernoulli/impl.h"

/* B_n mod p from the series x / (e^x - 1) mod p (A), for all the numbers
   (or the first count) of the iterator from nmax, with the test flags
   (the generic table kernel, split factors, the sieve by 3 or by 3 and 5
   or none) and threads chosen at random; a few numbers are also compared
   exactly with _bernoulli_fmpq_ui */

static void
check_range(flint_rand_t state, const nmod_poly_t A, ulong p, ulong pinv,
    slong nmax, slong count, int flags, int threads)
{
    fmpz_t numer, denom, num2, den2;
    bernoulli_rev_t iter;
    slong n;
    ulong m1, m2;

    if (flags < 0)
    {
        static const int sv[4] = { 0, BERNOULLI_REV_SIEVE_3,
            BERNOULLI_REV_SIEVE_15, BERNOULLI_REV_NO_SIEVE };
        flags = n_randint(state, 4) | sv[n_randint(state, 4)];
    }
    if (threads < 0)
        threads = 1 + n_randint(state, 4);

    flint_set_num_threads(threads);

    fmpz_init(numer);
    fmpz_init(denom);
    fmpz_init(num2);
    fmpz_init(den2);

    _bernoulli_rev_init_flags(iter, nmax, flags);

    for (n = nmax - (nmax % 2); n >= 0 && count > 0; n -= 2, count--)
    {
        bernoulli_rev_next(numer, denom, iter);

        m1 = fmpz_fdiv_ui(numer, p);
        m2 = fmpz_fdiv_ui(denom, p);
        m2 = n_invmod(m2, p);
        m1 = n_mulmod2_preinv(m1, m2, p, pinv);
        m2 = nmod_poly_get_coeff_ui(A, n);

        if (m1 != m2)
        {
            flint_printf("FAIL (mod p):\n");
            flint_printf("nmax = %wd, n = %wd, flags = %d, threads = %d\n",
                nmax, n, flags, threads);
            flint_printf("m1 = %wu mod %wu\n", m1, p);
            flint_printf("m2 = %wu mod %wu\n", m2, p);
            flint_abort();
        }

        if (n_randint(state, 100) == 0)
        {
            _bernoulli_fmpq_ui(num2, den2, n);
            if (!fmpz_equal(numer, num2) || !fmpz_equal(denom, den2))
            {
                flint_printf("FAIL (exact):\n");
                flint_printf("nmax = %wd, n = %wd, flags = %d, threads = %d\n",
                    nmax, n, flags, threads);
                flint_abort();
            }
        }
    }

    bernoulli_rev_clear(iter);

    fmpz_clear(numer);
    fmpz_clear(denom);
    fmpz_clear(num2);
    fmpz_clear(den2);

    flint_set_num_threads(1);
}

TEST_FUNCTION_START(bernoulli_rev, state)
{
    slong nmax, n, bound, iter;
    ulong p, pinv, m1;
    nmod_poly_t A;

    bound = 100000 * FLINT_MIN(1.0, 0.1 * flint_test_multiplier());
    bound = FLINT_MAX(bound, 10000);

    p = n_nextprime(UWORD(1) << (FLINT_BITS - 1), 0);
    pinv = n_preinvert_limb(p);

    nmod_poly_init(A, p);
    nmod_poly_set_coeff_ui(A, 1, 1);
    nmod_poly_exp_series(A, A, bound + 1);
    nmod_poly_shift_right(A, A, 1);
    nmod_poly_inv_series(A, A, bound + 1);

    m1 = 1;
    for (n = 0; n < A->length; n++)
    {
        A->coeffs[n] = n_mulmod2_preinv(A->coeffs[n], m1, p, pinv);
        m1 = n_mulmod2_preinv(m1, n + 1, p, pinv);
    }

    /* the first numbers from nmax up to the bound */
    for (nmax = 0; nmax <= bound; nmax = 1.5 * nmax + 2)
        check_range(state, A, p, pinv, nmax, (nmax < 8000) ? 4000 : 100, -1, -1);

    /* all the numbers from random nmax, reaching the threaded outputs
       (from about n = 2600) now and then */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        nmax = n_randint(state, (n_randint(state, 4) == 0) ? 9000 : 2000);
        check_range(state, A, p, pinv, nmax, WORD_MAX, -1, -1);
    }

    /* the threaded outputs with each kernel, with and without the sieves */
    for (iter = 0; iter < 8; iter++)
        check_range(state, A, p, pinv, 3000 + n_randint(state, 3000), WORD_MAX,
            (iter % 4) | ((iter < 4) ? BERNOULLI_REV_SIEVE_15 : BERNOULLI_REV_SIEVE_3),
            2 + n_randint(state, 3));

    nmod_poly_clear(A);

    TEST_FUNCTION_END(state);
}
