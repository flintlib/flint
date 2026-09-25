/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arb.h"
#include "fixed.h"

/* the cached log p_j and 2 arg(pi_j) are the exact floors at every
   prefix length, whether they come from the static tables or from the
   Machin sets by binary splitting (with followup series beyond) */
TEST_FUNCTION_START(fixed_machin_caches, state)
{
    slong iter;

    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        flint_set_num_threads(1 + n_randint(state, 4));
        int gaussian = n_randint(state, 2);
        slong num = 1 + n_randint(state, 64);
        slong nv = 1 + n_randint(state, (iter % 4 == 0) ? 200 : 80);
        slong nr = 1 + n_randint(state, nv);    /* prefix read */
        slong prec = FLINT_BITS * (nv + 4), j;
        n_primes_t iter_p;
        arb_t v, y;
        fmpq_t q;
        fmpz_t f;
        nn_ptr ref;

        arb_init(v); arb_init(y); fmpq_init(q); fmpz_init(f);
        ref = flint_malloc((nv + 1) * sizeof(ulong));

        if (gaussian)
        {
            _fixed_atan_gauss_clear();
            _fixed_atan_gauss_ensure(num, nv);
        }
        else
        {
            _fixed_log_primes_clear();
            _fixed_log_primes_ensure(num, nv);
        }

        n_primes_init(iter_p);
        for (j = 0; j < num; j++)
        {
            nn_srcptr e;

            if (gaussian)
            {
                e = _fixed_atan_gauss_entry(j, nr);
                fmpq_set_si(q, _fixed_gaussian_primes[2 * j + 1],
                    _fixed_gaussian_primes[2 * j]);
                arb_set_fmpq(y, q, prec);
                arb_atan(v, y, prec);
                arb_mul_2exp_si(v, v, 1);
            }
            else
            {
                e = _fixed_log_primes_entry(j, nr);
                arb_log_ui(v, n_primes_next(iter_p), prec);
            }

            arb_mul_2exp_si(v, v, FLINT_BITS * nr);
            arb_floor(v, v, prec + 64);
            if (!arb_get_unique_fmpz(f, v))
                continue;
            fmpz_get_ui_array(ref, nr + 1, f);

            if (mpn_cmp(ref, e, nr + 1) != 0)
            {
                flint_printf("FAIL: %s, num = %wd, nv = %wd, nr = %wd, "
                    "j = %wd\n", gaussian ? "gaussian" : "log", num, nv,
                    nr, j);
                flint_abort();
            }
        }
        n_primes_clear(iter_p);

        arb_clear(v); arb_clear(y); fmpq_clear(q); fmpz_clear(f);
        flint_free(ref);
    }

    flint_set_num_threads(1);
    TEST_FUNCTION_END(state);
}
