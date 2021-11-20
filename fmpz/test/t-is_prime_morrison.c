/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("is_prime_morrison....");
    fflush(stdout);

    for (i = 0; i < 30 * flint_test_multiplier(); i++)
    {
        fmpz_t p, F, R;
        mp_ptr pp1;
        slong num_pp1;
        double logd;
        ulong limit;

        fmpz_init(p);
        fmpz_init(F);
        fmpz_init(R);

        do {
           fmpz_randbits(p, state, n_randint(state, 330) + 2);
           fmpz_abs(p, p);
        } while (!fmpz_is_probabprime(p) || fmpz_cmp_ui(p, 2) == 0);

        logd = log(fmpz_get_d(p));
        limit = (ulong) (logd*logd*logd/10.0) + 2;

        pp1 = _nmod_vec_init((ulong) logd + 2);
        _fmpz_np1_trial_factors(p, pp1, &num_pp1, limit);

        result = fmpz_is_prime_morrison(F, R, p, pp1, num_pp1);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_print(p); printf("\n");
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(pp1);

        fmpz_clear(p);
        fmpz_clear(F);
        fmpz_clear(R);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
