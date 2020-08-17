/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int main(void)
{
    mp_ptr b1, b2;
    slong n;

    const slong maxn = 3000;

    FLINT_TEST_INIT(state);

    flint_printf("bell_number_nmod_vec....");
    fflush(stdout);    

    b1 = _nmod_vec_init(maxn);
    b2 = _nmod_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        nmod_t mod;
        mp_limb_t p;

        do {
            p = n_randtest_prime(state, 0);
        } while (p < n);

        nmod_init(&mod, p);

        arith_bell_number_nmod_vec_recursive(b1, n, mod);
        arith_bell_number_nmod_vec_series(b2, n, mod);

        if (!_nmod_vec_equal(b1, b2, n))
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n", n);
            abort();
        }
    }

    _nmod_vec_clear(b1);
    _nmod_vec_clear(b2);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
