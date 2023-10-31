/*
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "nmod_mat.h"
#include "fmpz.h"
#include "fmpz_mat.h"

TEST_FUNCTION_START(nmod_mat_det_howell, state)
{
    slong m, mod, rep;

    for (rep = 0; rep < 1000 * flint_test_multiplier(); rep++)
    {
        nmod_mat_t A;
        fmpz_mat_t B;
        mp_limb_t Adet;
        fmpz_t Bdet;
        ulong t;

	m = n_randint(state, 30);
        mod = n_randtest(state);
	mod += mod == 0;

        nmod_mat_init(A, m, m, mod);
        fmpz_mat_init(B, m, m);

        switch (rep % 3)
        {
            case 0:
                nmod_mat_randrank(A, state, m);
                nmod_mat_randops(A, n_randint(state, 2*m + 1), state);
                break;
            case 1:
                t = n_randint(state, m);
                t = FLINT_MIN(t, m);
                nmod_mat_randrank(A, state, t);
                nmod_mat_randops(A, n_randint(state, 2*m + 1), state);
                break;
            default:
                nmod_mat_randtest(A, state);
        }

        fmpz_mat_set_nmod_mat_unsigned(B, A);

        Adet = nmod_mat_det_howell(A);

        fmpz_init(Bdet);
        fmpz_mat_det_bareiss(Bdet, B);
        fmpz_mod_ui(Bdet, Bdet, mod);

        if (Adet != fmpz_get_ui(Bdet))
        {
            flint_printf("FAIL\n");
	    flint_printf("Adet = %wu, Bdet = %wu\n", Adet, fmpz_get_ui(Bdet));
            fflush(stdout);
            flint_abort();
        }

        nmod_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_clear(Bdet);
    }

    TEST_FUNCTION_END(state);
}
