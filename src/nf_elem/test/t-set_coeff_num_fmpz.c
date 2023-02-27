/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2018 Tommy Hofmann
                  2020 Julian Rüth

******************************************************************************/

#include <stdio.h>
#include "nf.h"
#include "nf_elem.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    flint_printf("set_coeff_num_fmpz....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * antic_test_multiplier(); i++)
    {
        nf_t nf;
        nf_elem_t a, b;
        fmpz_t d, d2;
        slong coeff;
        fmpq_t newcoeff, tempcoeff;

        nf_init_randtest(nf, state, 40, 200);

        fmpz_init(d);
        fmpz_init(d2);
        fmpq_init(tempcoeff);
        fmpq_init(newcoeff);

        nf_elem_init(a, nf);
        nf_elem_init(b, nf);
        nf_elem_randtest(a, state, 200, nf);
        nf_elem_set(b, a, nf);

        coeff = (slong) n_randint(state, fmpq_poly_degree(nf->pol));

	fmpz_randtest(d, state, 200);

        nf_elem_get_den(fmpq_denref(tempcoeff), a, nf);
        fmpz_set(fmpq_numref(tempcoeff), d);
        fmpq_canonicalise(tempcoeff);

        _nf_elem_set_coeff_num_fmpz(a, coeff, d, nf);
        nf_elem_get_coeff_fmpq(newcoeff, a, coeff, nf);

        result = fmpq_equal(newcoeff, tempcoeff);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); printf("\n");
            flint_printf("b = "); nf_elem_print_pretty(b, nf, "x"); printf("\n");
            flint_printf("coeff = %u\n", coeff);
            flint_printf("d = "); fmpz_print(d); printf("\n");
            flint_printf("newcoeff = "); fmpq_print(newcoeff); printf("\n");
            flint_printf("pol = "); fmpq_poly_print_pretty(nf->pol, "x"); printf("\n");
            abort();
        }

        nf_elem_clear(a, nf);
        nf_elem_clear(b, nf);
        
        nf_clear(nf);

        fmpz_clear(d);
        fmpz_clear(d2);

        fmpq_clear(tempcoeff);
        fmpq_clear(newcoeff);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
