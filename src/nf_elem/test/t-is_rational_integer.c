/*=============================================================================

    This file is part of Antic.

    Antic is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version. See <http://www.gnu.org/licenses/>.

=============================================================================*/
/******************************************************************************

    Copyright (C) 2018 Vincent Delecroix
                  2020 Julian Rüth

******************************************************************************/

#include "flint.h"
#include "fmpq_poly.h"
#include "nf.h"
#include "nf_elem.h"

int main(void)
{
    int i;
    flint_rand_t state;

    flint_printf("set_si_ui...");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        int is_int;
        int is_rat;

        fmpq_poly_t f;
        nf_t nf;
        nf_elem_t a;

        nf_init_randtest(nf, state, 20, 200);

        nf_elem_init(a, nf);

        fmpq_poly_init(f);

        fmpq_poly_randtest(f, state, fmpq_poly_degree(nf->pol) - 1, 200);
        nf_elem_set_fmpq_poly(a, f, nf);

        is_rat = fmpq_poly_length(f) <= 1;
        is_int = is_rat && fmpz_is_one(fmpq_poly_denref(f));

        if (nf_elem_is_rational(a, nf) != is_rat ||
            nf_elem_is_integer(a, nf) != is_int)
            {
                flint_printf("nf_elem_is_rational/nf_elem_is_integer wrong\n");
                flint_printf("nf = "); nf_print(nf); flint_printf("\n");
                flint_printf("f = "); fmpq_poly_print_pretty(f, "x"); flint_printf("\n");
                flint_printf("a = "); nf_elem_print_pretty(a, nf, "x"); flint_printf("\n");
                flint_abort();
            }

        nf_elem_clear(a, nf);
        fmpq_poly_clear(f);
        nf_clear(nf);
    }

    flint_randclear(state);
    flint_cleanup();
    flint_printf("PASS\n");
    return 0;
}
