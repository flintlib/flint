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

    Copyright (C) 2012 Sebastian Pancratz 
    Copyright (C) 2013 Mike Hansen

******************************************************************************/
#include <stdio.h>
#include <stdlib.h>

#include "fq.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"
#include "long_extras.h"

int
main(void)
{
    int i, j, k, result;
    flint_rand_t state;

    flint_printf("ctx_init... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 30; i++) {
        fmpz_t p;
        long d;
        fq_ctx_t ctx;
        
        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 50), 1));
        d = n_randint(state, 20) + 1;

        fq_ctx_init(ctx, p, d, "a");
        fq_ctx_clear(ctx);
    }

    for (i = 0; i < 30; i++) {
        fmpz_t p;
        long d;
        fq_ctx_t ctx_conway, ctx_mod;
        fmpz_mod_poly_t modulus;

        fq_t a, b, lhs, rhs;

        fmpz_init(p);
        fmpz_set_ui(p, n_randprime(state, 2 + n_randint(state, 3), 1));
        d = n_randint(state, 10) + 1;
        fq_ctx_init_conway(ctx_conway, p, d, "a");

        fmpz_mod_poly_init(modulus, p);
        fmpz_mod_poly_zero(modulus);
        for (j = 0; j < ctx_conway->len; j++)
            fmpz_mod_poly_set_coeff_fmpz(modulus, ctx_conway->j[j], &ctx_conway->a[j]);

        fq_ctx_init_modulus(ctx_mod, p, d, modulus, "a");

        fq_init(a);
        fq_init(b);
        fq_init(lhs);
        fq_init(rhs);

        for (k = 0; k < 30; k++)
        {

            fq_randtest(a, state, ctx_conway);
            fq_set(b, a);

            fq_mul(lhs, a, a, ctx_conway);
            fq_mul(rhs, b, b, ctx_mod);

            result = (fq_equal(lhs, rhs));
            if (!result)
            {
                flint_printf("FAIL:\n\n");
                flint_printf("a   = "), fq_print_pretty(a, ctx_conway), flint_printf("\n");
                flint_printf("b   = "), fq_print_pretty(b, ctx_mod), flint_printf("\n");
                flint_printf("lhs = "), fq_print_pretty(lhs, ctx_conway), flint_printf("\n");
                flint_printf("rhs = "), fq_print_pretty(rhs, ctx_mod), flint_printf("\n");
                abort();
            }

        }

        fq_clear(a);
        fq_clear(b);
        fq_clear(lhs);
        fq_clear(rhs);

        fmpz_mod_poly_clear(modulus);
        fq_ctx_clear(ctx_conway);
        fq_ctx_clear(ctx_mod);

    }


    flint_randclear(state);
    _fmpz_cleanup();
    flint_printf("PASS\n");

    return EXIT_SUCCESS;
}
