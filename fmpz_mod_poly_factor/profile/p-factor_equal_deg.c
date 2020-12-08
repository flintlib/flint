/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include "fmpz_mod_poly.h"
#include "profiler.h"

int main(void)
{
    slong i, j, d, r;
    flint_rand_t state;
    fmpz_t p, t;
    fmpz_mod_ctx_t ctx;
    fmpz_mod_poly_t a, b, c;
    fmpz_mod_poly_factor_t f;
    timeit_t timer;

    flint_randinit(state);
    fmpz_init_set_ui(p, 1);
    fmpz_init(t);
    fmpz_mod_ctx_init(ctx, p);
    fmpz_mod_poly_init(a, ctx);
    fmpz_mod_poly_init(b, ctx);
    fmpz_mod_poly_init(c, ctx);
    fmpz_mod_poly_factor_init(f, ctx);

    while (fmpz_bits(p) < 100)
    {
        if (fmpz_cmp_ui(p, 10) > 0)
        {
            fmpz_sqrt(t, p);
            fmpz_sqrt(t, t);
            fmpz_sqrt(t, t);
            fmpz_add_ui(t, t, 1);
            fmpz_mul(p, p, t);
        }
        fmpz_nextprime(p, p, 0);
        fmpz_mod_ctx_set_modulus(ctx, p);

        flint_printf("++++++++++ p = ");
        fmpz_print(p);
        flint_printf(" ++++++++++\n");

        for (d = 2; d <= 5; d++)
        {
            flint_printf("d = %wd: ", d);
            fmpz_mod_poly_randtest_monic_irreducible(a, state, d + 1, ctx);
            for (r = 1; r <= 10; r++)
            {
                slong reps = a->length + 1;
                reps = 1000000/reps/reps/fmpz_bits(p);

                fmpz_mod_poly_randtest_monic_irreducible(b, state, d + 1, ctx);
                fmpz_mod_poly_mul(a, a, b, ctx);

                timeit_start(timer);
                for (i = 0; i <= reps; i++)
                    fmpz_mod_poly_factor(f, a, ctx);
                timeit_stop(timer);

                flint_printf(" %06wd", timer->wall*100/reps);
                fflush(stdout);

                fmpz_mod_poly_one(c, ctx);
                for (i = 0; i < f->num; i++)
                {
                    if (fmpz_mod_poly_degree(f->poly + i, ctx) != d)
                    {
                        flint_printf("!!! oops !!!");
                        flint_abort();
                    }

                    for (j = 0; j < f->exp[i]; j++)
                        fmpz_mod_poly_mul(c, c, f->poly + i, ctx);
                }

                if (!fmpz_mod_poly_equal(c, a, ctx))
                {
                    flint_printf("!!! oops !!!");
                    flint_abort();
                }
            }

            flint_printf("\n");
        }
    }

    fmpz_mod_poly_clear(a, ctx);
    fmpz_mod_poly_clear(b, ctx);
    fmpz_mod_poly_clear(c, ctx);
    fmpz_mod_poly_factor_clear(f, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(p);
    fmpz_clear(t);
    flint_randclear(state);

    return 0;
}

