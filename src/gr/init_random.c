/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"
#include "long_extras.h"
#include "gr.h"

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)
{
    int which = n_randint(state, 100);

    if (which < 25)
        gr_ctx_init_fmpz(ctx);
    else if (which < 40)
        gr_ctx_init_nmod8(ctx, n_randtest(state) % 255 + 1);
    else if (which < 42)
        gr_ctx_init_nmod32(ctx, n_randtest(state) % UWORD(4294967295) + 1);
    else if (which < 45)
        gr_ctx_init_nmod(ctx, n_randtest_not_zero(state));
    else if (which < 50)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_randtest_not_zero(t, state, 100);
        fmpz_abs(t, t);
        gr_ctx_init_fmpz_mod(ctx, t);
        fmpz_clear(t);
    }
    else if (which < 55)
    {
        fmpz_t t;
        fmpz_init(t);

        switch (n_randint(state, 3))
        {
            case 0:
                fmpz_set_ui(t, n_randtest_prime(state, 0));
                gr_ctx_init_fq_nmod(ctx, t, 1 + n_randint(state, 4), NULL);
                break;
            case 1:
                fmpz_set_ui(t, n_randprime(state, 4, 0));
                gr_ctx_init_fq_zech(ctx, t, 1 + n_randint(state, 3), NULL);
                break;
            default:
                fmpz_randprime(t, state, 2 + n_randint(state, 100), 0);
                gr_ctx_init_fq(ctx, t, 1 + n_randint(state, 4), NULL);
        }

        fmpz_clear(t);
    }
    else if (which < 60)
        gr_ctx_init_fmpq(ctx);
    else if (which < 65)
        gr_ctx_init_real_arb(ctx, 2 + n_randint(state, 200));
    else if (which < 70)
        gr_ctx_init_complex_acb(ctx, 2 + n_randint(state, 200));
    else if (which == 75)
        gr_ctx_init_real_ca(ctx);
    else if (which == 76)
        gr_ctx_init_complex_ca(ctx);
    else if (which == 77)
        gr_ctx_init_real_algebraic_ca(ctx);
    else if (which == 78)
        gr_ctx_init_complex_algebraic_ca(ctx);
    else if (which == 79)
    {
        fmpz_poly_t g;
        fmpq_poly_t f;

        fmpz_poly_init(g);
        fmpq_poly_init(f);

        do
        {
            fmpz_poly_randtest_irreducible(g, state, 2 + n_randint(state, 5), 1 + n_randint(state, 10));
        } while (g->length < 2);

        fmpq_poly_set_fmpz_poly(f, g);
        fmpq_poly_scalar_div_ui(f, f, 1 + n_randtest(state) % 256);

        gr_ctx_init_nf(ctx, f);

        fmpz_poly_clear(g);
        fmpq_poly_clear(f);
    }

/*
slow -- but should be ok with degree limits

    else if (which == 98)
        gr_ctx_init_real_qqbar(ctx);
    else if (which == 99)
        gr_ctx_init_complex_qqbar(ctx);
    }
*/
    else
    {
        gr_ctx_init_fmpz(ctx);
    }
}
