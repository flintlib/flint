/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mpoly_factor.h"

void check_omega(slong lower, slong upper, const fmpz_mpoly_t p, const fmpz_mpoly_ctx_t ctx)
{
    slong i;
    fmpz_mpoly_t q;
    fmpz_mpoly_factor_t g, h;
    fmpz_t omega;

    fmpz_init(omega);
    fmpz_mpoly_factor_init(g, ctx);
    fmpz_mpoly_factor_init(h, ctx);
    fmpz_mpoly_init(q, ctx);

    if (!fmpz_mpoly_factor_zassenhaus(g, p, ctx))
    {
        flint_printf("check factorization 1 could be computed\n");
        flint_abort();
    }

    if (!fmpz_mpoly_factor(h, p, ctx))
    {
        flint_printf("check factorization 2 could be computed\n");
        flint_abort();
    }

    for (i = 0; i < g->num; i++)
    {
        if (g->poly[i].length < 1 || fmpz_sgn(g->poly[i].coeffs + 0) <= 0)
        {
            flint_printf("factorization is not unit normal\n");
            flint_abort();
        }
    }

    fmpz_zero(omega);
    for (i = 0; i < g->num; i++)
        fmpz_add(omega, omega, g->exp + i);

    if (fmpz_cmp_si(omega, lower) < 0 || fmpz_cmp_si(omega, upper) > 0)
    {
        flint_printf("factorization has wrong number of factors\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_expand(q, g, ctx);

    if (!fmpz_mpoly_equal(q, p, ctx))
    {
        flint_printf("factorization does not match original polynomial\n");
        flint_abort();        
    }

    fmpz_mpoly_factor_sort(g, ctx);
    fmpz_mpoly_factor_sort(h, ctx);
    if (fmpz_mpoly_factor_cmp(g, h, ctx) != 0)
    {
        flint_printf("factorizations do not match\n");
        flint_abort();        
    }

    for (i = 0; i < g->num; i++)
    {
        fmpz_mpoly_factor(h, g->poly + i, ctx);
        if (h->num != 1 || !fmpz_is_one(h->exp + 0))
        {
            flint_printf("FAIL:\nfactor is reducible\n");
            flint_abort();
        }
    }

    fmpz_mpoly_clear(q, ctx);
    fmpz_mpoly_factor_clear(g, ctx);
    fmpz_mpoly_factor_clear(h, ctx);
    fmpz_clear(omega);
    return;
}


int
main(void)
{
    slong i, j, tmul = 20;

    FLINT_TEST_INIT(state);

    flint_printf("factor_zassenhaus....");
    fflush(stdout);

    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a;
        fmpz_mpoly_factor_t f;
        const char * vars[] = {"x", "y", "z", "w"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_factor_init(f, ctx);

        fmpz_mpoly_set_str_pretty(a,
            "((x^2+y^2+z^2+2+w)*(x+1)*(y+2)*(z+3+w^3)+x*y*z*w)"
            "*((x^2+y^2+z^2+3+w)*(x^2+1)*(y^2+2)*(z^2+3+w^3)+x*y*z+ w)"
            "*((x^2+y^2+z^2+x+y+z+5+w)*(7+x+y^2)*(8+x+z^2+w^3)+x*y*z+ w)"
            "*((x^2+y^2+z^2+x+y+z+13+w)*(10+y+z^2+w^3)+x*y*z*w)"
        , vars, ctx);

        check_omega(4, 4, a, ctx);

        fmpz_mpoly_set_str_pretty(a,
            "(((x+y+z+2)*(x+y^2+1+z^2)+z+w^4)*((x+2*y+z^2)*(x^2+y+2+z^3)+z+w^6)+x*y*z*w)"
            "*(((x+y+z+7)*(x+y^2+11+z^2)+z+w^4)*((x+2*y+z^2+13)*(x^2+y+2+z^3+17)+z+w^6)+x*y*z*w)"
            "*(x+y+z+w+19)"
        , vars, ctx);

        check_omega(3, 3, a, ctx);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_factor_clear(f, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    for (i = 0; i < tmul * flint_test_multiplier(); i++)
    {
        slong lower;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, t;
        flint_bitcnt_t coeff_bits;
        slong nfacs, len;
        ulong expbound, powbound, pow;

        fmpz_mpoly_ctx_init_rand(ctx, state, 5);

        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(t, ctx);

        nfacs = 1 + (4 + n_randint(state, 5))/ctx->minfo->nvars;
        expbound = 2 + 20/nfacs/ctx->minfo->nvars;
        powbound = 1 + n_randint(state, 3);

        lower = 0;
        fmpz_mpoly_one(a, ctx);
        for (j = 0; j < nfacs; j++)
        {
            do {
                len = 1 + n_randint(state, 10);
                coeff_bits = 10 + n_randint(state, 300)/nfacs;
                fmpz_mpoly_randtest_bound(t, state, len, coeff_bits, expbound, ctx);
            } while (t->length == 0);
            pow = 1 + n_randint(state, powbound);
            if (!fmpz_mpoly_is_fmpz(t, ctx))
                lower += pow;
            fmpz_mpoly_pow_ui(t, t, pow, ctx);
            fmpz_mpoly_mul(a, a, t, ctx);
        }

        check_omega(lower, WORD_MAX, a, ctx);

        fmpz_mpoly_clear(t, ctx);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_ctx_clear(ctx);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
