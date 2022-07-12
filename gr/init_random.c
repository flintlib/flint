#include "gr.h"

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)
{
    int which = n_randint(state, 100);

    if (which < 30)
        gr_ctx_init_fmpz(ctx);
    else if (which < 45)
        gr_ctx_init_nmod8(ctx, (n_randtest(state) % 256) || 1);
    else if (which < 55)
    {
        fmpz_t t;
        fmpz_init(t);
        fmpz_randtest_not_zero(t, state, 100);
        fmpz_abs(t, t);
        gr_ctx_init_fmpz_mod(ctx, t);
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
