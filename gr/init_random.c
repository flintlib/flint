#include "gr.h"

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state)
{
    int which = n_randint(state, 100);

    if (which < 30)
        gr_ctx_init_fmpz(ctx);
    else if (which < 50)
        gr_ctx_init_nmod8(ctx, 1 + n_randint(state, 255));
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
