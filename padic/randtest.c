#include "padic.h"

void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
{
    rop[1] = n_randint(state, ctx->N);

    /* TODO:  Faster powering */
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, ctx->N - rop[1]);
        fmpz_randm(rop, state, pow);
        fmpz_clear(pow);
    }

    padic_normalise(rop, ctx);
}

void padic_randtest_not_zero(padic_t rop, flint_rand_t state, 
                             const padic_ctx_t ctx)
{
    /* TODO:  Convince yourself this terminates */
    do 
        padic_randtest(rop, state, ctx);
    while (padic_is_zero(rop, ctx));
}

