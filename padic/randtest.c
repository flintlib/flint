#include "padic.h"

void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
{
    long min = - (ctx->N / 10), max = ctx->N - 1;

    rop[1] = n_randint(state, max - min + 1) + min;

    {
        fmpz *pow;
        int alloc;

        _padic_ctx_pow_ui(&pow, &alloc, ctx->N - rop[1], ctx);
        fmpz_randm(rop, state, pow);
        if (alloc)
            fmpz_clear(pow);
    }

    padic_normalise(rop, ctx);
}

void padic_randtest_not_zero(padic_t rop, flint_rand_t state, 
                             const padic_ctx_t ctx)
{
    long i;

    padic_randtest(rop, state, ctx);

    for (i = 0; !padic_is_zero(rop, ctx) && i < 10; i++)
        padic_randtest(rop, state, ctx);

    if (padic_is_zero(rop, ctx))
        padic_one(rop, ctx);
}

