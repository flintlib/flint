#include "padic.h"

void padic_randtest(padic_t rop, flint_rand_t state, const padic_ctx_t ctx)
{
    long min, max;

    if (ctx->N > 0)
    {
        min = - ((ctx->N + 9) / 10);
        max = ctx->N;
    }
    else if (ctx->N < 0)
    {
        min = ctx->N - ((-ctx->N + 9) / 10);
        max = ctx->N;
    }
    else  /* ctx->N == 0 */
    {
        min = -10;
        max = 0;
    }

    padic_val(rop) = n_randint(state, max - min) + min;

    {
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(rop), ctx);
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

    for (i = 0; !_padic_is_zero(rop, ctx) && i < 10; i++)
        padic_randtest(rop, state, ctx);

    if (_padic_is_zero(rop, ctx))
    {
        fmpz_set_ui(padic_unit(rop), 1);
        padic_val(rop) = ctx->N - 1;
    }
}

