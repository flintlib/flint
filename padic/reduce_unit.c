#include <stdlib.h>
#include <assert.h>

#include "padic.h"

void _padic_reduce_unit(padic_t rop, const padic_ctx_t ctx)
{
    long e = ctx->N - rop[1];

    assert(e > 0);

    if (ctx->min <= e && e < ctx->max)
    {
        fmpz_mod(rop, rop, ctx->pow + (e - ctx->min));
    }
    else
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, e);
        fmpz_mod(rop, rop, pow);
        fmpz_clear(pow);
    }
}

