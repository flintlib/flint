#include <stdlib.h>
#include <assert.h>

#include "padic.h"

void _padic_reduce_unit(padic_t rop, const padic_ctx_t ctx)
{
    long e = ctx->N - padic_val(rop);

    assert(e > 0);

    if (ctx->min <= e && e < ctx->max)
    {
        fmpz_mod(padic_unit(rop), padic_unit(rop), ctx->pow + (e - ctx->min));
    }
    else
    {
        fmpz_t pow;

        fmpz_init(pow);
        fmpz_pow_ui(padic_unit(pow), ctx->p, e);
        fmpz_mod(padic_unit(rop), padic_unit(rop), pow);
        fmpz_clear(pow);
    }
}

