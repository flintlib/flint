#include "padic.h"

void padic_normalise(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(rop) && (rop[1] < ctx->N))
    {
        rop[1] += _fmpz_remove(rop, ctx->p, ctx->pinv);

        if (rop[1] < ctx->N)
            _padic_reduce_unit(rop, ctx);
        else
            padic_zero(rop, ctx);
    }
    else
        padic_zero(rop, ctx);
}

