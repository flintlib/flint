#include "padic.h"

void padic_normalise(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(rop) || rop[1] < ctx->N)
    {
        rop[1] += fmpz_remove(rop, rop, ctx->p);

        if (rop[1] < ctx->N)
            _padic_reduce_unit(rop, ctx);
        else
            padic_zero(rop);
    }
    else
        padic_zero(rop);
}

