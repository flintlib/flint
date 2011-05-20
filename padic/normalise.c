#include "padic.h"

void padic_normalise(padic_t rop, const padic_ctx_t ctx)
{
    if (!fmpz_is_zero(padic_unit(rop)) && (padic_val(rop) < ctx->N))
    {
        padic_val(rop) += _fmpz_remove(padic_unit(rop), ctx->p, ctx->pinv);

        if (padic_val(rop) < ctx->N)
            _padic_reduce_unit(rop, ctx);
        else
            padic_zero(rop, ctx);
    }
    else
        padic_zero(rop, ctx);
}

