#include "padic.h"

void padic_shift(padic_t rop, const padic_t op, long v, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op) || (padic_val(op) + v >= ctx->N))
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_set(padic_unit(rop), padic_unit(op));
        padic_val(rop) = padic_val(op) + v;
        _padic_reduce_unit(rop, ctx);
    }
}

