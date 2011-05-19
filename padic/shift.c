#include "padic.h"

void padic_shift(padic_t rop, const padic_t op, long v, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(op) || (op[1] + v >= ctx->N))
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_set(rop, op);
        rop[1] = op[1] + v;
        _padic_reduce_unit(rop, ctx);
    }
}

