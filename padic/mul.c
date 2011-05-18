#include "padic.h"

void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op1, ctx) || padic_is_zero(op2, ctx))
    {
        padic_zero(rop, ctx);
        return;
    }

    rop[1] = op1[1] + op2[1];

    if (rop[1] >= ctx->N)
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_mul(rop, op1, op2);
        _padic_reduce_unit(rop, ctx);
    }
}

