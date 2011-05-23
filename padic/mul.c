#include "padic.h"

void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (_padic_is_zero(op1) || _padic_is_zero(op2))
    {
        padic_zero(rop, ctx);
        return;
    }

    padic_val(rop) = padic_val(op1) + padic_val(op2);

    if (padic_val(rop) >= ctx->N)
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_mul(padic_unit(rop), padic_unit(op1), padic_unit(op2));
        _padic_reduce_unit(rop, ctx);
    }
}

