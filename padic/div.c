#include "padic.h"

void padic_div(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    padic_t inv;

    if (_padic_is_zero(op2))
    {
        printf("ERROR (padic_div).  op2 is zero.\n");
        abort();
    }

    if (padic_val(op1) - padic_val(op2) >= ctx->N)
    {
        padic_zero(rop, ctx);
        return;
    }

    padic_init(inv, ctx);

    _padic_inv(padic_unit(inv), padic_unit(op2), ctx->p, 
               ctx->N - padic_val(op1) + padic_val(op2));
    padic_val(inv) = - padic_val(op2);
    padic_mul(rop, op1, inv, ctx);

    padic_clear(inv, ctx);
}

