#include "padic.h"

/*
    Computes $u_1 p^{v_1} * u_2 p^{v_2} = (u_1 u_2) p^{v_1 + v_2}$.

    Assumes that \code{op1} , \code{op2} are non-zero. 
 */

void _padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
                const padic_ctx_t ctx)
{
    rop[1] = op1[1] + op2[1];

    if (rop[1] >= ctx->N)
    {
        padic_zero(rop, ctx);
        return;
    }

    fmpz_mul(rop, op1, op2);
}

void padic_mul(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    if (padic_is_zero(op1, ctx) || padic_is_zero(op2, ctx))
    {
        padic_zero(rop, ctx);
        return;
    }

    _padic_mul(rop, op1, op2, ctx);
    _padic_reduce_unit(rop, ctx);
}

