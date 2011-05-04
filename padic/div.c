#include "padic.h"

void padic_div(padic_t rop, const padic_t op1, const padic_t op2, 
               const padic_ctx_t ctx)
{
    padic_t inv;
    
    padic_init(inv, ctx);
    padic_inv(inv, op2, ctx);
    padic_mul(rop, op1, inv, ctx);
    padic_clear(inv, ctx);
}

