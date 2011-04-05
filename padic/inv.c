#include "padic.h"

void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;
    
    if (padic_is_zero(op, ctx))
    {
        printf("Exception (padic_inv).  Zero is not invertible.\n");
        abort();
    }

    /* op has no meaningful inverse to the given precision.. */
    /* TODO:  Decide what to do in this case */
    if (ctx->N + op[1] <= 0)
    {
        padic_zero(rop, ctx);
        return;
    }
    
    /* Unit part */
    fmpz_init(pow);
    fmpz_pow_ui(pow, ctx->p, ctx->N + op[1]);
    if (op[1] >= 0)
    {
        fmpz_invmod(rop, op, pow);
    }
    else
    {
        fmpz_mod(rop, op, pow);
        fmpz_invmod(rop, rop, pow);
    }
    fmpz_clear(pow);
    
    /* Valuation part */
    rop[1] = -op[1];
}

