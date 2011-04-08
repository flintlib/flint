#include "padic.h"

void padic_inv(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;
    
    if (padic_is_zero(op, ctx))
    {
        printf("Exception (padic_inv).  Zero is not invertible.\n");
        abort();
    }

    /*
        If x = u p^v has negative valuation with N <= -v 
        then there is no inverse of x defined modulo p^N.
     */
    if (ctx->N + op[1] <= 0)
    {
        padic_zero(rop, ctx);
        return;
    }

    /* Unit part */
    fmpz_init(pow);
    if (op[1] <= 0)
    {
        fmpz_pow_ui(pow, ctx->p, ctx->N + op[1]);
        fmpz_invmod(rop, op, pow);
    }
    else
    {
        fmpz_pow_ui(pow, ctx->p, ctx->N - op[1]);
        fmpz_mod(rop, op, pow);
        fmpz_invmod(rop, rop, pow);
    }
    fmpz_clear(pow);
    
    /* Valuation part */
    rop[1] = -op[1];
}

