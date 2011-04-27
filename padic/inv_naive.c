#include "padic.h"

void _padic_inv_naive(fmpz_t rop, const fmpz_t op, const fmpz_t p, long N)
{
    fmpz_t pow;

    fmpz_init(pow);
    fmpz_pow_ui(pow, p, N);
    fmpz_invmod(rop, op, pow);
    fmpz_clear(pow);
}

void padic_inv_naive(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;
    
    if (padic_is_zero(op, ctx))
    {
        printf("Exception (padic_inv_naive).  Zero is not invertible.\n");
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

    _padic_inv_naive(rop, op, ctx->p, ctx->N + op[1]);

    rop[1] = -op[1];
}

