#include "padic.h"

void padic_pow_si(padic_t rop, const padic_t op, long e, const padic_ctx_t ctx)
{
    fmpz_t pow;
    
    if (padic_is_zero(op, ctx) || e * op[1] >= ctx->N)
    {
        padic_zero(rop, ctx);
        return;
    }
    
    if (e > 0)
    {
        /* Form u' = u^e mod p^{N - v'} */
        rop[1] = e * op[1];

        /* TODO:  Compute power */
        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, ctx->N - rop[1]);
        fmpz_powm_ui(rop, op, e, pow);
        fmpz_clear(pow);
    }
    else if (e < 0)
    {
        _padic_inv(rop, op, ctx->p, (ctx->N - op[1] * e + (-e - 1)) / -e);
        rop[1] = e * op[1];

        fmpz_init(pow);
        fmpz_pow_ui(pow, ctx->p, ctx->N - rop[1]);
        fmpz_powm_ui(rop, rop, -e, pow);
        fmpz_clear(pow);
    }
    else
    {
        padic_one(rop, ctx);
    }
}
