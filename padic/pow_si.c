#include "padic.h"

void padic_pow_si(padic_t rop, const padic_t op, long e, const padic_ctx_t ctx)
{
    if (_padic_is_zero(op) || e * padic_val(op) >= ctx->N)
    {
        padic_zero(rop, ctx);
        return;
    }
    
    if (e > 0)
    {
        fmpz_t pow;
        int alloc;

        /* Form u' = u^e mod p^{N - v'} */
        padic_val(rop) = e * padic_val(op);

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(rop), ctx);
        fmpz_powm_ui(padic_unit(rop), padic_unit(op), e, pow);
        if (alloc)
            fmpz_clear(pow);
    }
    else if (e < 0)
    {
        fmpz_t pow;
        int alloc;

        _padic_inv(padic_unit(rop), padic_unit(op), 
                   ctx->p, (ctx->N - padic_val(op) * e + (-e - 1)) / -e);
        padic_val(rop) = e * padic_val(op);

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(rop), ctx);
        fmpz_powm_ui(padic_unit(rop), padic_unit(rop), -e, pow);
        if (alloc)
            fmpz_clear(pow);
    }
    else
    {
        padic_one(rop, ctx);
    }
}
