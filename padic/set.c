#include "padic.h"

void _padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_set(padic_unit(rop), padic_unit(op));
    padic_val(rop) = padic_val(op);
}

void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op), ctx);
        fmpz_mod(padic_unit(rop), padic_unit(op), pow);
        padic_val(rop) = padic_val(op);
        if (alloc)
            fmpz_clear(pow);
    }
}

