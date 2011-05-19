#include "padic.h"

void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;
    int alloc;

    if (padic_is_zero(op, ctx))
    {
        padic_zero(rop, ctx);
        return;
    }

    _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op), ctx);
    fmpz_sub(padic_unit(rop), pow, padic_unit(op));
    if (alloc)
        fmpz_clear(pow);

    padic_val(rop) = padic_val(op);
}

