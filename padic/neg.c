#include "padic.h"

void padic_neg(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t pow;

    if (padic_is_zero(op, ctx))
    {
        padic_zero(rop, ctx);
        return;
    }

    /* TODO:  Compute power */
    fmpz_init(pow);
    fmpz_pow_ui(pow, ctx->p, ctx->N - op[1]);
    fmpz_sub(rop, pow, op);
    fmpz_clear(pow);

    rop[1] = op[1];
}

