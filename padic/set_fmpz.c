#include "padic.h"

void _padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
{
    padic_val(rop) = fmpz_remove(padic_unit(rop), op, ctx->p);

    if (fmpz_sgn(op) < 0)
    {
        _padic_reduce_unit(rop, ctx);
    }
}

void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
{
    fmpz_set(padic_unit(rop), op);
    padic_val(rop) = 0;

    padic_normalise(rop, ctx);
}

