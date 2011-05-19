#include "padic.h"

void _padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(padic_unit(op)) || (padic_val(op) < 0))
    {
        fmpz_zero(rop);
    }
    else
    {
        if (padic_val(op) == 0)
        {
            fmpz_set(rop, padic_unit(op));
        }
        else  /* (padic_val(op) > 0) */
        {
            fmpz_pow_ui(rop, ctx->p, padic_val(op));
            fmpz_mul(rop, rop, padic_unit(op));
        }
    }
}

void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx) || (padic_val(op) < 0))
    {
        fmpz_zero(rop);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op), ctx);
        fmpz_mod(rop, padic_unit(op), pow);

        if (padic_val(op) > 0)
        {
            if (alloc)
                fmpz_clear(pow);
            _padic_ctx_pow_ui(pow, &alloc, padic_val(op), ctx);
            fmpz_mul(rop, rop, pow);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

