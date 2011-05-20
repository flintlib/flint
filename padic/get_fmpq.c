#include "padic.h"

void _padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(padic_unit(op)))
    {
        fmpq_set_si(rop, 0, 1);
    }
    else
    {
        fmpz_t pow;
        int alloc = 0;

        fmpz_set(fmpq_numref(rop), padic_unit(op));

        if (padic_val(op) == 0)
        {
            fmpz_set_ui(fmpq_denref(rop), 1);
        }
        else if (padic_val(op) > 0)
        {
            _padic_ctx_pow_ui(pow, &alloc, padic_val(op), ctx);
            fmpz_mul(fmpq_numref(rop), fmpq_numref(rop), pow);
            fmpz_set_ui(fmpq_denref(rop), 1);
        }
        else  /* padic_val(op) < 0 */
        {
            _padic_ctx_pow_ui(pow, &alloc, - padic_val(op), ctx);
            fmpz_set(fmpq_denref(rop), pow);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

void padic_get_fmpq(fmpq_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (padic_is_zero(op, ctx))
    {
        fmpq_set_si(rop, 0, 1);
    }
    else
    {
        fmpz_t pow;
        int alloc;

        _padic_ctx_pow_ui(pow, &alloc, ctx->N - padic_val(op), ctx);
        fmpz_mod(fmpq_numref(rop), padic_unit(op), pow);

        if (padic_val(op) == 0)
        {
            fmpz_set_ui(fmpq_denref(rop), 1);
        }
        else if (padic_val(op) > 0)
        {
            if (alloc)
                fmpz_clear(pow);
            _padic_ctx_pow_ui(pow, &alloc, padic_val(op), ctx);
            fmpz_mul(fmpq_numref(rop), fmpq_numref(rop), pow);
            fmpz_set_ui(fmpq_denref(rop), 1);
        }
        else  /* padic_val(op) < 0 */
        {
            if (alloc)
                fmpz_clear(pow);
            _padic_ctx_pow_ui(pow, &alloc, - padic_val(op), ctx);
            fmpz_set(fmpq_denref(rop), pow);
        }

        if (alloc)
            fmpz_clear(pow);
    }
}

