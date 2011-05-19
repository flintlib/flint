#include "padic.h"

void padic_set_fmpq(padic_t rop, const fmpq_t op, const padic_ctx_t ctx)
{
    if (fmpq_is_zero(op))
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpq_t t;

        fmpq_init(t);

        padic_val(rop)  = fmpz_remove(fmpq_numref(t), fmpq_numref(op), ctx->p);
        padic_val(rop) -= fmpz_remove(fmpq_denref(t), fmpq_denref(op), ctx->p);

        if (padic_val(rop) >= ctx->N)
        {
            padic_zero(rop, ctx);
            fmpq_clear(t);
            return;
        }

        _padic_inv(fmpq_denref(t), fmpq_denref(t), ctx->p, ctx->N - padic_val(rop));
        fmpz_mul(padic_unit(rop), fmpq_numref(t), fmpq_denref(t));
        _padic_reduce_unit(rop, ctx);
        fmpq_clear(t);
    }
}

