#include "padic.h"

void _padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_ui(t, op);
    _padic_set_fmpz(rop, t, ctx);
    fmpz_clear(t);
}

void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
{
    if (op == 0)
    {
        padic_zero(rop, ctx);
    }
    else if (fmpz_cmp_ui(ctx->p, op) > 0)
    {
        fmpz_set_ui(padic_unit(rop), op);
        padic_val(rop) = 0;
    }
    else
    {
        ulong p = fmpz_get_ui(ctx->p), q, r;

        /* Remove factors of p */
        padic_val(rop) = 0;
        r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        while (r == 0)
        {
            op = q;
            padic_val(rop)++;
            r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        }

        fmpz_set_ui(padic_unit(rop), op);
        _padic_reduce_unit(rop, ctx);
    }
}

