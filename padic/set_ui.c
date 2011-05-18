#include "padic.h"

void padic_set_ui(padic_t rop, ulong op, const padic_ctx_t ctx)
{
    if (op == 0)
    {
        padic_zero(rop, ctx);
    }
    else if (fmpz_cmp_ui(op, ctx->p) < 0)
    {
        fmpz_set_ui(rop, op);
        rop[1] = 0;
    }
    else
    {
        ulong p = fmpz_get_ui(ctx->p), q, r;

        /* Remove factors of p */
        rop[1] = 0;
        r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        while (r == 0)
        {
            op = q;
            rop[1] ++;
            r = n_divrem2_precomp(&q, op, p, ctx->pinv);
        }

        fmpz_set_ui(rop, op);

        _padic_reduce_unit(rop, ctx);
    }
}

