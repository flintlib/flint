#include "padic.h"

void padic_get_fmpz(fmpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(op))
    {
        fmpz_zero(rop);
    }
    else
    {
        if (op[1] == 0)
        {
            fmpz_set(rop, op);
        }
        else if (op[1] > 0)
        {
            fmpz_pow_ui(rop, ctx->p, op[1]);
            fmpz_mul(rop, rop, op);
        }
        else  /* op[1] < 0 */
        {
            fmpz_zero(rop);
        }
    }
}

