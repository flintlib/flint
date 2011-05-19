#include "padic.h"

void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(op))
    {
        mpz_set_ui(rop, 0);
    }
    else
    {
        if (op[1] == 0)
        {
            fmpz_get_mpz(rop, op);
        }
        else if (op[1] > 0)
        {
            fmpz_t x;

            fmpz_init(x);
            fmpz_pow_ui(x, ctx->p, op[1]);
            fmpz_mul(x, x, op);
            fmpz_get_mpz(rop, x);
            fmpz_clear(x);
        }
        else  /* op[1] < 0 */
        {
            mpz_set_ui(rop, 0);
        }
    }
}

