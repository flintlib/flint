#include "padic.h"

void padic_get_mpq(mpq_t rop, const padic_t op, const padic_ctx_t ctx)
{
    if (fmpz_is_zero(op))
    {
        mpq_set_si(rop, 0, 1);
    }
    else
    {
        if (op[1] == 0)
        {
            fmpz_get_mpz(mpq_numref(rop), op);
            mpz_set_ui(mpq_denref(rop), 1);
        }
        else if (op[1] > 0)
        {
            fmpz_t x;

            fmpz_init(x);
            fmpz_pow_ui(x, ctx->p, op[1]);
            fmpz_mul(x, x, op);
            fmpz_get_mpz(mpq_numref(rop), x);
            mpz_set_ui(mpq_denref(rop), 1);
            fmpz_clear(x);
        }
        else  /* op[1] < 0 */
        {
            fmpz_t x;

            fmpz_init(x);
            fmpz_pow_ui(x, ctx->p, -op[1]);
            fmpz_get_mpz(mpq_numref(rop), op);
            fmpz_get_mpz(mpq_denref(rop), x);
            fmpz_clear(x);
        }
    }
}

