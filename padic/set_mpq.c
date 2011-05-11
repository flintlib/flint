#include "padic.h"

void padic_set_mpq(padic_t rop, const mpq_t op, const padic_ctx_t ctx)
{

    if (mpz_sgn(mpq_numref(op)) == 0)
    {
        padic_zero(rop, ctx);
    }
    else
    {
        fmpz_t num, den;

        fmpz_init(num);
        fmpz_init(den);
        fmpz_set_mpz(num, mpq_numref(op));
        fmpz_set_mpz(den, mpq_denref(op));

        rop[1]  = fmpz_remove(num, num, ctx->p);
        rop[1] -= fmpz_remove(den, den, ctx->p);

        if (rop[1] >= ctx->N)
        {
            padic_zero(rop, ctx);
            fmpz_clear(num);
            fmpz_clear(den);
            return;
        }

        _padic_inv(den, den, ctx->p, ctx->N - rop[1]);
        fmpz_mul(rop, num, den);
        _padic_reduce_unit(rop, ctx);
        fmpz_clear(num);
        fmpz_clear(den);
    }
}

