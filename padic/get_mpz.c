#include "padic.h"

void padic_get_mpz(mpz_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    padic_get_fmpz(t, op, ctx);
    fmpz_get_mpz(rop, t);
    fmpz_clear(t);
}

