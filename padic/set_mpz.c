#include "padic.h"

void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
{
    fmpz_set_mpz(rop, op);
    rop[1] = 0;

    padic_normalise(rop, ctx);
}

