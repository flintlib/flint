#include "padic.h"

void _padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_mpz(t, op);
    _padic_set_fmpz(rop, t, ctx);
    fmpz_clear(t);
}

void padic_set_mpz(padic_t rop, const mpz_t op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_mpz(t, op);
    padic_set_fmpz(rop, t, ctx);
    fmpz_clear(t);
}

