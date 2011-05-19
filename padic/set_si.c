#include "padic.h"

void _padic_set_si(padic_t rop, long op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_si(t, op);
    _padic_set_fmpz(rop, t, ctx);
    fmpz_clear(t);
}

void padic_set_si(padic_t rop, long op, const padic_ctx_t ctx)
{
    fmpz_t t;

    fmpz_init(t);
    fmpz_set_si(t, op);
    padic_set_fmpz(rop, t, ctx);
    fmpz_clear(t);
}

