#include "padic.h"

/* TODO:  Improve this along the lines of set_ui */

void padic_set_si(padic_t rop, long op, const padic_ctx_t ctx)
{
    fmpz_t x;

    fmpz_init(x);
    fmpz_set_si(x, op);

    padic_set_fmpz(rop, x, ctx);

    fmpz_clear(x);
}

