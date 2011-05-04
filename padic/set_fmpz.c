#include "padic.h"

void padic_set_fmpz(padic_t rop, const fmpz_t op, const padic_ctx_t ctx)
{
    fmpz_set(rop, op);
    rop[1] = 0;

    padic_normalise(rop, ctx);
}

