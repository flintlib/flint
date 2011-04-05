#include "padic.h"

void padic_set(padic_t rop, const padic_t op, const padic_ctx_t ctx)
{
    fmpz_set(rop, op);
    rop[1] = op[1];
}

