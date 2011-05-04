#include "padic.h"

void padic_init(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_init(rop);
    rop[1] = 0;
}

