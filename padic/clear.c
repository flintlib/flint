#include "padic.h"

void padic_clear(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_clear(rop);
}

