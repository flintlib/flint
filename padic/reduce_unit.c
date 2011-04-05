#include <stdlib.h>
#include <assert.h>

#include "padic.h"

void _padic_reduce_unit(padic_t rop, const padic_ctx_t ctx)
{
    fmpz_t pow;

    assert(rop[1] < ctx->N);

    /* TODO:  Re-write this using the powers in the context */

    fmpz_init(pow);
    fmpz_pow_ui(pow, ctx->p, ctx->N - rop[1]);
    fmpz_mod(rop, rop, pow);
    fmpz_clear(pow);
}

