#include "padic.h"

void padic_ctx_clear(padic_ctx_t ctx)
{
    fmpz_clear(ctx->p);
    _fmpz_vec_clear(ctx->pow, ctx->max - ctx->min);
}

