#include "padic.h"

void padic_ctx_init(padic_ctx_t ctx, const fmpz_t p, long N, 
                    enum padic_print_mode mode)
{
    if (!(N > 0))
    {
        printf("Exception:  N = %ld in padic_ctx_init.\n", N);
        abort();
    }


    fmpz_init(ctx->p);
    fmpz_set(ctx->p, p);

    ctx->N = N;

    ctx->pinv = (!COEFF_IS_MPZ(*p)) ? n_precompute_inverse(fmpz_get_ui(p)) : 0;

    ctx->mode = mode;
}

