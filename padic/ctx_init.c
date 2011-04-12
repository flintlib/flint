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

    ctx->min = FLINT_MAX(1, 9 * N / 10);
    ctx->max = FLINT_MAX(N + 1, 11 * N / 10);

    {
        long i, len = ctx->max - ctx->min;

        ctx->pow = _fmpz_vec_init(len);

        fmpz_pow_ui(ctx->pow, p, ctx->min);
        for (i = 1; i < len; i++)
            fmpz_mul(ctx->pow + i, ctx->pow + (i - 1), p);
    }

    ctx->mode = mode;
}

