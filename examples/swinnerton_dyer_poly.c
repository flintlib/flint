/* This file is public domain. Author: Fredrik Johansson. */

#include "flint/profiler.h"
#include "ca.h"
#include "ca_vec.h"

void
_ca_poly_mullow(ca_ptr res, ca_srcptr x, slong xlen,
    ca_srcptr y, slong ylen, slong len, ca_ctx_t ctx)
{
    slong i, j;
    ca_t t;

    for (i = 0; i < len; i++)
        ca_zero(res + i, ctx);

    ca_init(t, ctx);

    for (i = 0; i < xlen; i++)
    {
        for (j = 0; j < FLINT_MIN(ylen, len - i); j++)
        {
            ca_mul(t, x + i, y + j, ctx);
            ca_add(res + i + j, res + i + j, t, ctx);
        }
    }

    ca_clear(t, ctx);
}

void
swinnerton_dyer_poly(ca_ptr T, ulong n, slong trunc, ca_ctx_t ctx)
{
    ca_ptr square_roots, tmp1, tmp2, tmp3;
    ca_t one;
    slong i, j, k, N;

    N = WORD(1) << n;
    trunc = FLINT_MIN(trunc, N + 1);

    ca_init(one, ctx);
    ca_one(one, ctx);

    square_roots = _ca_vec_init(n, ctx);
    tmp1 = flint_malloc((N/2 + 1) * sizeof(ca_struct));
    tmp2 = flint_malloc((N/2 + 1) * sizeof(ca_struct));
    tmp3 = _ca_vec_init(N, ctx);

    for (i = 0; i < n; i++)
        ca_sqrt_ui(square_roots + i, n_nth_prime(i + 1), ctx);

    /* Build linear factors */
    for (i = 0; i < N; i++)
    {
        ca_zero(T + i, ctx);

        for (j = 0; j < n; j++)
        {
            if ((i >> j) & 1)
                ca_add(T + i, T + i, square_roots + j, ctx);
            else
                ca_sub(T + i, T + i, square_roots + j, ctx);
        }
    }

    /* For each level... */
    for (i = 0; i < n; i++)
    {
        slong stride = UWORD(1) << i;

        for (j = 0; j < N; j += 2*stride)
        {
            for (k = 0; k < stride; k++)
            {
                tmp1[k] = T[j + k];
                tmp2[k] = T[j + stride + k];
            }

            tmp1[stride] = *one;
            tmp2[stride] = *one;

            _ca_poly_mullow(tmp3, tmp1, stride + 1, tmp2, stride + 1,
                FLINT_MIN(2 * stride, trunc), ctx);
            _ca_vec_set(T + j, tmp3, FLINT_MIN(2 * stride, trunc), ctx);
        }
    }

    ca_one(T + N, ctx);
    _ca_vec_clear(square_roots, n, ctx);
    flint_free(tmp1);
    flint_free(tmp2);
    _ca_vec_clear(tmp3, UWORD(1) << n, ctx);
    ca_clear(one, ctx);
}

int main(int argc, char *argv[])
{
    ca_ctx_t ctx;
    ca_ptr poly;
    slong i, n, N;

    if (argc < 2)
    {
        flint_printf("usage: build/examples/sdpoly n\n");
        return 1;
    }

    n = atol(argv[1]);
    if (n < 0 || n > 20)
        flint_abort();
    N = (1 << n);

    TIMEIT_ONCE_START
    ca_ctx_init(ctx);
    poly = _ca_vec_init(N + 1, ctx);
    swinnerton_dyer_poly(poly, n, N + 1, ctx);

    for (i = 0; i <= N; i++)
    {
        ca_print(poly + i, ctx);
        flint_printf("\n");
    }

    _ca_vec_clear(poly, N + 1, ctx);
    ca_ctx_clear(ctx);

    flint_printf("\n");
    TIMEIT_ONCE_STOP
    SHOW_MEMORY_USAGE

    flint_cleanup();
    return EXIT_SUCCESS;
}
