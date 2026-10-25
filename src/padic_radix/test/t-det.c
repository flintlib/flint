/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "gr.h"
#include "gr_mat.h"
#include "padic_radix.h"

static slong
_fmpz_val(const fmpz_t x, ulong p)
{
    fmpz_t t;
    slong v = 0;
    fmpz_init_set(t, x);
    while (!fmpz_is_zero(t) && fmpz_divisible_ui(t, p))
    {
        fmpz_divexact_ui(t, t, p);
        v++;
    }
    fmpz_clear(t);
    return v;
}

TEST_FUNCTION_START(padic_radix_det, state)
{
    slong iter;

    /* generic consistency (det(AB) = det(A) det(B)) */
    for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, 100));
        gr_mat_test_det((gr_method_mat_unary_op_get_scalar) padic_radix_mat_det, state, 10, 8, ctx);
        gr_mat_test_det((gr_method_mat_unary_op_get_scalar) padic_radix_mat_det_lu_full, state, 10, 8, ctx);
        gr_mat_test_det((gr_method_mat_unary_op_get_scalar) padic_radix_mat_det_levels, state, 10, 8, ctx);
        gr_ctx_clear(ctx);
    }

    /* matrices U diag(p^a_i) V known to O(p^N): the value must be correct
       to the claimed precision, and full pivoting attains the optimal
       precision N + min v(adj(A)) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        fmpz_mat_t U, V, M, Mi;
        fmpz_t det, den, t;
        gr_mat_t A;
        padic_radix_t d, e, z;
        ulong p;
        slong n, N, amax, i, j, opt, k;

        p = n_randint(state, 2) ? n_nth_prime(1 + n_randint(state, 5)) : n_randprime(state, 2 + n_randint(state, 30), 1);
        n = n_randint(state, n_randint(state, 4) ? 12 : 40);
        N = 1 + n_randint(state, 40);
        amax = n_randint(state, 4);

        gr_ctx_init_padic_radix(ctx, p, N + n * (amax + 1) * 64 + 10, PADIC_RADIX_PREC_INF, 0);

        fmpz_mat_init(U, n, n);
        fmpz_mat_init(V, n, n);
        fmpz_mat_init(M, n, n);
        fmpz_mat_init(Mi, n, n);
        fmpz_init(det);
        fmpz_init(den);
        fmpz_init(t);
        gr_mat_init(A, n, n, ctx);
        padic_radix_init(d, ctx);
        padic_radix_init(e, ctx);
        padic_radix_init(z, ctx);

        fmpz_mat_randbits(U, state, 8);
        fmpz_mat_randbits(V, state, 8);
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                fmpz_ui_pow_ui(t, p, n_randint(state, amax + 1));
                fmpz_mul(fmpz_mat_entry(U, i, j), fmpz_mat_entry(U, i, j), t);
            }
        fmpz_mat_mul(M, U, V);
        fmpz_mat_det(det, M);

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
            {
                padic_radix_struct * x = gr_mat_entry_ptr(A, i, j, ctx);
                GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(x, fmpz_mat_entry(M, i, j), ctx));
                x->N = N;
                GR_MUST_SUCCEED(_padic_radix_finalize(x, ctx));
            }

        GR_MUST_SUCCEED(padic_radix_exact_set_fmpz(e, det, ctx));

        for (k = 0; k < 4; k++)
        {
            int status;

            if (k == 0)
                status = padic_radix_mat_det_lu_full(d, A, ctx);
            else if (k == 1)
                status = padic_radix_mat_det(d, A, ctx);
            else if (k == 2)
                status = gr_mat_det_lu(d, A, ctx);   /* via the pivot method */
            else
                status = padic_radix_mat_det_levels(d, A, ctx);

            if (status == GR_UNABLE && k == 2)
                continue;

            GR_MUST_SUCCEED(status);
            GR_MUST_SUCCEED(padic_radix_sub(z, d, e, ctx));

            if (z->u.size != 0)
            {
                flint_printf("FAIL: value (k = %wd)\n", k);
                gr_ctx_println(ctx);
                flint_printf("d = "); gr_println(d, ctx);
                flint_printf("det = "); fmpz_print(det); flint_printf("\n");
                flint_abort();
            }

            if ((k == 0 || k == 3) && !fmpz_is_zero(det))
            {
                fmpz_mat_inv(Mi, den, M);
                opt = WORD_MAX;
                for (i = 0; i < n; i++)
                    for (j = 0; j < n; j++)
                        if (!fmpz_is_zero(fmpz_mat_entry(Mi, i, j)))
                        {
                            fmpz_mul(t, fmpz_mat_entry(Mi, i, j), det);
                            fmpz_divexact(t, t, den);
                            opt = FLINT_MIN(opt, _fmpz_val(t, p));
                        }

                /* the differential bound is the true precision when the
                   largest elementary divisor p^(a_n) is small enough
                   compared to p^N */
                if (n >= 1 && opt != WORD_MAX && _fmpz_val(det, p) - opt < N)
                {
                    opt = opt + N;

                    if (d->N != opt)
                    {
                        flint_printf("FAIL: precision %wd, optimal %wd\n", d->N, opt);
                        gr_ctx_println(ctx);
                        fmpz_mat_print_pretty(M); flint_printf("\n");
                        flint_printf("d = "); gr_println(d, ctx);
                        flint_abort();
                    }
                }
            }
        }

        fmpz_mat_clear(U);
        fmpz_mat_clear(V);
        fmpz_mat_clear(M);
        fmpz_mat_clear(Mi);
        fmpz_clear(det);
        fmpz_clear(den);
        fmpz_clear(t);
        gr_mat_clear(A, ctx);
        padic_radix_clear(d, ctx);
        padic_radix_clear(e, ctx);
        padic_radix_clear(z, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
