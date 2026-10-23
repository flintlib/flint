/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "fq_zech.h"
#include "fq_zech_vec.h"
#include "fq_zech_poly.h"

TEST_FUNCTION_START(fq_zech_poly_mulmid_univariate, state)
{
    slong iter;

    for (iter = 0; iter < 400 * flint_test_multiplier(); iter++)
    {
        fq_zech_ctx_t ctx;
        fq_zech_struct *a, *b, *c, *d;
        slong len1, len2, nlo, nhi, i, maxlen;
        int squaring;

        switch (n_randint(state, 8))
        {
            case 0:
                fq_zech_ctx_init_random_ui(ctx, 2, 1 + n_randint(state, 12), "a");
                break;
            case 1:
                fq_zech_ctx_init_random_ui(ctx, 3, 1 + n_randint(state, 7), "a");
                break;
            case 2:
                fq_zech_ctx_init_random_ui(ctx, n_randprime(state, 7, 1), 1 + n_randint(state, 2), "a");
                break;
            case 3:
                fq_zech_ctx_init_random_ui(ctx, n_randprime(state, 16, 1), 1, "a");
                break;
            default:
                fq_zech_ctx_init_randtest(ctx, state, n_randint(state, 4));
                break;
        }

        maxlen = n_randint(state, 5) ? 30 : 300;
        len1 = 1 + n_randint(state, maxlen);
        len2 = 1 + n_randint(state, maxlen);
        squaring = n_randint(state, 5) == 0;
        if (squaring)
            len2 = len1;

        nhi = 1 + n_randint(state, len1 + len2 - 1);
        nlo = n_randint(state, 4) ? n_randint(state, nhi) : 0;

        a = _fq_zech_vec_init(len1, ctx);
        b = _fq_zech_vec_init(len2, ctx);
        c = _fq_zech_vec_init(nhi, ctx);
        d = _fq_zech_vec_init(nhi, ctx);

        for (i = 0; i < len1; i++)
        {
            if (n_randint(state, 4) == 0)
                fq_zech_zero(a + i, ctx);
            else
                fq_zech_rand(a + i, state, ctx);
        }

        for (i = 0; i < len2; i++)
        {
            if (n_randint(state, 4) == 0)
                fq_zech_zero(b + i, ctx);
            else
                fq_zech_rand(b + i, state, ctx);
        }

        /* zero padding */
        if (n_randint(state, 4) == 0)
            for (i = n_randint(state, len1); i < len1; i++)
                fq_zech_zero(a + i, ctx);
        if (n_randint(state, 4) == 0)
            for (i = n_randint(state, len2); i < len2; i++)
                fq_zech_zero(b + i, ctx);

        if (squaring)
        {
            _fq_zech_poly_mullow_classical(d, a, len1, a, len1, nhi, ctx);
            _fq_zech_poly_mulmid_univariate(c, a, len1, a, len1, nlo, nhi, ctx);
        }
        else
        {
            _fq_zech_poly_mullow_classical(d, a, len1, b, len2, nhi, ctx);
            _fq_zech_poly_mulmid_univariate(c, a, len1, b, len2, nlo, nhi, ctx);
        }

        if (!_fq_zech_vec_equal(c, d + nlo, nhi - nlo, ctx))
        {
            flint_printf("FAIL\n");
            flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd, squaring = %d\n", len1, len2, nlo, nhi, squaring);
            fq_zech_ctx_print(ctx);
            flint_abort();
        }

        /* aliasing */
        if (nlo == 0 && !squaring)
        {
            fq_zech_struct * t = _fq_zech_vec_init(FLINT_MAX(len1, nhi), ctx);

            _fq_zech_vec_set(t, a, len1, ctx);
            _fq_zech_poly_mulmid_univariate(t, t, len1, b, len2, 0, nhi, ctx);

            if (!_fq_zech_vec_equal(t, d, nhi, ctx))
            {
                flint_printf("FAIL (aliasing)\n");
                flint_printf("len1 = %wd, len2 = %wd, nhi = %wd\n", len1, len2, nhi);
                fq_zech_ctx_print(ctx);
                flint_abort();
            }

            _fq_zech_vec_clear(t, FLINT_MAX(len1, nhi), ctx);
        }

        _fq_zech_vec_clear(a, len1, ctx);
        _fq_zech_vec_clear(b, len2, ctx);
        _fq_zech_vec_clear(c, nhi, ctx);
        _fq_zech_vec_clear(d, nhi, ctx);

        fq_zech_ctx_clear(ctx);
    }

    /* Fixed fields with 2^17 <= q - 1 < 2^18 (a separate range of the
       crossover between classical and univariate multiplication). The
       tables are built once per field. */
    {
        const ulong ps[2] = { 409, 3 };
        const slong ds[2] = { 2, 11 };
        slong k;

        for (k = 0; k < 2; k++)
        {
            fq_zech_ctx_t ctx;
            fq_zech_poly_t A, B, C, D;

            fq_zech_ctx_init_random_ui(ctx, ps[k], ds[k], "a");

            if (ctx->qm1 < (UWORD(1) << 17) || ctx->qm1 >= (UWORD(1) << 18))
            {
                flint_printf("FAIL (field size)\n");
                flint_abort();
            }

            fq_zech_poly_init(A, ctx);
            fq_zech_poly_init(B, ctx);
            fq_zech_poly_init(C, ctx);
            fq_zech_poly_init(D, ctx);

            for (iter = 0; iter < 10 * flint_test_multiplier(); iter++)
            {
                slong len1, len2, nlo, nhi, n;
                fq_zech_struct * c;

                len1 = 1 + n_randint(state, 80);
                len2 = 1 + n_randint(state, 80);

                fq_zech_poly_randtest_not_zero(A, state, len1, ctx);
                fq_zech_poly_randtest_not_zero(B, state, len2, ctx);
                len1 = A->length;
                len2 = B->length;

                /* full product, through the automatic algorithm selection */
                fq_zech_poly_mul(C, A, B, ctx);
                fq_zech_poly_mul_classical(D, A, B, ctx);

                if (!fq_zech_poly_equal(C, D, ctx))
                {
                    flint_printf("FAIL (fixed field, mul)\n");
                    flint_printf("len1 = %wd, len2 = %wd\n", len1, len2);
                    fq_zech_ctx_print(ctx);
                    flint_abort();
                }

                /* truncated product */
                n = 1 + n_randint(state, len1 + len2 - 1);
                fq_zech_poly_mullow(C, A, B, n, ctx);
                fq_zech_poly_mullow_classical(D, A, B, n, ctx);

                if (!fq_zech_poly_equal(C, D, ctx))
                {
                    flint_printf("FAIL (fixed field, mullow)\n");
                    flint_printf("len1 = %wd, len2 = %wd, n = %wd\n", len1, len2, n);
                    fq_zech_ctx_print(ctx);
                    flint_abort();
                }

                /* squaring */
                fq_zech_poly_sqr(C, A, ctx);
                fq_zech_poly_mul_classical(D, A, A, ctx);

                if (!fq_zech_poly_equal(C, D, ctx))
                {
                    flint_printf("FAIL (fixed field, sqr)\n");
                    flint_printf("len1 = %wd\n", len1);
                    fq_zech_ctx_print(ctx);
                    flint_abort();
                }

                /* middle product computed directly */
                nhi = 1 + n_randint(state, len1 + len2 - 1);
                nlo = n_randint(state, nhi);

                fq_zech_poly_mullow_classical(D, A, B, nhi, ctx);
                fq_zech_poly_fit_length(D, nhi, ctx);
                for (n = D->length; n < nhi; n++)
                    fq_zech_zero(D->coeffs + n, ctx);

                c = _fq_zech_vec_init(nhi - nlo, ctx);
                _fq_zech_poly_mulmid_univariate(c, A->coeffs, len1, B->coeffs, len2, nlo, nhi, ctx);

                if (!_fq_zech_vec_equal(c, D->coeffs + nlo, nhi - nlo, ctx))
                {
                    flint_printf("FAIL (fixed field, mulmid)\n");
                    flint_printf("len1 = %wd, len2 = %wd, nlo = %wd, nhi = %wd\n", len1, len2, nlo, nhi);
                    fq_zech_ctx_print(ctx);
                    flint_abort();
                }

                _fq_zech_vec_clear(c, nhi - nlo, ctx);
            }

            fq_zech_poly_clear(A, ctx);
            fq_zech_poly_clear(B, ctx);
            fq_zech_poly_clear(C, ctx);
            fq_zech_poly_clear(D, ctx);
            fq_zech_ctx_clear(ctx);
        }
    }

    TEST_FUNCTION_END(state);
}
