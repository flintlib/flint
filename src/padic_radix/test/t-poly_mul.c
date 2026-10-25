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
#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "padic_radix.h"

TEST_FUNCTION_START(padic_radix_poly_mul, state)
{
    slong iter;

    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx;
        slong len1, len2, nlo, nhi, i, k, e;
        padic_radix_struct * a, * b, * r1, * r2;
        int s1, s2, alg;

        gr_ctx_init_padic_radix_randtest(ctx, state, 1 + n_randint(state, n_randint(state, 4) ? 100 : 2000));
        e = radix_limb_exponent(PADIC_RADIX_CTX_RADIX(ctx));

        len1 = 1 + n_randint(state, n_randint(state, 4) ? 20 : 150);
        len2 = 1 + n_randint(state, n_randint(state, 4) ? 20 : 150);
        nhi = 1 + n_randint(state, len1 + len2 - 1);
        nlo = n_randint(state, nhi);

        a = gr_heap_init_vec(len1, ctx);
        b = gr_heap_init_vec(len2, ctx);
        r1 = gr_heap_init_vec(nhi, ctx);
        r2 = gr_heap_init_vec(nhi, ctx);

        for (i = 0; i < len1 + len2; i++)
        {
            padic_radix_struct * x = (i < len1) ? a + i : b + (i - len1);
            padic_radix_randtest(x, state, ctx);
            if (x->u.size != 0 && n_randint(state, 2))
            {
                slong rel = (x->N == PADIC_RADIX_EXACT) ? 0 : x->N - x->v;
                x->v = (slong) n_randint(state, 6 * e + 1) - 3 * e;
                if (x->N != PADIC_RADIX_EXACT)
                    x->N = x->v + rel;
                if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
                    GR_IGNORE(padic_radix_zero(x, ctx));
            }
        }

        /* multivariate-like structure: blocks of m coefficients with exact
           zero padding, exact leading ones, and block-dependent precisions */
        if (n_randint(state, 3) == 0)
        {
            slong m = 2 + n_randint(state, 8), shiftN = 0;
            for (i = 0; i < len1 + len2; i++)
            {
                padic_radix_struct * x = (i < len1) ? a + i : b + (i - len1);
                slong pos = (i < len1 ? i : i - len1) % m;
                if (pos == 0)
                    shiftN = n_randint(state, 3 * e + 1);
                if (pos == m - 1 || n_randint(state, 6) == 0)
                    GR_IGNORE(padic_radix_zero(x, ctx));
                else if (pos == m - 2)
                    GR_IGNORE(padic_radix_one(x, ctx));
                else if (x->u.size != 0 && x->N != PADIC_RADIX_EXACT && x->N - x->v > shiftN + 1)
                {
                    x->N -= shiftN;
                    if (_padic_radix_finalize(x, ctx) != GR_SUCCESS)
                        GR_IGNORE(padic_radix_zero(x, ctx));
                }
            }
        }

        /* the packed kernel gives exactly the delayed dot products (and the
           products, for coefficients with a single term) */
        s1 = _padic_radix_poly_mulmid_packed(r1, a, len1, b, len2, nlo, nhi, ctx);

        if (s1 != GR_UNABLE)
        {
            s2 = GR_SUCCESS;
            for (k = nlo; k < nhi; k++)
            {
                slong ilo = FLINT_MAX(0, k - len2 + 1), ihi = FLINT_MIN(k, len1 - 1);
                if (ihi == ilo)
                    s2 |= padic_radix_mul(r2 + k - nlo, a + ilo, b + k - ilo, ctx);
                else
                    s2 |= padic_radix_dot_strided_delayed(r2 + k - nlo, NULL, 0,
                        a + ilo, 1, b + k - ilo, -1, ihi - ilo + 1, ctx);
            }

            if (s1 == GR_SUCCESS && s2 == GR_SUCCESS)
                for (k = 0; k < nhi - nlo; k++)
                    if (r1[k].v != r2[k].v || r1[k].N != r2[k].N ||
                        !radix_integer_equal(&r1[k].u, &r2[k].u, PADIC_RADIX_CTX_RADIX(ctx)))
                    {
                        flint_printf("FAIL: packed mulmid vs delayed dot, coefficient %wd\n", k + nlo);
                        gr_ctx_println(ctx);
                        flint_printf("r1 = "); gr_println(r1 + k, ctx);
                        flint_printf("r2 = "); gr_println(r2 + k, ctx);
                        flint_abort();
                    }
        }

        /* Kronecker substitution and products modulo B^S: exactly the
           packed product (the min-plus convolutions are exact by default),
           and with PADIC_RADIX_FAST_BOUNDS never more precision, with the
           same values to that precision */
        for (alg = 0; alg < 3; alg++)
        {
            gr_ctx_t ctx2;
            gr_ctx_struct * cc = ctx;
            int fast = n_randint(state, 4) == 0;

            if (fast)
            {
                gr_ctx_init_padic_radix(ctx2, PADIC_RADIX_CTX_RADIX(ctx)->b.n,
                    PADIC_RADIX_CTX_PREC_REL(ctx), PADIC_RADIX_CTX_PREC_ABS(ctx),
                    PADIC_RADIX_CTX_FLAGS(ctx) | PADIC_RADIX_FAST_BOUNDS);
                cc = ctx2;
            }

            if (alg == 0)
                s1 = _padic_radix_poly_mulmid_KS(r1, a, len1, b, len2, nlo, nhi, cc);
            else if (alg == 1)
                s1 = _padic_radix_poly_mulmid_mod(r1, a, len1, b, len2, nlo, nhi, cc);
            else
                s1 = _padic_radix_poly_mulmid(r1, a, len1, b, len2, nlo, nhi, cc);

            if (s1 == GR_SUCCESS)
            {
                padic_radix_t d;
                padic_radix_init(d, ctx);
                s2 = _padic_radix_poly_mulmid_packed(r2, a, len1, b, len2, nlo, nhi, ctx);
                if (s2 == GR_SUCCESS)
                    for (k = 0; k < nhi - nlo; k++)
                    {
                        int bad;

                        if (alg == 2)
                        {
                            /* the dispatcher may use the classical algorithm,
                               whose dot products can differ in precision */
                            GR_IGNORE(padic_radix_sub(d, r1 + k, r2 + k, ctx));
                            bad = (d->u.size != 0 && d->v < FLINT_MIN(r1[k].N, r2[k].N));
                        }
                        else if (!fast)
                            bad = (r1[k].v != r2[k].v || r1[k].N != r2[k].N ||
                                !radix_integer_equal(&r1[k].u, &r2[k].u, PADIC_RADIX_CTX_RADIX(ctx)));
                        else
                        {
                            GR_IGNORE(padic_radix_sub(d, r1 + k, r2 + k, ctx));
                            bad = (r1[k].N > r2[k].N && r2[k].N != PADIC_RADIX_EXACT)
                                || (r1[k].N == PADIC_RADIX_EXACT && r2[k].N != PADIC_RADIX_EXACT)
                                || (d->u.size != 0 && d->v < FLINT_MIN(r1[k].N, r2[k].N));
                        }

                        if (bad)
                        {
                            flint_printf("FAIL: algorithm %d (fast = %d) vs packed, coefficient %wd\n", alg, fast, k + nlo);
                            gr_ctx_println(ctx);
                            flint_printf("r1 = "); gr_println(r1 + k, ctx);
                            flint_printf("r2 = "); gr_println(r2 + k, ctx);
                            flint_abort();
                        }
                    }
                padic_radix_clear(d, ctx);
            }

            if (fast)
                gr_ctx_clear(ctx2);
        }

        gr_heap_clear_vec(a, len1, ctx);
        gr_heap_clear_vec(b, len2, ctx);
        gr_heap_clear_vec(r1, nhi, ctx);
        gr_heap_clear_vec(r2, nhi, ctx);
        gr_ctx_clear(ctx);
    }

    TEST_FUNCTION_END(state);
}
