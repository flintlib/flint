/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "test_helpers.h"
#include "arf.h"
#include "gr.h"
#include "gr_vec.h"
#include "gr_special.h"
#include "dfloat.h"

/* ulp-nonoverlapping: |x[i+1]| <= ulp(x[i]), zeros only at the end
   (a nonfinite head is accepted as is) */
static int
is_canonical(const double * x, int n)
{
    int i;
    if (!(fabs(x[0]) <= DBL_MAX))
        return 1;
    for (i = 0; i + 1 < n; i++)
    {
        if (x[i] == 0.0)
        {
            for (; i < n; i++)
                if (x[i] != 0.0)
                    return 0;
            return 1;
        }
        if (fabs(x[i + 1]) > ldexp(1.0, ilogb(x[i]) - 52))
            return 0;
    }
    return 1;
}

TEST_FUNCTION_START(canonical, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    /* canonicalise: exact and canonical */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        int n = 1 + n_randint(state, DFLOAT_MAX_N), i;
        double x[DFLOAT_MAX_N], y[DFLOAT_MAX_N];
        arf_t a, b;

        /* valid (weakly nonoverlapping) expansions: random ones, and
           the results of the kernels on random ones */
        _dfloat_randtest(x, n, state);
        if (n_randint(state, 2))
        {
            double z[DFLOAT_MAX_N], e = 0.0;
            _dfloat_randtest(z, n, state);
            if (_dfloat_abs_sum(x, n) < 1e100 && _dfloat_abs_sum(z, n) < 1e100)
            {
                switch (n_randint(state, 3))
                {
                    case 0: _dfloat_add_n(y, &e, n, x, z); break;
                    case 1: _dfloat_sub_n(y, &e, n, x, z); break;
                    default: _dfloat_mul_n(y, &e, n, x, z); break;
                }
                for (i = 0; i < n; i++)
                    x[i] = y[i];
            }
        }
        for (i = 0; i < n; i++)
            if (!(fabs(x[i]) < 1e300))
                x[i] = 0.0;

        _dfloat_canonicalise(y, n, x);

        arf_init(a); arf_init(b);
        _dfloat_get_arf(a, x, n);
        _dfloat_get_arf(b, y, n);
        if (!arf_equal(a, b) || !is_canonical(y, n))
        {
            flint_printf("FAIL: canonicalise (n = %d)\n", n);
            for (i = 0; i < n; i++) flint_printf("%.17g ", x[i]);
            flint_printf("\n");
            for (i = 0; i < n; i++) flint_printf("%.17g ", y[i]);
            flint_printf("\n");
            flint_abort();
        }
        arf_clear(a); arf_clear(b);
    }

    /* contrived inputs: interior zeros, overlapping and unsorted tails */
    {
        const double cases[][4] = {
            {1.0, 0.0, 0x1p-120, 0.0},
            {1.0, 0.0, 0.0, 0x1p-170},
            {1.0, 0x1p-110, 0x1p-111, 0x1p-112},
            {1.0, 0x1p-53, 0x1p-53, 0x1p-53},
            {1.0, -0x1p-53, 0x1p-106, -0x1p-159},
            {1.0, 0x1p-100, 0x1p-60, 0x1p-160},
            {0.0, 0x1p-100, 0.0, 0x1p-160},
            {1.0, 0x1p-52, 0x1p-104, 0x1p-156},
            {0x1p-1000, 0x1p-1060, 0x1p-1074, 0.0},
            {3.0, 0x1.8p-52, 0x1.8p-104, 0x1.8p-156},
        };
        slong c;
        int n, i;
        for (c = 0; c < (slong) (sizeof(cases) / sizeof(cases[0])); c++)
        {
            for (n = 1; n <= DFLOAT_MAX_N; n++)
            {
                double y[DFLOAT_MAX_N];
                arf_t a, b;
                _dfloat_canonicalise(y, n, cases[c]);
                arf_init(a); arf_init(b);
                _dfloat_get_arf(a, cases[c], n);
                _dfloat_get_arf(b, y, n);
                if (!arf_equal(a, b) || !is_canonical(y, n))
                {
                    flint_printf("FAIL: canonicalise case %wd (n = %d)\n", c, n);
                    for (i = 0; i < n; i++) flint_printf("%.17g ", y[i]);
                    flint_printf("\n");
                    flint_abort();
                }
                arf_clear(a); arf_clear(b);
            }
        }
    }

    /* strong contexts: same values as the weak ones, canonical form */
    for (iter = 0; iter < 3000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t wctx, sctx;
        int n = 1 + n_randint(state, DFLOAT_MAX_N), ball = n_randint(state, 2), which, i;
        slong len = n_randint(state, 12);
        gr_ptr x, y, rw, rs;
        arf_t a, b;
        int sw, ss;

        GR_MUST_SUCCEED(gr_ctx_init_dfloat(wctx, n, ball));
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(sctx, n, ball | DFLOAT_STRONG));

        x = gr_heap_init_vec(len + 1, wctx);
        y = gr_heap_init_vec(len + 1, wctx);
        rw = gr_heap_init_vec(len + 1, wctx);
        rs = gr_heap_init_vec(len + 1, wctx);
        for (i = 0; i <= len; i++)
        {
            GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(x, i, wctx->sizeof_elem), state, wctx));
            GR_MUST_SUCCEED(gr_randtest(GR_ENTRY(y, i, wctx->sizeof_elem), state, wctx));
        }

        which = n_randint(state, 10);
        switch (which)
        {
            case 0: sw = gr_add(rw, x, y, wctx); ss = gr_add(rs, x, y, sctx); break;
            case 1: sw = gr_sub(rw, x, y, wctx); ss = gr_sub(rs, x, y, sctx); break;
            case 2: sw = gr_mul(rw, x, y, wctx); ss = gr_mul(rs, x, y, sctx); break;
            case 3: sw = gr_sqr(rw, x, wctx); ss = gr_sqr(rs, x, sctx); break;
            case 4: sw = gr_div(rw, x, y, wctx); ss = gr_div(rs, x, y, sctx); break;
            case 5: sw = gr_sqrt(rw, x, wctx); ss = gr_sqrt(rs, x, sctx); break;
            case 6: sw = gr_exp(rw, x, wctx); ss = gr_exp(rs, x, sctx); break;
            case 7: sw = _gr_vec_mul(rw, x, y, len, wctx); ss = _gr_vec_mul(rs, x, y, len, sctx); break;
            case 8: sw = _gr_vec_add(rw, x, y, len, wctx); ss = _gr_vec_add(rs, x, y, len, sctx); break;
            default: sw = _gr_vec_dot(rw, x, 0, x, y, len, wctx); ss = _gr_vec_dot(rs, x, 0, x, y, len, sctx); break;
        }

        if (sw != ss)
        {
            flint_printf("FAIL: status (n = %d, ball = %d, which = %d)\n", n, ball, which);
            flint_abort();
        }

        arf_init(a); arf_init(b);
        for (i = 0; i < ((which >= 7 && which <= 8) ? len : 1); i++)
        {
            const double * dw = GR_ENTRY(rw, i, wctx->sizeof_elem);
            const double * ds = GR_ENTRY(rs, i, wctx->sizeof_elem);
            _dfloat_get_arf(a, dw, n);
            _dfloat_get_arf(b, ds, n);
            if (sw == GR_SUCCESS && (!(arf_equal(a, b) || (arf_is_nan(a) && arf_is_nan(b)))
                    || (ball && !(dw[n] == ds[n] || (dw[n] != dw[n] && ds[n] != ds[n])))
                    || !is_canonical(ds, n)))
            {
                flint_printf("FAIL: strong (n = %d, ball = %d, which = %d, i = %wd)\n", n, ball, which, i);
                for (i = 0; i < n + ball; i++) flint_printf("%.17g ", dw[i]);
                flint_printf("\n");
                for (i = 0; i < n + ball; i++) flint_printf("%.17g ", ds[i]);
                flint_printf("\n");
                flint_abort();
            }
        }
        arf_clear(a); arf_clear(b);

        gr_heap_clear_vec(x, len + 1, wctx);
        gr_heap_clear_vec(y, len + 1, wctx);
        gr_heap_clear_vec(rw, len + 1, wctx);
        gr_heap_clear_vec(rs, len + 1, wctx);
        gr_ctx_clear(wctx);
        gr_ctx_clear(sctx);
    }

    TEST_FUNCTION_END(state);
}
