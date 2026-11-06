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
#include "acb.h"
#include "gr.h"
#include "gr_vec.h"
#include "dfloat.h"

/* The plain complex vector kernels used by the eigensolvers, on data
   with the structure that sends the SIMD code into its special cases
   (exact zeros, real and imaginary values, exact small integers whose
   combinations cancel, nonfinite values):

   - _gr_vec_mul_scalar, _gr_vec_addmul_scalar, _gr_vec_submul_scalar
     agree bitwise with the elementwise gr_mul and gr_add / gr_sub;
   - _gr_vec_apply_2x2 gives every pair the same result whatever its
     position (a block of four, a partial block, a single pair: the
     SIMD and the scalar code), is accurate against acb for finite
     data, and with the strong overlay gives the canonical form of the
     plain result;
   - gr_abs2, gr_cmp and gr_cmpabs agree with acb. */

/* d[0..n-1] the real part, d[n..2n-1] the imaginary part */
static void
_t_cv_part(double * d, int n, flint_rand_t state, int kind)
{
    int k;
    for (k = 0; k < n; k++)
        d[k] = 0.0;
    switch (kind)
    {
        case 0:
            break;
        case 1:
            d[0] = (double) ((slong) n_randint(state, 7) - 3);
            break;
        case 2:
            d[0] = (n_randint(state, 2) ? INFINITY : -INFINITY);
            break;
        case 3:
            d[0] = NAN;
            break;
        default:
        {
            /* a normalized expansion with random components */
            double s = ldexp(1.0, (int) n_randint(state, 40) - 20);
            d[0] = s * (1.0 + (double) n_randint(state, 1 << 30) * 0x1p-30) * (n_randint(state, 2) ? 1 : -1);
            for (k = 1; k < n; k++)
                d[k] = ldexp(d[k - 1], -53 - (int) n_randint(state, 3)) * (n_randint(state, 2) ? 0.75 : -0.625);
            break;
        }
    }
}

static void
_t_cv_rand(gr_ptr x, int n, flint_rand_t state, int allow_special)
{
    double * d = x;
    int kre, kim;
    ulong r = n_randint(state, 16);
    /* mostly generic, otherwise one of the special structures */
    kre = (r < 8) ? 4 : (int) n_randint(state, 5);
    kim = (r < 8) ? 4 : (int) n_randint(state, 5);
    if (!allow_special || n_randint(state, 8) != 0)
    {
        if (kre == 2 || kre == 3) kre = 4;
        if (kim == 2 || kim == 3) kim = 1;
    }
    _t_cv_part(d, n, state, kre);
    _t_cv_part(d + n, n, state, kim);
}

static int
_t_cv_equal(gr_srcptr a, gr_srcptr b, int n)
{
    const double * x = a, * y = b;
    int k;
    for (k = 0; k < 2 * n; k++)
        if (x[k] != y[k] && !(x[k] != x[k] && y[k] != y[k]))
            return 0;
        else if (x[k] == 0.0 && y[k] == 0.0 && signbit(x[k]) != signbit(y[k]))
            return 0;
    return 1;
}

static int
_t_cv_finite(gr_srcptr a, int n)
{
    const double * x = a;
    int k;
    for (k = 0; k < 2 * n; k++)
        if (!(fabs(x[k]) <= 1e300))
            return 0;
    return 1;
}

TEST_FUNCTION_START(complex_vec, state)
{
    slong iter;

    DFLOAT_TEST_SKIP_IF_UNSUPPORTED(state);

    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t ctx, sctx, CC;
        int n, which;
        slong len, i, stride, sz;
        gr_ptr x, y, r1, r2, c, t, M;
        acb_t a, b, e, m;

        n = 1 + n_randint(state, DFLOAT_MAX_N);
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(ctx, n, DFLOAT_COMPLEX));
        GR_MUST_SUCCEED(gr_ctx_init_dfloat(sctx, n, DFLOAT_COMPLEX | DFLOAT_STRONG));
        gr_ctx_init_complex_acb(CC, 53 * n + 200);
        sz = ctx->sizeof_elem;
        which = n_randint(state, 3);
        len = n_randint(state, 4) ? n_randint(state, 14) : n_randint(state, 40);
        stride = 1 + n_randint(state, 3);

        x = gr_heap_init_vec(len * stride + 1, ctx);
        y = gr_heap_init_vec(len * stride + 1, ctx);
        r1 = gr_heap_init_vec(len * stride + 1, ctx);
        r2 = gr_heap_init_vec(len * stride + 1, ctx);
        M = gr_heap_init_vec(4, ctx);
        c = gr_heap_init(ctx);
        t = gr_heap_init(ctx);
        acb_init(a); acb_init(b); acb_init(e); acb_init(m);

        for (i = 0; i < len * stride; i++)
        {
            _t_cv_rand(GR_ENTRY(x, i, sz), n, state, 1);
            _t_cv_rand(GR_ENTRY(y, i, sz), n, state, 1);
        }

        if (which == 0)
        {
            /* products by a scalar: bitwise the elementwise operations */
            int op = n_randint(state, 3), alias = n_randint(state, 2);
            _t_cv_rand(c, n, state, 0);
            GR_MUST_SUCCEED(_gr_vec_set(r1, y, len, ctx));
            GR_MUST_SUCCEED(_gr_vec_set(r2, y, len, ctx));
            if (alias)
                GR_MUST_SUCCEED(_gr_vec_set(x, y, len, ctx));
            for (i = 0; i < len; i++)
            {
                gr_ptr ri = GR_ENTRY(r2, i, sz);
                GR_MUST_SUCCEED(gr_mul(t, GR_ENTRY(x, i, sz), c, ctx));
                if (op == 0) GR_MUST_SUCCEED(gr_set(ri, t, ctx));
                else if (op == 1) GR_MUST_SUCCEED(gr_add(ri, ri, t, ctx));
                else GR_MUST_SUCCEED(gr_sub(ri, ri, t, ctx));
            }
            if (op == 0) GR_MUST_SUCCEED(_gr_vec_mul_scalar(r1, alias ? r1 : x, len, c, ctx));
            else if (op == 1) GR_MUST_SUCCEED(_gr_vec_addmul_scalar(r1, alias ? r1 : x, len, c, ctx));
            else GR_MUST_SUCCEED(_gr_vec_submul_scalar(r1, alias ? r1 : x, len, c, ctx));
            for (i = 0; i < len; i++)
            {
                if (!_t_cv_equal(GR_ENTRY(r1, i, sz), GR_ENTRY(r2, i, sz), n))
                {
                    flint_printf("FAIL: scalar product (n = %d, op = %d, alias = %d, len = %wd, i = %wd)\n", n, op, alias, len, i);
                    gr_println(GR_ENTRY(x, i, sz), ctx); gr_println(c, ctx);
                    gr_println(GR_ENTRY(r1, i, sz), ctx); gr_println(GR_ENTRY(r2, i, sz), ctx);
                    flint_abort();
                }
            }
        }
        else if (which == 1)
        {
            /* 2 x 2 transformations */
            int kind = n_randint(state, 4), k;

            for (k = 0; k < 4; k++)
                _t_cv_rand(GR_ENTRY(M, k, sz), n, state, n_randint(state, 4) == 0);
            if (kind == 1 || kind == 2)
            {
                /* real diagonal, or real; as in the plane rotations */
                for (k = 0; k < 4; k++)
                    if (kind == 2 || k == 0 || k == 3)
                        _t_cv_part((double *) GR_ENTRY(M, k, sz) + n, n, state, 0);
            }

            GR_MUST_SUCCEED(_gr_vec_set(r1, x, len * stride, ctx));
            GR_MUST_SUCCEED(_gr_vec_set(r2, y, len * stride, ctx));
            GR_MUST_SUCCEED(_gr_vec_apply_2x2(r1, r2, stride, len, M, ctx));

            for (i = 0; i < len; i++)
            {
                gr_ptr u = GR_ENTRY(x, i * stride, sz), v = GR_ENTRY(y, i * stride, sz);
                gr_ptr ru = GR_ENTRY(r1, i * stride, sz), rv = GR_ENTRY(r2, i * stride, sz);
                gr_ptr p, q;
                int ok;

                /* the pair alone, and in the last place of a partial
                   block preceded by copies of other pairs */
                p = gr_heap_init_vec(4, ctx);
                q = gr_heap_init_vec(4, ctx);
                for (k = 0; k < 2; k++)
                {
                    slong l = k ? 1 + n_randint(state, 3) : 1;
                    slong s;
                    for (s = 0; s < l - 1; s++)
                    {
                        GR_MUST_SUCCEED(gr_set(GR_ENTRY(p, s, sz), GR_ENTRY(x, n_randint(state, len * stride), sz), ctx));
                        GR_MUST_SUCCEED(gr_set(GR_ENTRY(q, s, sz), GR_ENTRY(y, n_randint(state, len * stride), sz), ctx));
                    }
                    GR_MUST_SUCCEED(gr_set(GR_ENTRY(p, l - 1, sz), u, ctx));
                    GR_MUST_SUCCEED(gr_set(GR_ENTRY(q, l - 1, sz), v, ctx));
                    GR_MUST_SUCCEED(_gr_vec_apply_2x2(p, q, 1, l, M, ctx));
                    if (!_t_cv_equal(GR_ENTRY(p, l - 1, sz), ru, n) || !_t_cv_equal(GR_ENTRY(q, l - 1, sz), rv, n))
                    {
                        flint_printf("FAIL: apply_2x2 position (n = %d, kind = %d, len = %wd, stride = %wd, i = %wd, l = %wd)\n", n, kind, len, stride, i, l);
                        for (s = 0; s < 4; s++) gr_println(GR_ENTRY(M, s, sz), ctx);
                        gr_println(u, ctx); gr_println(v, ctx);
                        gr_println(ru, ctx); gr_println(GR_ENTRY(p, l - 1, sz), ctx);
                        gr_println(rv, ctx); gr_println(GR_ENTRY(q, l - 1, sz), ctx);
                        flint_abort();
                    }
                }

                /* the strong overlay: the same values */
                GR_MUST_SUCCEED(gr_set(p, u, ctx));
                GR_MUST_SUCCEED(gr_set(q, v, ctx));
                GR_MUST_SUCCEED(_gr_vec_apply_2x2(p, q, 1, 1, M, sctx));
                if (_t_cv_finite(ru, n) && _t_cv_finite(rv, n)
                    && (gr_equal(p, ru, ctx) != T_TRUE || gr_equal(q, rv, ctx) != T_TRUE))
                {
                    flint_printf("FAIL: apply_2x2 strong (n = %d, i = %wd)\n", n, i);
                    gr_println(ru, ctx); gr_println(p, ctx);
                    gr_println(rv, ctx); gr_println(q, ctx);
                    flint_abort();
                }

                gr_heap_clear_vec(p, 4, ctx);
                gr_heap_clear_vec(q, 4, ctx);

                /* accuracy: |r - exact| <= 2^(-53 n + 8) (|Mj0| |u| + |Mj1| |v|) */
                ok = _t_cv_finite(u, n) && _t_cv_finite(v, n) && _t_cv_finite(ru, n) && _t_cv_finite(rv, n);
                for (k = 0; k < 4; k++)
                    ok = ok && _t_cv_finite(GR_ENTRY(M, k, sz), n);
                if (ok)
                {
                    acb_t mu, mv, mm;
                    arb_t err, bnd, tt;
                    int row;
                    acb_init(mu); acb_init(mv); acb_init(mm);
                    arb_init(err); arb_init(bnd); arb_init(tt);
                    GR_MUST_SUCCEED(gr_set_other(mu, u, ctx, CC));
                    GR_MUST_SUCCEED(gr_set_other(mv, v, ctx, CC));
                    for (row = 0; row < 2; row++)
                    {
                        GR_MUST_SUCCEED(gr_set_other(mm, GR_ENTRY(M, 2 * row, sz), ctx, CC));
                        acb_mul(a, mm, mu, 53 * n + 200);
                        acb_abs(bnd, mm, 53 * n + 200);
                        acb_abs(tt, mu, 53 * n + 200);
                        arb_mul(bnd, bnd, tt, 53 * n + 200);
                        GR_MUST_SUCCEED(gr_set_other(mm, GR_ENTRY(M, 2 * row + 1, sz), ctx, CC));
                        acb_mul(b, mm, mv, 53 * n + 200);
                        acb_add(a, a, b, 53 * n + 200);
                        acb_abs(tt, mm, 53 * n + 200);
                        acb_abs(err, mv, 53 * n + 200);
                        arb_mul(tt, tt, err, 53 * n + 200);
                        arb_add(bnd, bnd, tt, 53 * n + 200);
                        arb_mul_2exp_si(bnd, bnd, -53 * n + 8);
                        GR_MUST_SUCCEED(gr_set_other(b, row ? rv : ru, ctx, CC));
                        acb_sub(b, b, a, 53 * n + 200);
                        acb_abs(err, b, 53 * n + 200);
                        /* (an absolute allowance for underflow) */
                        arb_set_d(tt, 1e-280);
                        arb_add(bnd, bnd, tt, 53 * n + 200);
                        if (!arb_le(err, bnd))
                        {
                            flint_printf("FAIL: apply_2x2 accuracy (n = %d, kind = %d, i = %wd, row = %d)\n", n, kind, i, row);
                            for (k = 0; k < 4; k++) gr_println(GR_ENTRY(M, k, sz), ctx);
                            gr_println(u, ctx); gr_println(v, ctx);
                            gr_println(row ? rv : ru, ctx);
                            arb_printd(err, 10); flint_printf("\n"); arb_printd(bnd, 10); flint_printf("\n");
                            flint_abort();
                        }
                    }
                    acb_clear(mu); acb_clear(mv); acb_clear(mm);
                    arb_clear(err); arb_clear(bnd); arb_clear(tt);
                }
            }
        }
        else
        {
            /* abs2, cmp, cmpabs against acb */
            for (i = 0; i < FLINT_MIN(len, 6); i++)
            {
                gr_ptr u = GR_ENTRY(x, i, sz), v = GR_ENTRY(y, i, sz);
                int st, cmp, ref;
                arb_t s1, s2;

                if (!_t_cv_finite(u, n) || !_t_cv_finite(v, n))
                    continue;
                if (n_randint(state, 4) == 0)
                    GR_MUST_SUCCEED(gr_set(v, u, ctx));
                if (n_randint(state, 4) == 0)
                    GR_MUST_SUCCEED(gr_neg(v, u, ctx));
                if (n_randint(state, 4) == 0)
                    GR_MUST_SUCCEED(gr_conj(v, u, ctx));

                arb_init(s1); arb_init(s2);
                GR_MUST_SUCCEED(gr_set_other(a, u, ctx, CC));
                GR_MUST_SUCCEED(gr_set_other(b, v, ctx, CC));

                /* abs2: real, within 2^(-53n+4) relative */
                GR_MUST_SUCCEED(gr_abs2(t, u, ctx));
                GR_MUST_SUCCEED(gr_set_other(e, t, ctx, CC));
                acb_abs(s1, a, 53 * n + 200);
                arb_sqr(s1, s1, 53 * n + 200);
                arb_sub(s2, acb_realref(e), s1, 53 * n + 200);
                arb_abs(s2, s2);
                arb_mul_2exp_si(s1, s1, -53 * n + 4);
                if (!arb_is_zero(acb_imagref(e)) || !arb_le(s2, s1))
                {
                    flint_printf("FAIL: abs2 (n = %d)\n", n);
                    gr_println(u, ctx); gr_println(t, ctx);
                    flint_abort();
                }

                /* cmp: only for real values */
                st = gr_cmp(&cmp, u, v, ctx);
                if (arb_is_zero(acb_imagref(a)) && arb_is_zero(acb_imagref(b)))
                {
                    ref = arf_cmp(arb_midref(acb_realref(a)), arb_midref(acb_realref(b)));
                    if (st != GR_SUCCESS || cmp != ref)
                    {
                        flint_printf("FAIL: cmp real (n = %d, %d, %d, %d)\n", n, st, cmp, ref);
                        gr_println(u, ctx); gr_println(v, ctx);
                        flint_abort();
                    }
                }
                else if (st != GR_DOMAIN)
                {
                    flint_printf("FAIL: cmp nonreal (n = %d, %d)\n", n, st);
                    gr_println(u, ctx); gr_println(v, ctx);
                    flint_abort();
                }

                /* cmpabs: the sign of |u| - |v| (exactly decided for
                   |u|^2 - |v|^2 as an exact acb computation) */
                st = gr_cmpabs(&cmp, u, v, ctx);
                acb_get_mid(e, a);
                acb_get_mid(m, b);
                arb_sqr(s1, acb_realref(e), ARF_PREC_EXACT);
                arb_addmul(s1, acb_imagref(e), acb_imagref(e), ARF_PREC_EXACT);
                arb_sqr(s2, acb_realref(m), ARF_PREC_EXACT);
                arb_addmul(s2, acb_imagref(m), acb_imagref(m), ARF_PREC_EXACT);
                ref = arf_cmp(arb_midref(s1), arb_midref(s2));
                if (st != GR_SUCCESS || cmp != ref)
                {
                    flint_printf("FAIL: cmpabs (n = %d, %d, %d, %d)\n", n, st, cmp, ref);
                    gr_println(u, ctx); gr_println(v, ctx);
                    flint_abort();
                }
                arb_clear(s1); arb_clear(s2);
            }
        }

        gr_heap_clear_vec(x, len * stride + 1, ctx);
        gr_heap_clear_vec(y, len * stride + 1, ctx);
        gr_heap_clear_vec(r1, len * stride + 1, ctx);
        gr_heap_clear_vec(r2, len * stride + 1, ctx);
        gr_heap_clear_vec(M, 4, ctx);
        gr_heap_clear(c, ctx);
        gr_heap_clear(t, ctx);
        acb_clear(a); acb_clear(b); acb_clear(e); acb_clear(m);
        gr_ctx_clear(ctx);
        gr_ctx_clear(sctx);
        gr_ctx_clear(CC);
    }

    TEST_FUNCTION_END(state);
}
