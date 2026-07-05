/*
    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* generated using Claude Opus 4.8 */

#include "test_helpers.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_ore_poly.h"

static const ore_algebra_t ds_diff_algs[2] = {
    ORE_ALGEBRA_DERIVATIVE, ORE_ALGEBRA_EULER_DERIVATIVE
};

static const ore_algebra_t ds_shift_algs[4] = {
    ORE_ALGEBRA_FORWARD_SHIFT, ORE_ALGEBRA_BACKWARD_SHIFT,
    ORE_ALGEBRA_FORWARD_DIFFERENCE, ORE_ALGEBRA_BACKWARD_DIFFERENCE
};

/* Returns 0 only if a == x^net * b is provably false */
static int
check_scaled(const gr_ore_poly_t a, const gr_ore_poly_t b, slong net,
             gr_ore_poly_ctx_t octx)
{
    int ok = 1, status = GR_SUCCESS;
    gr_ctx_struct * cctx = GR_ORE_POLY_ELEM_CTX(octx);
    slong sz = cctx->sizeof_elem, i;
    const gr_ore_poly_struct * lo = (net >= 0) ? b : a;
    const gr_ore_poly_struct * hi = (net >= 0) ? a : b;
    gr_ptr xp;
    gr_ore_poly_t scaled;

    GR_TMP_INIT(xp, cctx);
    gr_ore_poly_init(scaled, octx);

    status |= gr_gen(xp, cctx);
    status |= gr_pow_ui(xp, xp, (ulong) FLINT_ABS(net), cctx);
    gr_ore_poly_fit_length(scaled, lo->length, octx);
    for (i = 0; i < lo->length; i++)
        status |= gr_mul(GR_ENTRY(scaled->coeffs, i, sz), xp, GR_ENTRY(lo->coeffs, i, sz), cctx);
    _gr_ore_poly_set_length(scaled, lo->length, octx);
    _gr_ore_poly_normalise(scaled, octx);

    if (status == GR_SUCCESS && gr_ore_poly_equal(scaled, hi, octx) == T_FALSE)
        ok = 0;

    gr_ore_poly_clear(scaled, octx);
    GR_TMP_CLEAR(xp, cctx);
    return ok;
}

/* Returns 0 only if a == S^net * b is provably false, where S is the forward
   shift operator. */
static int
check_shifted(const gr_ore_poly_t a, const gr_ore_poly_t b, slong net,
              gr_ore_poly_ctx_t octx)
{
    int ok = 1, status = GR_SUCCESS;
    slong var = GR_ORE_POLY_ORE_DATA(octx)->base_var;
    gr_ctx_struct * base = GR_ORE_POLY_ELEM_CTX(octx);
    gr_ore_poly_ctx_t fs_ctx;
    gr_ore_poly_t af, bf, s, lhs, rhs;
    slong ea, eb, la, lb, c;

    gr_ore_poly_ctx_init(fs_ctx, base, var, ORE_ALGEBRA_FORWARD_SHIFT);
    gr_ore_poly_init(af, fs_ctx);
    gr_ore_poly_init(bf, fs_ctx);
    gr_ore_poly_init(s, fs_ctx);
    gr_ore_poly_init(lhs, fs_ctx);
    gr_ore_poly_init(rhs, fs_ctx);

    status |= gr_ore_poly_convert(af, &ea, a, fs_ctx, octx);
    status |= gr_ore_poly_convert(bf, &eb, b, fs_ctx, octx);

    la = ea;
    lb = net + eb;
    c = FLINT_MAX(0, -FLINT_MIN(la, lb));

    status |= gr_ore_poly_gen(s, fs_ctx);
    status |= gr_pow_ui(lhs, s, (ulong) (la + c), fs_ctx);
    status |= gr_ore_poly_mul(lhs, lhs, af, fs_ctx);
    status |= gr_pow_ui(rhs, s, (ulong) (lb + c), fs_ctx);
    status |= gr_ore_poly_mul(rhs, rhs, bf, fs_ctx);

    if (status == GR_SUCCESS && gr_ore_poly_equal(lhs, rhs, fs_ctx) == T_FALSE)
        ok = 0;

    gr_ore_poly_clear(af, fs_ctx);
    gr_ore_poly_clear(bf, fs_ctx);
    gr_ore_poly_clear(s, fs_ctx);
    gr_ore_poly_clear(lhs, fs_ctx);
    gr_ore_poly_clear(rhs, fs_ctx);
    gr_ore_poly_ctx_clear(fs_ctx);
    return ok;
}

/* Applies a recurrence R = sum_k q_k(nu) S^k to the coefficients of a
   polynomial. */
static int
apply_forward_recurrence(gr_poly_t s, const gr_ore_poly_t R, const gr_poly_t f, gr_ctx_t cctx)
{
    gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);
    slong bsz = cctx->sizeof_elem, ssz = sctx->sizeof_elem;
    slong df = f->length;
    int status = GR_SUCCESS;
    gr_ptr mval, qval, term;

    GR_TMP_INIT3(mval, qval, term, sctx);
    gr_poly_fit_length(s, df, sctx);
    for (slong n = 0; n < df; n++)
    {
        gr_ptr sm = GR_ENTRY(s->coeffs, n, ssz);
        status |= gr_zero(sm, sctx);
        status |= gr_set_si(mval, n, sctx);
        for (slong k = 0; k < R->length && n + k < df; k++)
        {
            const gr_poly_struct * qk = (const gr_poly_struct *) GR_ENTRY(R->coeffs, k, bsz);
            status |= gr_poly_evaluate(qval, qk, mval, sctx);
            status |= gr_mul(term, qval, GR_ENTRY(f->coeffs, n + k, ssz), sctx);
            status |= gr_add(sm, sm, term, sctx);
        }
    }
    _gr_poly_set_length(s, df, sctx);
    _gr_poly_normalise(s, sctx);
    GR_TMP_CLEAR3(mval, qval, term, sctx);
    return status;
}

TEST_GR_FUNCTION_START(gr_ore_poly_differential_to_shift, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_ore_poly_ctx_t ctx_d, ctx_s;
        ore_algebra_t diff_alg, shift_alg;
        slong p1, p2, ngens, var;
        int status = GR_SUCCESS;
        int expect_success = 1;

        switch (n_randint(state, 8))
        {
            case 0:
                gr_ctx_init_random(cctx, state);
                break;
            case 1:
                gr_ctx_init_random_mpoly(cctx, state);
                break;
            default:
                gr_ctx_init_random_poly(cctx, state);
                break;
        }

        ngens = 0;
        status |= gr_ctx_ngens(&ngens, cctx);
        var = n_randint(state, ngens);

        diff_alg = ds_diff_algs[n_randint(state, 2)];
        shift_alg = ds_shift_algs[n_randint(state, 4)];

        if (!n_randint(state, 16))
        {
            diff_alg = ORE_ALGEBRA_COMMUTATIVE;
            expect_success = 0;
        }
        if (!n_randint(state, 16))
        {
            shift_alg = ORE_ALGEBRA_COMMUTATIVE;
            expect_success = 0;
        }

        gr_ore_poly_ctx_init(ctx_d, cctx, var, diff_alg);
        gr_ore_poly_ctx_init(ctx_s, cctx, var, shift_alg);

        gr_ore_poly_t op, res_s, op2;
        gr_ore_poly_init(op, ctx_d);
        gr_ore_poly_init(res_s, ctx_s);
        gr_ore_poly_init(op2, ctx_d);

        status |= gr_ore_poly_randtest(op, state, 1 + n_randint(state, 5), ctx_d);

        status |= gr_ore_poly_differential_to_shift(res_s, &p1, op, ctx_s, ctx_d);
        status |= gr_ore_poly_shift_to_differential(op2, &p2, res_s, ctx_d, ctx_s);

        expect_success = (expect_success && var == 0
            && cctx->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_FMPZ
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_CC_ACB
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_NMOD));
        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        /* round trip */

        if (status == GR_SUCCESS && !check_scaled(op2, op, p1 - p2, ctx_d))
        {
            flint_printf("FAIL: round trip\n");
            flint_abort();
        }

        /* a forward-shift recurrence R from op satisfies
           (R a)_m = coeff_{m-pf}(op . f) for the coefficient sequence of f */

        if (cctx->which_ring == GR_CTX_GR_POLY && shift_alg == ORE_ALGEBRA_FORWARD_SHIFT)
        {
            gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);
            gr_ptr f = gr_heap_init(cctx);
            gr_ptr g = gr_heap_init(cctx);
            gr_poly_t s, eg;

            gr_poly_init(s, sctx);
            gr_poly_init(eg, sctx);

            status |= gr_randtest(f, state, cctx);
            status |= gr_ore_poly_apply(g, op, f, ctx_d);
            status |= apply_forward_recurrence(s, res_s, (gr_poly_struct *) f, cctx);
            if (p1 >= 0)
                status |= gr_poly_shift_left(eg, (gr_poly_struct *) g, p1, sctx);
            else
                status |= gr_poly_shift_right(eg, (gr_poly_struct *) g, -p1, sctx);

            if (status == GR_SUCCESS && gr_poly_equal(s, eg, sctx) == T_FALSE)
            {
                flint_printf("FAIL: action\n");
                flint_abort();
            }

            gr_poly_clear(s, sctx);
            gr_poly_clear(eg, sctx);
            gr_heap_clear(f, cctx);
            gr_heap_clear(g, cctx);
        }

        gr_ore_poly_clear(op, ctx_d);
        gr_ore_poly_clear(res_s, ctx_s);
        gr_ore_poly_clear(op2, ctx_d);
        gr_ore_poly_ctx_clear(ctx_d);
        gr_ore_poly_ctx_clear(ctx_s);
        gr_ctx_clear(cctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}

TEST_GR_FUNCTION_START(gr_ore_poly_shift_to_differential, state, count_success, count_domain, count_unable)
{
    for (slong iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        gr_ctx_t cctx;
        gr_ore_poly_ctx_t ctx_d, ctx_s;
        ore_algebra_t diff_alg, shift_alg;
        slong p1, p2, ngens, var;
        int status = GR_SUCCESS;
        int expect_success = 1;

        switch (n_randint(state, 8))
        {
            case 0:
                gr_ctx_init_random(cctx, state);
                break;
            case 1:
                gr_ctx_init_random_mpoly(cctx, state);
                break;
            default:
                gr_ctx_init_random_poly(cctx, state);
                break;
        }

        ngens = 0;
        status |= gr_ctx_ngens(&ngens, cctx);
        var = n_randint(state, ngens);

        diff_alg = ds_diff_algs[n_randint(state, 2)];
        shift_alg = ds_shift_algs[n_randint(state, 4)];


        if (!n_randint(state, 16))
        {
            diff_alg = ORE_ALGEBRA_COMMUTATIVE;
            expect_success = 0;
        }
        if (!n_randint(state, 16))
        {
            shift_alg = ORE_ALGEBRA_COMMUTATIVE;
            expect_success = 0;
        }

        gr_ore_poly_ctx_init(ctx_d, cctx, var, diff_alg);
        gr_ore_poly_ctx_init(ctx_s, cctx, var, shift_alg);

        gr_ore_poly_t op, res_d, op2;
        gr_ore_poly_init(op, ctx_s);
        gr_ore_poly_init(res_d, ctx_d);
        gr_ore_poly_init(op2, ctx_s);

        status |= gr_ore_poly_randtest(op, state, 1 + n_randint(state, 5), ctx_s);

        status |= gr_ore_poly_shift_to_differential(res_d, &p1, op, ctx_d, ctx_s);
        status |= gr_ore_poly_differential_to_shift(op2, &p2, res_d, ctx_s, ctx_d);

        expect_success = (expect_success && var == 0
            && cctx->which_ring == GR_CTX_GR_POLY
            && (POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_FMPZ
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_CC_ACB
                || POLYNOMIAL_ELEM_CTX(cctx)->which_ring == GR_CTX_NMOD));
        if (expect_success && status != GR_SUCCESS)
        {
            flint_printf("FAIL: unexpected failure\n");
            flint_abort();
        }

        count_success += (status == GR_SUCCESS);
        count_domain += ((status & GR_DOMAIN) != 0);
        count_unable += ((status & GR_UNABLE) != 0);

        /* round trip recovers S^(p1-p2) * op */

        if (status == GR_SUCCESS && !check_shifted(op2, op, p1 - p2, ctx_s))
        {
            flint_printf("FAIL: round trip\n");
            flint_abort();
        }

        if (cctx->which_ring == GR_CTX_GR_POLY && shift_alg == ORE_ALGEBRA_FORWARD_SHIFT)
        {
            gr_ctx_struct * sctx = POLYNOMIAL_ELEM_CTX(cctx);
            gr_ptr g = gr_heap_init(cctx);
            gr_poly_t s, eg;

            gr_poly_init(s, sctx);
            gr_poly_init(eg, sctx);

            slong flen = op->length + FLINT_ABS(p1) + 2 + n_randint(state, 5);
            gr_ptr f = gr_heap_init(cctx);
            status |= gr_poly_randtest((gr_poly_struct *) f, state, flen, sctx);

            status |= gr_ore_poly_apply(g, res_d, f, ctx_d);
            status |= apply_forward_recurrence(s, op, (gr_poly_struct *) f, cctx);
            if (p1 >= 0)
                status |= gr_poly_shift_left(eg, (gr_poly_struct *) g, p1, sctx);
            else
                status |= gr_poly_shift_right(eg, (gr_poly_struct *) g, -p1, sctx);

            if (status == GR_SUCCESS && gr_poly_equal(s, eg, sctx) == T_FALSE)
            {
                flint_printf("FAIL: action\n");
                flint_abort();
            }

            gr_poly_clear(s, sctx);
            gr_poly_clear(eg, sctx);
            gr_heap_clear(f, cctx);
            gr_heap_clear(g, cctx);
        }

        gr_ore_poly_clear(op, ctx_s);
        gr_ore_poly_clear(res_d, ctx_d);
        gr_ore_poly_clear(op2, ctx_s);
        gr_ore_poly_ctx_clear(ctx_d);
        gr_ore_poly_ctx_clear(ctx_s);
        gr_ctx_clear(cctx);
    }

    TEST_GR_FUNCTION_END(state, count_success, count_domain, count_unable);
}
