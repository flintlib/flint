/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"
#include "fmpz_mpoly_factor.h"
#include "fmpz_mpoly_q.h"
#include "gr.h"
#include "gr/impl.h"
#include "gr_vec.h"
#include "gr_generic.h"

typedef struct
{
    fmpz_mpoly_ctx_t mctx;
    char ** vars;
}
_gr_fmpz_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpz_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

static int _gr_fmpz_mpoly_q_factored_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Fraction field of multivariate polynomials over Integer ring (fmpz)");
    status |= gr_stream_write(out, " in ");
    status |= gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->minfo->nvars);
    status |= gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_LEX)
        status |= gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGLEX)
        status |= gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGREVLEX)
        status |= gr_stream_write(out, ", degrevlex order");
    status |= gr_stream_write(out, ", factored denominators");
    return status;
}

/* Some methods are identical to their fmpz_mpoly counterparts */
#define _gr_fmpz_mpoly_q_factored_ctx_clear _gr_fmpz_mpoly_ctx_clear
#define _gr_fmpz_mpoly_q_factored_ctx_set_gen_names _gr_fmpz_mpoly_ctx_set_gen_names
#define _gr_fmpz_mpoly_q_factored_ctx_ngens _gr_fmpz_mpoly_ctx_ngens
#define _gr_fmpz_mpoly_q_factored_ctx_gen_name _gr_fmpz_mpoly_ctx_gen_name

static void
_gr_fmpz_mpoly_q_factored_init(fmpz_mpoly_q_factored_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_init(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mpoly_q_factored_clear(fmpz_mpoly_q_factored_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mpoly_q_factored_swap(fmpz_mpoly_q_factored_t x, fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_swap(x, y, MPOLYNOMIAL_MCTX(ctx));
}

static void
_gr_fmpz_mpoly_q_factored_set_shallow(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    *res = *x;
}

static int
_gr_fmpz_mpoly_q_factored_randtest(fmpz_mpoly_q_factored_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong bits, length;
    flint_bitcnt_t exp_bound;

    if (n_randint(state, 10) != 0)
        bits = 10;
    else
        bits = 100;

    /* Split state-mutating calls so consumption order is architecture-independent. */
    length = n_randint(state, 5);
    exp_bound = 1 + n_randint(state, 3);
    fmpz_mpoly_q_factored_randtest(res, state, length, bits, exp_bound, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_randtest_small(fmpz_mpoly_q_factored_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong length;
    flint_bitcnt_t exp_bound;

    length = n_randint(state, 3);
    exp_bound = 1 + n_randint(state, 3);
    fmpz_mpoly_q_factored_randtest(res, state, length, 3, exp_bound, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static slong
_gr_fmpz_mpoly_q_factored_length(const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    slong i, len = x->num.length;

    for (i = 0; i < x->den.num; i++)
        len += x->den.poly[i].length;

    return len;
}

static int
_gr_fmpz_mpoly_q_factored_write(gr_stream_t out, fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    return gr_stream_write_free(out, fmpz_mpoly_q_factored_get_str_pretty(x,
        (const char **) MPOLYNOMIAL_CTX(ctx)->vars, MPOLYNOMIAL_MCTX(ctx)));
}

static truth_t
_gr_fmpz_mpoly_q_factored_equal(const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_factored_equal(x, y, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpz_mpoly_q_factored_is_zero(const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_factored_is_zero(x, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static truth_t
_gr_fmpz_mpoly_q_factored_is_one(const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_factored_is_one(x, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

static int
_gr_fmpz_mpoly_q_factored_zero(fmpz_mpoly_q_factored_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_one(fmpz_mpoly_q_factored_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;

    n = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_factored_gen(((fmpz_mpoly_q_factored_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx));

    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set(res, x, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set_si(fmpz_mpoly_q_factored_t res, slong v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set_si(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set_ui(fmpz_mpoly_q_factored_t res, ulong v, gr_ctx_t ctx)
{
    fmpz_mpoly_set_ui(fmpz_mpoly_q_factored_numref(res), v, MPOLYNOMIAL_MCTX(ctx));
    fmpz_mpoly_factor_one(fmpz_mpoly_q_factored_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_t v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set_fmpq(fmpz_mpoly_q_factored_t res, const fmpq_t v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_set_other(fmpz_mpoly_q_factored_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx->which_ring == GR_CTX_FMPZ_MPOLY_Q_FACTORED ||
        x_ctx->which_ring == GR_CTX_FMPZ_MPOLY_Q ||
        x_ctx->which_ring == GR_CTX_FMPZ_MPOLY)
    {
        if (MPOLYNOMIAL_MCTX(ctx)->minfo->nvars == MPOLYNOMIAL_MCTX(x_ctx)->minfo->nvars &&
            MPOLYNOMIAL_MCTX(ctx)->minfo->ord == MPOLYNOMIAL_MCTX(x_ctx)->minfo->ord)
        {
            if (x_ctx->which_ring == GR_CTX_FMPZ_MPOLY_Q_FACTORED)
                fmpz_mpoly_q_factored_set(res, x, MPOLYNOMIAL_MCTX(ctx));
            else if (x_ctx->which_ring == GR_CTX_FMPZ_MPOLY_Q)
                fmpz_mpoly_q_factored_set_fmpz_mpoly_q(res, x, MPOLYNOMIAL_MCTX(ctx));
            else
                fmpz_mpoly_q_factored_set_fmpz_mpoly(res, x, MPOLYNOMIAL_MCTX(ctx));

            return GR_SUCCESS;
        }
    }

    return gr_generic_set_other(res, x, x_ctx, ctx);
}

static int
_gr_fmpz_mpoly_q_factored_neg(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_neg(res, x, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_add(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_add(res, x, y, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_add_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_add_si(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_add_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_add_fmpz(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_add_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_add_fmpq(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_sub(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_sub(res, x, y, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_sub_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_sub_si(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_sub_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_sub_fmpz(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_sub_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_sub_fmpq(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_mul(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_mul(res, x, y, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_mul_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_mul_si(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_mul_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_mul_fmpz(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_mul_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_mul_fmpq(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_div(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_mpoly_q_factored_t y, gr_ctx_t ctx)
{
    if (fmpz_mpoly_q_factored_is_zero(y, MPOLYNOMIAL_MCTX(ctx)))
        return GR_DOMAIN;

    fmpz_mpoly_q_factored_div(res, x, y, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_div_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpz_mpoly_q_factored_div_si(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_div_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_is_zero(c))
        return GR_DOMAIN;

    fmpz_mpoly_q_factored_div_fmpz(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_div_fmpq(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpq_t c, gr_ctx_t ctx)
{
    if (fmpq_is_zero(c))
        return GR_DOMAIN;

    fmpz_mpoly_q_factored_div_fmpq(res, x, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static truth_t
_gr_fmpz_mpoly_q_factored_is_invertible(const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_factored_is_zero(x, MPOLYNOMIAL_MCTX(ctx)) ? T_FALSE : T_TRUE;
}

static int
_gr_fmpz_mpoly_q_factored_inv(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, gr_ctx_t ctx)
{
    if (fmpz_mpoly_q_factored_is_zero(x, MPOLYNOMIAL_MCTX(ctx)))
        return GR_DOMAIN;

    fmpz_mpoly_q_factored_inv(res, x, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_pow_fmpz(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const fmpz_t e, gr_ctx_t ctx)
{
    const fmpz_mpoly_ctx_struct * mctx = MPOLYNOMIAL_MCTX(ctx);

    if (fmpz_mpoly_q_factored_is_zero(x, mctx))
    {
        if (fmpz_sgn(e) < 0)
            return GR_DOMAIN;

        if (fmpz_is_zero(e))
            fmpz_mpoly_q_factored_one(res, mctx);
        else
            fmpz_mpoly_q_factored_zero(res, mctx);

        return GR_SUCCESS;
    }

    /* +/- 1 */
    if (x->den.num == 0 && fmpz_is_one(x->den.constant) &&
        fmpz_mpoly_is_fmpz(fmpz_mpoly_q_factored_numref(x), mctx) &&
        fmpz_is_pm1(x->num.coeffs))
    {
        if (fmpz_is_one(x->num.coeffs) || fmpz_is_even(e))
            fmpz_mpoly_q_factored_one(res, mctx);
        else
            fmpz_mpoly_q_factored_set_si(res, -1, mctx);

        return GR_SUCCESS;
    }

    /* Anything else would certainly exhaust memory. */
    if (!fmpz_fits_si(e) || FLINT_ABS(fmpz_get_si(e)) > WORD(1) << (FLINT_BITS / 2))
        return GR_UNABLE;

    fmpz_mpoly_q_factored_pow_fmpz(res, x, e, mctx);
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_pow_ui(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, ulong e, gr_ctx_t ctx)
{
    fmpz_t f;
    int status;
    fmpz_init_set_ui(f, e);
    status = _gr_fmpz_mpoly_q_factored_pow_fmpz(res, x, f, ctx);
    fmpz_clear(f);
    return status;
}

static int
_gr_fmpz_mpoly_q_factored_pow_si(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, slong e, gr_ctx_t ctx)
{
    fmpz_t f;
    int status;
    fmpz_init_set_si(f, e);
    status = _gr_fmpz_mpoly_q_factored_pow_fmpz(res, x, f, ctx);
    fmpz_clear(f);
    return status;
}

static int
_gr_fmpz_mpoly_q_factored_numerator(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const gr_ctx_t ctx)
{
    fmpz_mpoly_q_factored_set_fmpz_mpoly(res, fmpz_mpoly_q_factored_numref(x), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

static int
_gr_fmpz_mpoly_q_factored_denominator(fmpz_mpoly_q_factored_t res, const fmpz_mpoly_q_factored_t x, const gr_ctx_t ctx)
{
    fmpz_mpoly_t t;
    int status = GR_SUCCESS;

    fmpz_mpoly_init(t, MPOLYNOMIAL_MCTX(ctx));
    if (fmpz_mpoly_factor_expand(t, fmpz_mpoly_q_factored_denref(x), MPOLYNOMIAL_MCTX(ctx)))
        fmpz_mpoly_q_factored_set_fmpz_mpoly(res, t, MPOLYNOMIAL_MCTX(ctx));
    else
        status = GR_UNABLE;
    fmpz_mpoly_clear(t, MPOLYNOMIAL_MCTX(ctx));

    return status;
}

int _gr_fmpz_mpoly_q_factored_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mpoly_q_factored_methods;

gr_method_tab_input _gr_fmpz_mpoly_q_factored_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE,        (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE,     (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_SET_GEN_NAMES,               (gr_funcptr) _gr_fmpz_mpoly_q_factored_ctx_set_gen_names},
    {GR_METHOD_CTX_NGENS,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_ctx_ngens},
    {GR_METHOD_CTX_GEN_NAME, (gr_funcptr) _gr_fmpz_mpoly_q_factored_ctx_gen_name},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mpoly_q_factored_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mpoly_q_factored_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpz_mpoly_q_factored_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_randtest_small},
    {_GR_METHOD_LENGTH,     (gr_funcptr) _gr_fmpz_mpoly_q_factored_length},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpz_mpoly_q_factored_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpz_mpoly_q_factored_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpz_mpoly_q_factored_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_is_one},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_fmpz_mpoly_q_factored_gens},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_fmpz_mpoly_q_factored_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_fmpq},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_set_other},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_add},
    {GR_METHOD_ADD_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_add_si},
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_add_fmpq},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_sub},
    {GR_METHOD_SUB_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_sub_si},
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_sub_fmpq},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_mul},
    {GR_METHOD_MUL_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_mul_si},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_mul_fmpq},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_div},
    {GR_METHOD_DIV_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_si},
    {GR_METHOD_DIV_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_fmpz},
    {GR_METHOD_DIV_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_fmpq},
    {GR_METHOD_DIVEXACT,         (gr_funcptr) _gr_fmpz_mpoly_q_factored_div},
    {GR_METHOD_DIVEXACT_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_si},
    {GR_METHOD_DIVEXACT_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_fmpz},
    {GR_METHOD_DIVEXACT_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_div_fmpq},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mpoly_q_factored_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_is_invertible},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_pow_ui},
    {GR_METHOD_POW_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_factored_pow_si},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_factored_pow_fmpz},
    {GR_METHOD_NUMERATOR,   (gr_funcptr) _gr_fmpz_mpoly_q_factored_numerator},
    {GR_METHOD_DENOMINATOR, (gr_funcptr) _gr_fmpz_mpoly_q_factored_denominator},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_mpoly_q_factored(gr_ctx_t ctx, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_FMPZ_MPOLY_Q_FACTORED;
    ctx->sizeof_elem = sizeof(fmpz_mpoly_q_factored_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);
    MPOLYNOMIAL_CTX(ctx)->vars = NULL;

    ctx->methods = _gr_fmpz_mpoly_q_factored_methods;

    if (!_gr_fmpz_mpoly_q_factored_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mpoly_q_factored_methods, _gr_fmpz_mpoly_q_factored_methods_input);
        _gr_fmpz_mpoly_q_factored_methods_initialized = 1;
    }
}
