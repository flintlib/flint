/*
    Copyright (C) 2023 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include "gr.h"
#include "gr_vec.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "fmpz_mpoly_factor.h"

typedef struct
{
    fmpz_mpoly_ctx_t mctx;
}
_gr_fmpz_mpoly_ctx_t;

#define MPOLYNOMIAL_CTX(ring_ctx) ((_gr_fmpz_mpoly_ctx_t *)(GR_CTX_DATA_AS_PTR(ring_ctx)))
#define MPOLYNOMIAL_MCTX(ring_ctx) (MPOLYNOMIAL_CTX(ring_ctx)->mctx)

int _gr_fmpz_mpoly_q_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Fraction field of multivariate polynomials over Integer ring (fmpz)");
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, MPOLYNOMIAL_MCTX(ctx)->minfo->nvars);
    gr_stream_write(out, " variables");
    if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (MPOLYNOMIAL_MCTX(ctx)->minfo->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void
_gr_fmpz_mpoly_q_ctx_clear(gr_ctx_t ctx)
{
    fmpz_mpoly_ctx_clear(MPOLYNOMIAL_MCTX(ctx));
    flint_free(GR_CTX_DATA_AS_PTR(ctx));
}

void
_gr_fmpz_mpoly_q_init(fmpz_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_init(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_q_clear(fmpz_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_clear(res, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_q_swap(fmpz_mpoly_q_t poly1, fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mpoly_q_swap(poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
}

void
_gr_fmpz_mpoly_q_set_shallow(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly, gr_ctx_t ctx)
{
    *res = *poly;
}

int
_gr_fmpz_mpoly_q_randtest(fmpz_mpoly_q_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong bits;

    if (n_randint(state, 10) != 0)
        bits = 10;
    else
        bits = 100;

    fmpz_mpoly_q_randtest(res, state, n_randint(state, 5), bits, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_randtest_small(fmpz_mpoly_q_t res, flint_rand_t state, gr_ctx_t ctx)
{
    fmpz_mpoly_q_randtest(res, state, n_randint(state, 3), 3, 1 + n_randint(state, 3), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_write(gr_stream_t out, fmpz_mpoly_q_t f, gr_ctx_t ctx)
{
    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(f), MPOLYNOMIAL_MCTX(ctx)))
    {
        gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), NULL, MPOLYNOMIAL_MCTX(ctx)));
    }
    else if (fmpz_mpoly_is_fmpz(fmpz_mpoly_q_denref(f), MPOLYNOMIAL_MCTX(ctx)))
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), NULL, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")/");
        gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(f), NULL, MPOLYNOMIAL_MCTX(ctx)));
    }
    else
    {
        gr_stream_write(out, "(");
        gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(f), NULL, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")/(");
        gr_stream_write_free(out, fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(f), NULL, MPOLYNOMIAL_MCTX(ctx)));
        gr_stream_write(out, ")");
    }

    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_q_equal(const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_equal(poly1, poly2, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mpoly_q_is_zero(const fmpz_mpoly_q_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_is_zero(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fmpz_mpoly_q_is_one(const fmpz_mpoly_q_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_is_one(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mpoly_q_zero(fmpz_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_zero(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_one(fmpz_mpoly_q_t res, gr_ctx_t ctx)
{
    fmpz_mpoly_q_one(res, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_q_gens(gr_vec_t res, gr_ctx_t ctx)
{
    slong i, n;

    n = MPOLYNOMIAL_MCTX(ctx)->minfo->nvars;

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        fmpz_mpoly_q_gen(((fmpz_mpoly_q_struct *) res->entries) + i, i, MPOLYNOMIAL_MCTX(ctx));

    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t mat, gr_ctx_t ctx)
{
    fmpz_mpoly_q_set(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_set_si(fmpz_mpoly_q_t res, slong v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_set_si(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_set_ui(fmpz_mpoly_q_t res, ulong v, gr_ctx_t ctx)
{
    fmpz_mpoly_set_ui(fmpz_mpoly_q_numref(res), v, MPOLYNOMIAL_MCTX(ctx));
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_set_fmpz(fmpz_mpoly_q_t res, const fmpz_t v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_set_fmpz(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_set_fmpq(fmpz_mpoly_q_t res, const fmpq_t v, gr_ctx_t ctx)
{
    fmpz_mpoly_q_set_fmpq(res, v, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t mat, gr_ctx_t ctx)
{
    fmpz_mpoly_q_neg(res, mat, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mpoly_q_add(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_add_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_add_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_fmpz_mpoly_q_add_ui(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_add_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_fmpz_mpoly_q_add_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_add_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_add_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_add_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mpoly_q_sub(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_sub_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_sub_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_fmpz_mpoly_q_sub_ui(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_sub_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_fmpz_mpoly_q_sub_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_sub_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_sub_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_sub_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    fmpz_mpoly_q_mul(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_mul_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_mul_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_fmpz_mpoly_q_mul_ui(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_mul_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_fmpz_mpoly_q_mul_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_mul_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_mul_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    fmpz_mpoly_q_mul_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    if (fmpz_mpoly_q_is_zero(poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_DOMAIN;

    fmpz_mpoly_q_div(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_div_si(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, slong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpz_mpoly_q_div_si(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

/*
int
_gr_fmpz_mpoly_q_div_ui(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    if (c == 0)
        return GR_DOMAIN;

    fmpz_mpoly_q_div_ui(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}
*/

int
_gr_fmpz_mpoly_q_div_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_is_zero(c))
        return GR_DOMAIN;

    fmpz_mpoly_q_div_fmpz(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_div_fmpq(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpq_t c, gr_ctx_t ctx)
{
    if (fmpq_is_zero(c))
        return GR_DOMAIN;

    fmpz_mpoly_q_div_fmpq(res, poly1, c, MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fmpz_mpoly_q_is_invertible(const fmpz_mpoly_q_t c, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_is_zero(c, MPOLYNOMIAL_MCTX(ctx)) ? T_FALSE : T_TRUE;
}

int
_gr_fmpz_mpoly_q_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t c, gr_ctx_t ctx)
{
    if (!fmpz_mpoly_q_is_zero(c, MPOLYNOMIAL_MCTX(ctx)))
    {
        fmpz_mpoly_q_inv(res, c, MPOLYNOMIAL_MCTX(ctx));
        return GR_SUCCESS;
    }

    return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_q_pow_ui(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, ulong c, gr_ctx_t ctx)
{
    if (fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(poly1), c, MPOLYNOMIAL_MCTX(ctx)) &&
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(poly1), c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;
}

int
_gr_fmpz_mpoly_q_pow_fmpz(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_t c, gr_ctx_t ctx)
{
    if (fmpz_sgn(c) < 0)
    {
        int status;

        status = gr_inv(res, poly1, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_t e;
            fmpz_init(e);
            fmpz_neg(e, c);
            status = _gr_fmpz_mpoly_q_pow_fmpz(res, res, e, ctx);
            fmpz_clear(e);
        }

        return status;
    }

    if (fmpz_mpoly_pow_fmpz(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(poly1), c, MPOLYNOMIAL_MCTX(ctx)) &&
        fmpz_mpoly_pow_fmpz(fmpz_mpoly_q_denref(res), fmpz_mpoly_q_denref(poly1), c, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_UNABLE;

}

int
_gr_fmpz_mpoly_q_numerator(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const gr_ctx_t ctx)
{
    fmpz_mpoly_set(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_numref(x), MPOLYNOMIAL_MCTX(ctx));
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fmpz_mpoly_q_denominator(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const gr_ctx_t ctx)
{
    fmpz_mpoly_set(fmpz_mpoly_q_numref(res), fmpz_mpoly_q_denref(x), MPOLYNOMIAL_MCTX(ctx));
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), MPOLYNOMIAL_MCTX(ctx));
    return GR_SUCCESS;
}


/*
truth_t
_gr_fmpz_mpoly_q_is_square(const fmpz_mpoly_q_t poly, gr_ctx_t ctx)
{
    return fmpz_mpoly_q_is_square(poly, MPOLYNOMIAL_MCTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fmpz_mpoly_q_sqrt(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly, gr_ctx_t ctx)
{
    if (fmpz_mpoly_q_sqrt(res, poly, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;
    else
        return GR_DOMAIN;
}

int
_gr_fmpz_mpoly_q_gcd(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t poly1, const fmpz_mpoly_q_t poly2, gr_ctx_t ctx)
{
    if (fmpz_mpoly_q_gcd(res, poly1, poly2, MPOLYNOMIAL_MCTX(ctx)))
        return GR_SUCCESS;

    return GR_DOMAIN;
}
*/

/*
int
_gr_fmpz_mpoly_q_factor(fmpz_mpoly_q_t c, gr_vec_t factors, gr_vec_t exponents, gr_srcptr x, int flags, gr_ctx_t ctx)
{
    fmpz_mpoly_q_factor_t fac;
    gr_ctx_t ZZ;
    slong i;
    int status = GR_SUCCESS;

    fmpz_mpoly_q_factor_init(fac, MPOLYNOMIAL_MCTX(ctx));

    if (fmpz_mpoly_q_factor(fac, x, MPOLYNOMIAL_MCTX(ctx)))
    {
        fmpz_mpoly_q_set_fmpz(c, fac->constant, MPOLYNOMIAL_MCTX(ctx));

        gr_ctx_init_fmpz(ZZ);

        gr_vec_set_length(factors, fac->num, ctx);
        gr_vec_set_length(exponents, fac->num, ZZ);

        for (i = 0; i < fac->num; i++)
        {
            fmpz_mpoly_q_swap((fmpz_mpoly_q_struct *) (factors->entries) + i, fac->poly + i, MPOLYNOMIAL_MCTX(ctx));
            fmpz_swap((fmpz *) (exponents->entries) + i, fac->exp + i);
        }

        gr_ctx_clear(ZZ);
    }
    else
    {
        status = GR_UNABLE;
    }

    fmpz_mpoly_q_factor_clear(fac, MPOLYNOMIAL_MCTX(ctx));

    return status;
}
*/

int _gr_fmpz_mpoly_q_methods_initialized = 0;

gr_static_method_table _gr_fmpz_mpoly_q_methods;

gr_method_tab_input _gr_fmpz_mpoly_q_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) _gr_fmpz_mpoly_q_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) _gr_fmpz_mpoly_q_ctx_clear},
    {GR_METHOD_CTX_IS_RING,                     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING,         (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,          (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,                    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,                   (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,    (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_THREADSAFE,               (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,        (gr_funcptr) _gr_fmpz_mpoly_q_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) _gr_fmpz_mpoly_q_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) _gr_fmpz_mpoly_q_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) _gr_fmpz_mpoly_q_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_fmpz_mpoly_q_randtest},
    {GR_METHOD_RANDTEST_SMALL,    (gr_funcptr) _gr_fmpz_mpoly_q_randtest_small},
    {GR_METHOD_WRITE,       (gr_funcptr) _gr_fmpz_mpoly_q_write},
    {GR_METHOD_ZERO,        (gr_funcptr) _gr_fmpz_mpoly_q_zero},
    {GR_METHOD_ONE,         (gr_funcptr) _gr_fmpz_mpoly_q_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) _gr_fmpz_mpoly_q_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) _gr_fmpz_mpoly_q_is_one},
    {GR_METHOD_GENS,        (gr_funcptr) _gr_fmpz_mpoly_q_gens},
    {GR_METHOD_EQUAL,       (gr_funcptr) _gr_fmpz_mpoly_q_equal},
    {GR_METHOD_SET,         (gr_funcptr) _gr_fmpz_mpoly_q_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_set_fmpq},
    {GR_METHOD_NEG,         (gr_funcptr) _gr_fmpz_mpoly_q_neg},
    {GR_METHOD_ADD,         (gr_funcptr) _gr_fmpz_mpoly_q_add},
    {GR_METHOD_ADD_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_add_si},
/*    {GR_METHOD_ADD_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_add_ui}, */
    {GR_METHOD_ADD_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_add_fmpz},
    {GR_METHOD_ADD_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_add_fmpq},
    {GR_METHOD_SUB,         (gr_funcptr) _gr_fmpz_mpoly_q_sub},
    {GR_METHOD_SUB_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_sub_si},
/*    {GR_METHOD_SUB_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_sub_ui}, */
    {GR_METHOD_SUB_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_sub_fmpz},
    {GR_METHOD_SUB_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_sub_fmpq},
    {GR_METHOD_MUL,         (gr_funcptr) _gr_fmpz_mpoly_q_mul},
    {GR_METHOD_MUL_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_mul_si},
/*    {GR_METHOD_MUL_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_mul_ui}, */
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_mul_fmpq},
    {GR_METHOD_DIV,         (gr_funcptr) _gr_fmpz_mpoly_q_div},
    {GR_METHOD_DIV_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_div_si},
/*    {GR_METHOD_DIV_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_div_ui}, */
    {GR_METHOD_DIV_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_div_fmpz},
    {GR_METHOD_DIV_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_div_fmpq},
    {GR_METHOD_DIVEXACT,         (gr_funcptr) _gr_fmpz_mpoly_q_div},
    {GR_METHOD_DIVEXACT_SI,      (gr_funcptr) _gr_fmpz_mpoly_q_div_si},
/*    {GR_METHOD_DIVEXACT_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_div_ui}, */
    {GR_METHOD_DIVEXACT_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_div_fmpz},
    {GR_METHOD_DIVEXACT_FMPQ,    (gr_funcptr) _gr_fmpz_mpoly_q_div_fmpq},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fmpz_mpoly_q_inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fmpz_mpoly_q_is_invertible},
    {GR_METHOD_POW_UI,      (gr_funcptr) _gr_fmpz_mpoly_q_pow_ui},
    {GR_METHOD_POW_FMPZ,    (gr_funcptr) _gr_fmpz_mpoly_q_pow_fmpz},
    {GR_METHOD_NUMERATOR,      (gr_funcptr) _gr_fmpz_mpoly_q_numerator},
    {GR_METHOD_DENOMINATOR,      (gr_funcptr) _gr_fmpz_mpoly_q_denominator},

/*
    {GR_METHOD_SQRT,        (gr_funcptr) _gr_fmpz_mpoly_q_sqrt},
    {GR_METHOD_IS_SQUARE,   (gr_funcptr) _gr_fmpz_mpoly_q_is_square},
    {GR_METHOD_GCD,         (gr_funcptr) _gr_fmpz_mpoly_q_gcd},
    {GR_METHOD_FACTOR,      (gr_funcptr) _gr_fmpz_mpoly_q_factor},
*/
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_fmpz_mpoly_q(gr_ctx_t ctx, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_FMPZ_MPOLY_Q;
    ctx->sizeof_elem = sizeof(fmpz_mpoly_q_struct);
    GR_CTX_DATA_AS_PTR(ctx) = flint_malloc(sizeof(_gr_fmpz_mpoly_ctx_t));
    ctx->size_limit = WORD_MAX;

    fmpz_mpoly_ctx_init(MPOLYNOMIAL_MCTX(ctx), nvars, ord);

    ctx->methods = _gr_fmpz_mpoly_q_methods;

    if (!_gr_fmpz_mpoly_q_methods_initialized)
    {
        gr_method_tab_init(_gr_fmpz_mpoly_q_methods, _gr_fmpz_mpoly_q_methods_input);
        _gr_fmpz_mpoly_q_methods_initialized = 1;
    }
}
