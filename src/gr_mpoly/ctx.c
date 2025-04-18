/*
    Copyright (C) 2022, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "mpoly.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_mpoly.h"

int gr_mpoly_ctx_write(gr_stream_t out, gr_mpoly_ctx_t ctx)
{
    gr_stream_write(out, "Ring of multivariate polynomials over ");
    gr_ctx_write(out, GR_MPOLY_CCTX(ctx));
    gr_stream_write(out, " in ");
    gr_stream_write_si(out, GR_MPOLY_NVARS(ctx));
    gr_stream_write(out, " variables");
    if (GR_MPOLY_MCTX(ctx)->ord == ORD_LEX)
        gr_stream_write(out, ", lex order");
    else if (GR_MPOLY_MCTX(ctx)->ord == ORD_DEGLEX)
        gr_stream_write(out, ", deglex order");
    else if (GR_MPOLY_MCTX(ctx)->ord == ORD_DEGREVLEX)
        gr_stream_write(out, ", degrevlex order");
    return GR_SUCCESS;
}

void
gr_mpoly_ctx_clear(gr_mpoly_ctx_t ctx)
{
    if (GR_MPOLY_VARS(ctx) != NULL)
    {
        slong i;
        for (i = 0; i < GR_MPOLY_NVARS(ctx); i++)
            flint_free(GR_MPOLY_VARS(ctx)[i]);
        flint_free(GR_MPOLY_VARS(ctx));
    }

    mpoly_ctx_clear(GR_MPOLY_MCTX(ctx));
    flint_free(GR_MPOLY_MCTX(ctx));
}

int
gr_mpoly_ctx_set_gen_names(gr_mpoly_ctx_t ctx, const char ** s)
{
    slong i, nvars, len;

    nvars = GR_MPOLY_NVARS(ctx);

    if (GR_MPOLY_VARS(ctx) == NULL)
    {
        GR_MPOLY_VARS(ctx) = flint_malloc(nvars * sizeof(char *));
        for (i = 0; i < nvars; i++)
            GR_MPOLY_VARS(ctx)[i] = NULL;
    }

    for (i = 0; i < nvars; i++)
    {
        len = strlen(s[i]);
        GR_MPOLY_VARS(ctx)[i] = flint_realloc(GR_MPOLY_VARS(ctx)[i], len + 1);
        memcpy(GR_MPOLY_VARS(ctx)[i], s[i], len + 1);
    }

    return GR_SUCCESS;
}

truth_t
gr_mpoly_ctx_is_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_zero_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_zero_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_commutative_ring(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_integral_domain(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(GR_MPOLY_CCTX(ctx));
}

truth_t
gr_mpoly_ctx_is_field(gr_mpoly_ctx_t ctx)
{
    if (GR_MPOLY_NVARS(ctx) == 0)
        return gr_ctx_is_field(GR_MPOLY_CCTX(ctx));
    else
        return T_FALSE;
}

truth_t
gr_mpoly_ctx_is_threadsafe(gr_mpoly_ctx_t ctx)
{
    return gr_ctx_is_threadsafe(GR_MPOLY_CCTX(ctx));
}

int
gr_mpoly_gens(gr_vec_t res, gr_mpoly_ctx_t ctx)
{
    slong i, n;
    int status = GR_SUCCESS;

    n = GR_MPOLY_NVARS(ctx);

    gr_vec_set_length(res, n, ctx);
    for (i = 0; i < n; i++)
        status |= gr_mpoly_gen(((gr_mpoly_struct *) res->entries) + i, i, ctx);

    return status;
}

int
gr_mpoly_gens_recursive(gr_vec_t vec, gr_mpoly_ctx_t ctx)
{
    int status;
    gr_vec_t vec1;
    slong i, n, m;

    /* Get generators of the element ring */
    gr_vec_init(vec1, 0, GR_MPOLY_CCTX(ctx));
    status = gr_gens_recursive(vec1, GR_MPOLY_CCTX(ctx));
    n = vec1->length;

    m = GR_MPOLY_NVARS(ctx);

    gr_vec_set_length(vec, n + m, ctx);

    /* Promote to polynomials */
    for (i = 0; i < n; i++)
        status |= gr_mpoly_set_scalar(gr_vec_entry_ptr(vec, i, ctx),
                gr_vec_entry_srcptr(vec1, i, GR_MPOLY_CCTX(ctx)), ctx);

    for (i = 0; i < m; i++)
        status |= gr_mpoly_gen(((gr_mpoly_struct *) vec->entries) + n + i, i, ctx);

    gr_vec_clear(vec1, GR_MPOLY_CCTX(ctx));

    return status;
}

/* FIXME: this may inappropriately return GR_DOMAIN for nonconstant
   polynomials non-integral domains. See AbstractAlgebra. */
int
gr_mpoly_inv(gr_mpoly_t res, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    if (poly->length == 0)
    {
        if (gr_ctx_is_zero_ring(GR_MPOLY_CCTX(ctx)) == T_TRUE)
            return gr_mpoly_zero(res, ctx);
        else
            return GR_DOMAIN;
    }
    else if (poly->length == 1)
    {
        slong N;
        gr_ptr c;
        int status;

        N = mpoly_words_per_exp(poly->bits, GR_MPOLY_MCTX(ctx));
        if (!mpoly_monomial_is_zero(poly->exps + N*0, N))
            return GR_DOMAIN;

        /* todo: avoid the temporary */
        GR_TMP_INIT(c, GR_MPOLY_CCTX(ctx));
        status = gr_inv(c, poly->coeffs, GR_MPOLY_CCTX(ctx));
        status |= gr_mpoly_set_scalar(res, c, ctx);
        GR_TMP_CLEAR(c, GR_MPOLY_CCTX(ctx));
        return status;
    }
    else
    {
        if (gr_is_zero(poly->coeffs, GR_MPOLY_CCTX(ctx)) == T_FALSE)
            return GR_DOMAIN;
        else
            return GR_UNABLE;
    }
}


static int _gr_mpoly_remove_zeros(gr_mpoly_t A, gr_mpoly_ctx_t ctx)
{
    mpoly_ctx_struct * mctx = GR_MPOLY_MCTX(ctx);
    gr_ctx_struct * cctx = GR_MPOLY_CCTX(ctx);
    slong i, N;
    slong Alen, Blen;
    ulong * Aexp;
    fmpz * Acoeff;
    int status = GR_SUCCESS;
    slong sz = cctx->sizeof_elem;

    N = mpoly_words_per_exp(A->bits, mctx);

    Blen = A->length;
    Aexp = A->exps;
    Acoeff = A->coeffs;

    Alen = 0;
    for (i = 0; i < Blen; i++)
    {
        if (i != Alen)
            mpoly_monomial_set(Aexp + N*Alen, Aexp + N*i, N);

        Alen += (gr_is_zero(GR_ENTRY(Acoeff, Alen, sz), cctx) != T_TRUE);
    }

    A->length = Alen;

    return status;
}

int
gr_mpoly_canonical_associate(gr_mpoly_t res, gr_mpoly_t u, const gr_mpoly_t poly, gr_mpoly_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong len = poly->length;

    if (len == 0)
    {
        status = gr_mpoly_zero(res, ctx);
        if (u != NULL)
            status |= gr_mpoly_one(u, ctx);
    }
    else
    {
        gr_ptr c;

        if (res != poly)
            status |= gr_mpoly_set(res, poly, ctx);

        FLINT_ASSERT(len == res->length);

        GR_TMP_INIT(c, GR_MPOLY_CCTX(ctx));

        status |= gr_canonical_associate(res->coeffs, c, res->coeffs, GR_MPOLY_CCTX(ctx));
        status |= _gr_vec_mul_scalar(GR_ENTRY(res->coeffs, 1, GR_MPOLY_CCTX(ctx)->sizeof_elem),
                                     GR_ENTRY(res->coeffs, 1, GR_MPOLY_CCTX(ctx)->sizeof_elem),
                                     len - 1, c, GR_MPOLY_CCTX(ctx));
        /* todo: can skip over some rings */
        status |= _gr_mpoly_remove_zeros(res, ctx);

        if (u != NULL)
            status |= gr_mpoly_set_scalar(u, c, ctx);

        GR_TMP_CLEAR(c, GR_MPOLY_CCTX(ctx));
    }

    return status;
}


int _gr_mpoly_methods_initialized = 0;

gr_static_method_table _gr_mpoly_methods;

gr_method_tab_input _gr_mpoly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) gr_mpoly_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) gr_mpoly_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_mpoly_ctx_is_ring},
    {GR_METHOD_CTX_IS_ZERO_RING,     (gr_funcptr) gr_mpoly_ctx_is_zero_ring},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_mpoly_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_mpoly_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_mpoly_ctx_is_field},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) gr_mpoly_ctx_is_threadsafe},
    {GR_METHOD_CTX_SET_GEN_NAMES,       (gr_funcptr) gr_mpoly_ctx_set_gen_names},
    {GR_METHOD_INIT,        (gr_funcptr) gr_mpoly_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) gr_mpoly_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) gr_mpoly_swap},
    {GR_METHOD_SET_SHALLOW, (gr_funcptr) gr_mpoly_set_shallow},
    {GR_METHOD_RANDTEST,    (gr_funcptr) _gr_mpoly_randtest_default},
    {_GR_METHOD_LENGTH,     (gr_funcptr) gr_mpoly_length},
    {GR_METHOD_WRITE,       (gr_funcptr) gr_mpoly_write},
    {GR_METHOD_GENS,        (gr_funcptr) gr_mpoly_gens},
    {GR_METHOD_GENS_RECURSIVE,       (gr_funcptr) gr_mpoly_gens_recursive},
    {GR_METHOD_ZERO,        (gr_funcptr) gr_mpoly_zero},
    {GR_METHOD_ONE,         (gr_funcptr) gr_mpoly_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) gr_mpoly_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) gr_mpoly_is_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) gr_mpoly_equal},
    {GR_METHOD_SET,         (gr_funcptr) gr_mpoly_set},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) gr_mpoly_set_other},
    {GR_METHOD_SET_UI,      (gr_funcptr) gr_mpoly_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) gr_mpoly_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) gr_mpoly_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) gr_mpoly_set_fmpq},
    {GR_METHOD_SET_STR,     (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,         (gr_funcptr) gr_mpoly_neg},
    {GR_METHOD_ADD,         (gr_funcptr) gr_mpoly_add},
    {GR_METHOD_SUB,         (gr_funcptr) gr_mpoly_sub},
    {GR_METHOD_MUL,         (gr_funcptr) gr_mpoly_mul},
    {GR_METHOD_MUL_UI,      (gr_funcptr) gr_mpoly_mul_ui},
    {GR_METHOD_MUL_SI,      (gr_funcptr) gr_mpoly_mul_si},
    {GR_METHOD_MUL_FMPZ,    (gr_funcptr) gr_mpoly_mul_fmpz},
    {GR_METHOD_MUL_FMPQ,    (gr_funcptr) gr_mpoly_mul_fmpq},
    {GR_METHOD_INV,         (gr_funcptr) gr_mpoly_inv},
    {GR_METHOD_CANONICAL_ASSOCIATE,         (gr_funcptr) gr_mpoly_canonical_associate},
    {0,                     (gr_funcptr) NULL},
};

/* todo: first arg as gr_mpoly_ctx_t */
void
gr_mpoly_ctx_init(gr_mpoly_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord)
{
    ctx->which_ring = GR_CTX_GR_MPOLY;
    ctx->sizeof_elem = sizeof(gr_mpoly_struct);
    ctx->size_limit = WORD_MAX;

    /* by reference */
    GR_MPOLY_CCTX(ctx) = base_ring;

    /* allocated here */
    GR_MPOLY_MCTX(ctx) = flint_malloc(sizeof(mpoly_ctx_struct));
    mpoly_ctx_init(GR_MPOLY_MCTX(ctx), nvars, ord);

    GR_MPOLY_VARS(ctx) = NULL;

    ctx->methods = _gr_mpoly_methods;

    if (!_gr_mpoly_methods_initialized)
    {
        gr_method_tab_init(_gr_mpoly_methods, _gr_mpoly_methods_input);
        _gr_mpoly_methods_initialized = 1;
    }
}

void
gr_mpoly_ctx_init_rand(gr_mpoly_ctx_t ctx, flint_rand_t state, gr_ctx_t base_ring, slong max_nvars)
{
    gr_ctx_init_gr_mpoly(ctx, base_ring, n_randint(state, max_nvars + 1), mpoly_ordering_randtest(state));
}
