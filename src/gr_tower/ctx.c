/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    gr context for the top field F_n of a tower. Elements have the same
    representation as elements of the top quotient ring context; most
    operations are delegated to it. The zero test, equality test and
    inversion go through the dynamic refinement machinery, which makes
    them complete: the field is a genuine field, and the tower is refined
    whenever an operation exposes a factorization of a defining polynomial.
*/

#include <string.h>
#include "fmpq.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "gr_tower.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

/* the context data: the tower, and the level (number of algebraic
   generators) and number of transcendental generators the context was
   made for; a context stays valid when algebraic steps are appended to
   the tower (its elements then live in a lower field, which persists),
   but not when the tower is rebuilt with a new base field; the nested
   context its elements belong to is recorded, so that they can still be
   cleared when the context has become stale (the tower keeps superseded
   nested contexts until it is cleared) */
typedef struct
{
    gr_tower_struct * T;
    slong level;
    slong num_trans;
    ulong structure_version;
    gr_ctx_struct * top;
}
tower_field_ctx_struct;

#define TCTX(ctx) ((tower_field_ctx_struct *) (ctx)->data)
/* the context is stale after a rebuild of the tower with a new base
   field, or a change of the order of its generators */
#define VALID(ctx) (TOWER(ctx)->num_trans == TCTX(ctx)->num_trans && TOWER(ctx)->structure_version == TCTX(ctx)->structure_version)
#define TOWER(ctx) (TCTX(ctx)->T)
#define TOP(ctx) gr_tower_field_at(TOWER(ctx), TCTX(ctx)->level)
#define LEVEL(ctx) (TCTX(ctx)->level)

static int
_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k;
    gr_tower_struct * T = TOWER(ctx);

    status |= gr_stream_write(out, "Tower field ");
    status |= gr_ctx_write(out, T->base);

    for (k = 1; k <= LEVEL(ctx); k++)
    {
        status |= gr_stream_write(out, (k == 1) ? "(" : ", ");
        status |= gr_stream_write(out, GR_TOWER_STEP(T, k - 1)->name);
    }

    if (LEVEL(ctx) > 0)
        status |= gr_stream_write(out, ")");

    return status;
}

static void _ctx_clear(gr_ctx_t ctx) { }

static truth_t _ctx_is_field(gr_ctx_t ctx) { return gr_ctx_is_field(TOWER(ctx)->base); }
static truth_t _ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(TOWER(ctx)->base); }
static truth_t _ctx_is_rational_vector_space(gr_ctx_t ctx) { return gr_ctx_is_rational_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_real_vector_space(gr_ctx_t ctx) { return gr_ctx_is_real_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_complex_vector_space(gr_ctx_t ctx) { return gr_ctx_is_complex_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(TOWER(ctx)->base); }

/* Representations are canonical only relative to the current tower;
   is_zero and equal always canonicalize, so the field behaves canonically. */
static truth_t _ctx_is_canonical(gr_ctx_t ctx) { return T_FALSE; }

/* Refinement mutates the tower; concurrent use requires external locking for now. */
static truth_t _ctx_is_threadsafe(gr_ctx_t ctx) { return T_FALSE; }

static gr_ptr _ctx_base(gr_ctx_t ctx) { return TOWER(ctx)->base; }

static int
_ctx_ngens(slong * ngens, gr_ctx_t ctx)
{
    *ngens = LEVEL(ctx);
    return GR_SUCCESS;
}

static int
_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    gr_tower_struct * T = TOWER(ctx);
    size_t len;

    if (i < 0 || i >= LEVEL(ctx))
        return GR_DOMAIN;

    len = strlen(GR_TOWER_STEP(T, i)->name);
    *name = flint_malloc(len + 1);
    memcpy(*name, GR_TOWER_STEP(T, i)->name, len + 1);
    return GR_SUCCESS;
}

static void _init(gr_ptr x, gr_ctx_t ctx) { if (!VALID(ctx)) flint_throw(FLINT_ERROR, "gr_tower_field: the tower was rebuilt (transcendental generator adjoined, or generators reordered) after this context was created\n"); gr_init(x, TOP(ctx)); }
static void
_clear(gr_ptr x, gr_ctx_t ctx)
{
    if (VALID(ctx))
        gr_clear(x, TOP(ctx));
    else if (TOWER(ctx)->keep_retired)
        gr_clear(x, TCTX(ctx)->top);   /* (a superseded context, kept by the tower) */
    /* (else: a stale element is leaked rather than misinterpreted) */
}
static void _swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) { gr_swap(x, y, TOP(ctx)); }
static void _set_shallow(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) { gr_set_shallow(x, y, TOP(ctx)); }
static int _randtest(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_randtest(x, state, TOP(ctx)); }
static int _write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_write(out, x, TOP(ctx)); }
static int _zero(gr_ptr x, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_zero(x, TOP(ctx)); }
static int _one(gr_ptr x, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_one(x, TOP(ctx)); }
static int _set(gr_ptr x, gr_srcptr y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_set(x, y, TOP(ctx)); }
static int _set_si(gr_ptr x, slong y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_set_si(x, y, TOP(ctx)); }
static int _set_ui(gr_ptr x, ulong y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_set_ui(x, y, TOP(ctx)); }
static int _set_fmpz(gr_ptr x, const fmpz_t y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_set_fmpz(x, y, TOP(ctx)); }
static int _set_fmpq(gr_ptr x, const fmpq_t y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_set_fmpq(x, y, TOP(ctx)); }
static int _neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_neg(res, x, TOP(ctx)); }
static int _add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_add(res, x, y, TOP(ctx)); }
static int _sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_sub(res, x, y, TOP(ctx)); }
static int _mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_mul(res, x, y, TOP(ctx)); }
static int _sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_sqr(res, x, TOP(ctx)); }
static int _mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_mul_si(res, x, y, TOP(ctx)); }
static int _mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_mul_fmpz(res, x, y, TOP(ctx)); }
static int _mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { if (!VALID(ctx)) return GR_UNABLE; return gr_mul_fmpq(res, x, y, TOP(ctx)); }

static int
_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (!VALID(ctx) || (x_ctx->which_ring == GR_CTX_GR_TOWER_FIELD && !VALID(x_ctx)))
        return GR_UNABLE;
    if (x_ctx == ctx)
        return gr_set(res, x, TOP(ctx));
    if (x_ctx->which_ring == GR_CTX_GR_TOWER_FIELD && TOWER(x_ctx) == TOWER(ctx) && LEVEL(x_ctx) <= LEVEL(ctx))
        return gr_tower_promote(res, x, LEVEL(x_ctx), LEVEL(ctx), TOWER(ctx));
    return gr_set_other(res, x, x_ctx, TOP(ctx));
}

static int
_gens(gr_vec_t vec, gr_ctx_t ctx)
{
    gr_tower_struct * T = TOWER(ctx);
    gr_ctx_struct * top = TOP(ctx);
    slong k;
    int status = GR_SUCCESS;
    gr_ptr t;

    if (!VALID(ctx))
        return GR_UNABLE;

    gr_vec_set_length(vec, LEVEL(ctx), ctx);

    for (k = 1; k <= LEVEL(ctx) && status == GR_SUCCESS; k++)
    {
        /* generator of F_k, coerced up to F_n */
        gr_ctx_struct * Fk = gr_tower_field_at(T, k);
        GR_TMP_INIT(t, Fk);
        status |= gr_gen(t, Fk);
        status |= gr_set_other(gr_vec_entry_ptr(vec, k - 1, ctx), t, Fk, top);
        GR_TMP_CLEAR(t, Fk);
    }

    return status;
}

static truth_t _is_zero(gr_srcptr x, gr_ctx_t ctx) { return gr_tower_is_zero_at(x, LEVEL(ctx), TOWER(ctx)); }
static truth_t _equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return gr_tower_equal_at(x, y, LEVEL(ctx), TOWER(ctx)); }

static truth_t
_is_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t res;
    GR_TMP_INIT(t, TOP(ctx));
    if (gr_sub_ui(t, x, 1, TOP(ctx)) == GR_SUCCESS)
        res = gr_tower_is_zero_at(t, LEVEL(ctx), TOWER(ctx));
    else
        res = T_UNKNOWN;
    GR_TMP_CLEAR(t, TOP(ctx));
    return res;
}

static truth_t
_is_neg_one(gr_srcptr x, gr_ctx_t ctx)
{
    gr_ptr t;
    truth_t res;
    GR_TMP_INIT(t, TOP(ctx));
    if (gr_add_ui(t, x, 1, TOP(ctx)) == GR_SUCCESS)
        res = gr_tower_is_zero_at(t, LEVEL(ctx), TOWER(ctx));
    else
        res = T_UNKNOWN;
    GR_TMP_CLEAR(t, TOP(ctx));
    return res;
}

static int _inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return gr_tower_inv_at(res, x, LEVEL(ctx), TOWER(ctx)); }
static int _div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return gr_tower_div_at(res, x, y, LEVEL(ctx), TOWER(ctx)); }

static truth_t
_is_invertible(gr_srcptr x, gr_ctx_t ctx)
{
    return truth_not(gr_tower_is_zero_at(x, LEVEL(ctx), TOWER(ctx)));
}

/* Division-free determinant (see lazy_dense.c). */
static int
_mat_det(gr_ptr res, const gr_mat_t A, gr_ctx_t ctx)
{
    if (gr_mat_nrows(A, ctx) <= 4)
        return gr_mat_det_cofactor(res, A, ctx);
    return gr_mat_det_berkowitz(res, A, ctx);
}

/* Factorization and roots of polynomials over the field (Trager's
   method; GR_UNABLE when the steps cannot be proven irreducible or the
   norms are too large). */
static int
_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, int flags, gr_ctx_t ctx)
{
    gr_ptr lc;
    int status;

    if (!VALID(ctx))
        return GR_UNABLE;

    GR_TMP_INIT(lc, TOP(ctx));
    status = gr_tower_poly_factor(lc, fac, mult, f, LEVEL(ctx), TOWER(ctx));
    if (status == GR_SUCCESS)
        status = gr_poly_set_scalar((gr_poly_struct *) c, lc, ctx);
    GR_TMP_CLEAR(lc, TOP(ctx));
    return status;
}

static int
_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f, int flags, gr_ctx_t ctx)
{
    if (!VALID(ctx))
        return GR_UNABLE;
    return gr_tower_poly_roots(roots, mult, f, LEVEL(ctx), TOWER(ctx));
}

int _gr_tower_field_methods_initialized = 0;

gr_static_method_table _gr_tower_field_methods;

gr_method_tab_input _gr_tower_field_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) _ctx_is_field},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _ctx_is_complex_vector_space},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) _ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) _ctx_is_finite_characteristic},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) _ctx_is_exact},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) _ctx_is_canonical},
    {GR_METHOD_CTX_NGENS,               (gr_funcptr) _ctx_ngens},
    {GR_METHOD_POLY_GCD,                (gr_funcptr) _gr_poly_gcd_subresultant},
    {GR_METHOD_POLY_XGCD,               (gr_funcptr) _gr_poly_xgcd_subresultant},
    {GR_METHOD_POLY_RESULTANT,          (gr_funcptr) _gr_poly_resultant_subresultant},
    {GR_METHOD_POLY_FACTOR,             (gr_funcptr) _poly_factor},
    {GR_METHOD_POLY_ROOTS,              (gr_funcptr) _poly_roots},
    {GR_METHOD_CTX_GEN_NAME,            (gr_funcptr) _ctx_gen_name},
    {GR_METHOD_CTX_BASE,                (gr_funcptr) _ctx_base},

    {GR_METHOD_INIT,            (gr_funcptr) _init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _write},
    {GR_METHOD_ZERO,            (gr_funcptr) _zero},
    {GR_METHOD_ONE,             (gr_funcptr) _one},
    {GR_METHOD_GENS,            (gr_funcptr) _gens},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _is_one},
    {GR_METHOD_IS_NEG_ONE,      (gr_funcptr) _is_neg_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _equal},
    {GR_METHOD_SET,             (gr_funcptr) _set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _set_other},
    {GR_METHOD_SET_STR,         (gr_funcptr) gr_generic_set_str_balance_additions},
    {GR_METHOD_NEG,             (gr_funcptr) _neg},
    {GR_METHOD_ADD,             (gr_funcptr) _add},
    {GR_METHOD_SUB,             (gr_funcptr) _sub},
    {GR_METHOD_MUL,             (gr_funcptr) _mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _mul_fmpz},
    {GR_METHOD_MUL_FMPQ,        (gr_funcptr) _mul_fmpq},
    {GR_METHOD_SQR,             (gr_funcptr) _sqr},
    {GR_METHOD_INV,             (gr_funcptr) _inv},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _is_invertible},
    {GR_METHOD_DIV,             (gr_funcptr) _div},
    {GR_METHOD_MAT_DET,         (gr_funcptr) _mat_det},
    {0,                         (gr_funcptr) NULL},
};

int
gr_tower_field_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)
{
    return gr_tower_get_acb(res, x, prec, TOWER(ctx));
}

void
gr_ctx_init_tower_field(gr_ctx_t ctx, gr_tower_t T)
{
    ctx->which_ring = GR_CTX_GR_TOWER_FIELD;
    ctx->sizeof_elem = gr_tower_field(T)->sizeof_elem;
    ctx->size_limit = WORD_MAX;
    TCTX(ctx)->T = T;
    TCTX(ctx)->level = T->length;
    TCTX(ctx)->num_trans = T->num_trans;
    TCTX(ctx)->structure_version = T->structure_version;
    TCTX(ctx)->top = gr_tower_field(T);
    FLINT_ASSERT(sizeof(tower_field_ctx_struct) <= GR_CTX_STRUCT_DATA_BYTES);

    ctx->methods = _gr_tower_field_methods;

    if (!_gr_tower_field_methods_initialized)
    {
        gr_method_tab_init(_gr_tower_field_methods, _gr_tower_field_methods_input);
        _gr_tower_field_methods_initialized = 1;
    }
}

POP_OPTIONS
