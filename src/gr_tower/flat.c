/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Flat representation of the top field of a tower: elements are
    multivariate rational functions (fmpz_mpoly_q) in the generators,
    with numerator and denominator reduced modulo the triangular set
    {m_1, ..., m_n} (with denominators cleared, so that m_k has leading
    term c_k x_k^{d_k} with c_k a positive integer). Since the leading
    monomials are pairwise coprime, the triangular set is a Groebner
    basis for the lex order with x_n > ... > x_1, and reduced numerators
    are canonical up to scaling.

    Compared with the nested representation, inversion is free
    (swap numerator and denominator), sparse expressions stay sparse,
    but denominators are not canonical and equality testing requires a
    cross-multiplication.
*/

#include <string.h>
#include "fmpq.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "acb.h"
#include "gr.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "gr_tower.h"
#include "gr_tower/impl.h"

PUSH_OPTIONS
OPTIMIZE_OSIZE

typedef struct
{
    gr_tower_flat_struct F;
    char ** vars;      /* (variable names for printing, rebuilt when the layout changes) */
    slong nvars;
}
flat_ctx_struct;

#define FLAT(ctx) ((flat_ctx_struct *) GR_CTX_DATA_AS_PTR(ctx))
#define FCORE(ctx) (&FLAT(ctx)->F)
#define MCTX(ctx) (FCORE(ctx)->mctx)
#define TOWER(ctx) (FCORE(ctx)->T)
#define VAR(ctx, k) GR_TOWER_FLAT_VAR(FCORE(ctx), k)

static int
_reduce(fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    return gr_tower_flat_reduce(x, FCORE(ctx));
}

/* -------------------------------------------------------------------- */
/* context                                                               */
/* -------------------------------------------------------------------- */

static int
_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong k;
    gr_tower_struct * T = TOWER(ctx);

    status |= gr_stream_write(out, "Tower field (flat) ");
    status |= gr_ctx_write(out, T->consts);
    for (k = 1; k <= T->num_trans; k++)
    {
        status |= gr_stream_write(out, (k == 1) ? "(" : ", ");
        status |= gr_stream_write(out, GR_TOWER_TRANS(T, k - 1)->name);
    }
    for (k = 1; k <= T->length; k++)
    {
        status |= gr_stream_write(out, (k == 1 && T->num_trans == 0) ? "(" : ", ");
        status |= gr_stream_write(out, GR_TOWER_STEP(T, k - 1)->name);
    }
    if (T->length + T->num_trans > 0)
        status |= gr_stream_write(out, ")");
    return status;
}

static void
_ctx_clear(gr_ctx_t ctx)
{
    flat_ctx_struct * F = FLAT(ctx);
    slong i;

    if (F->vars != NULL)
    {
        for (i = 0; i < F->nvars; i++)
            flint_free(F->vars[i]);
        flint_free(F->vars);
    }
    gr_tower_flat_clear(&F->F);
    flint_free(F);
}

static truth_t _ctx_is_field(gr_ctx_t ctx) { return gr_ctx_is_field(TOWER(ctx)->base); }
static truth_t _ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(TOWER(ctx)->base); }
static truth_t _ctx_is_rational_vector_space(gr_ctx_t ctx) { return gr_ctx_is_rational_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_real_vector_space(gr_ctx_t ctx) { return gr_ctx_is_real_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_complex_vector_space(gr_ctx_t ctx) { return gr_ctx_is_complex_vector_space(TOWER(ctx)->base); }
static truth_t _ctx_is_canonical(gr_ctx_t ctx) { return T_FALSE; }
static truth_t _ctx_is_threadsafe(gr_ctx_t ctx) { return T_FALSE; }
/* (as for the nested context: the base field K(t_1, ..., t_r)) */
static gr_ptr _ctx_base(gr_ctx_t ctx) { return TOWER(ctx)->base; }

static int
_ctx_ngens(slong * ngens, gr_ctx_t ctx)
{
    *ngens = TOWER(ctx)->num_trans + TOWER(ctx)->length;
    return GR_SUCCESS;
}

static int
_ctx_gen_name(char ** name, slong i, gr_ctx_t ctx)
{
    gr_tower_struct * T = TOWER(ctx);
    const char * s;
    size_t len;

    if (i < 0 || i >= T->num_trans + T->length)
        return GR_DOMAIN;

    s = (i < T->num_trans) ? GR_TOWER_TRANS(T, i)->name : GR_TOWER_STEP(T, i - T->num_trans)->name;
    len = strlen(s);
    *name = flint_malloc(len + 1);
    memcpy(*name, s, len + 1);
    return GR_SUCCESS;
}

/* -------------------------------------------------------------------- */
/* elements                                                              */
/* -------------------------------------------------------------------- */

static void _init(fmpz_mpoly_q_t x, gr_ctx_t ctx) { fmpz_mpoly_q_init(x, MCTX(ctx)); }
static void _clear(fmpz_mpoly_q_t x, gr_ctx_t ctx) { fmpz_mpoly_q_clear(x, MCTX(ctx)); }
static void _swap(fmpz_mpoly_q_t x, fmpz_mpoly_q_t y, gr_ctx_t ctx) { fmpz_mpoly_q_swap(x, y, MCTX(ctx)); }
static void _set_shallow(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_ctx_t ctx) { *res = *x; }

/* the variable names of the current layout of the polynomial context */
static void
_vars_update(flat_ctx_struct * F)
{
    gr_tower_struct * T = F->F.T;
    slong k;

    if (F->vars != NULL)
    {
        for (k = 0; k < F->nvars; k++)
            flint_free(F->vars[k]);
        flint_free(F->vars);
    }
    F->nvars = F->F.cap;
    F->vars = flint_malloc(sizeof(char *) * F->nvars);
    for (k = 0; k < F->nvars; k++)
    {
        F->vars[k] = flint_malloc(2);
        strcpy(F->vars[k], "x");
    }
    for (k = 0; k < T->num_gens; k++)
    {
        const char * name = GR_TOWER_GEN(T, k)->name;
        slong v = GR_TOWER_FLAT_VAR_D(&F->F, k);
        flint_free(F->vars[v]);
        F->vars[v] = flint_malloc(strlen(name) + 1);
        strcpy(F->vars[v], name);
    }
}

static int
_write(gr_stream_t out, const fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    flat_ctx_struct * F = FLAT(ctx);
    char * s;
    int status = GR_SUCCESS;

    _vars_update(F);

    if (fmpz_mpoly_is_one(fmpz_mpoly_q_denref(x), MCTX(ctx)))
    {
        s = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(x), (const char **) F->vars, MCTX(ctx));
        status |= gr_stream_write_free(out, s);
    }
    else
    {
        s = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_numref(x), (const char **) F->vars, MCTX(ctx));
        status |= gr_stream_write(out, "(");
        status |= gr_stream_write_free(out, s);
        status |= gr_stream_write(out, ")/(");
        s = fmpz_mpoly_get_str_pretty(fmpz_mpoly_q_denref(x), (const char **) F->vars, MCTX(ctx));
        status |= gr_stream_write_free(out, s);
        status |= gr_stream_write(out, ")");
    }

    return status;
}

static int _zero(fmpz_mpoly_q_t res, gr_ctx_t ctx) { fmpz_mpoly_q_zero(res, MCTX(ctx)); return GR_SUCCESS; }
static int _one(fmpz_mpoly_q_t res, gr_ctx_t ctx) { fmpz_mpoly_q_one(res, MCTX(ctx)); return GR_SUCCESS; }
static int _set(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_ctx_t ctx) { fmpz_mpoly_q_set(res, x, MCTX(ctx)); return GR_SUCCESS; }
static int _set_si(fmpz_mpoly_q_t res, slong c, gr_ctx_t ctx) { fmpz_mpoly_q_set_si(res, c, MCTX(ctx)); return GR_SUCCESS; }
static int _set_ui(fmpz_mpoly_q_t res, ulong c, gr_ctx_t ctx) { fmpz_t t; fmpz_init_set_ui(t, c); fmpz_mpoly_q_set_fmpz(res, t, MCTX(ctx)); fmpz_clear(t); return GR_SUCCESS; }
static int _set_fmpz(fmpz_mpoly_q_t res, const fmpz_t c, gr_ctx_t ctx) { fmpz_mpoly_q_set_fmpz(res, c, MCTX(ctx)); return GR_SUCCESS; }
static int _set_fmpq(fmpz_mpoly_q_t res, const fmpq_t c, gr_ctx_t ctx) { fmpz_mpoly_q_set_fmpq(res, c, MCTX(ctx)); return GR_SUCCESS; }

static int
_gen_i(fmpz_mpoly_q_t res, slong i, gr_ctx_t ctx)
{
    slong r = TOWER(ctx)->num_trans;

    /* generators are ordered t_1, ..., t_r, a_1, ..., a_n */
    if (i < 0 || i >= r + TOWER(ctx)->length)
        return GR_DOMAIN;
    if (i < r)
        fmpz_mpoly_gen(fmpz_mpoly_q_numref(res), GR_TOWER_FLAT_TVAR(FCORE(ctx), i + 1), MCTX(ctx));
    else
        fmpz_mpoly_gen(fmpz_mpoly_q_numref(res), VAR(ctx, i - r + 1), MCTX(ctx));
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), MCTX(ctx));
    return _reduce(res, ctx);
}

static int
_gens(gr_vec_t vec, gr_ctx_t ctx)
{
    slong i, n = TOWER(ctx)->length + TOWER(ctx)->num_trans;
    int status = GR_SUCCESS;
    gr_vec_set_length(vec, n, ctx);
    for (i = 0; i < n; i++)
        status |= _gen_i(gr_vec_entry_ptr(vec, i, ctx), i, ctx);
    return status;
}

/* nested top field element -> flat */
int
gr_tower_flat_set_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_tower_flat_set_nested_at(res, x, TOWER(ctx)->length, FCORE(ctx));
}

int
gr_tower_flat_get_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    return gr_tower_flat_get_nested_at(res, x, TOWER(ctx)->length, FCORE(ctx));
}

static int
_set_other(fmpz_mpoly_q_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
        return _set(res, x, ctx);

    if (x_ctx == gr_tower_field(TOWER(ctx)) ||
        (x_ctx->which_ring == GR_CTX_GR_TOWER_FIELD && GR_CTX_DATA_AS_PTR(x_ctx) == TOWER(ctx)))
        return gr_tower_flat_set_nested(res, x, ctx);

    {
        fmpq_t c;
        int status;
        fmpq_init(c);
        status = gr_get_fmpq(c, x, x_ctx);
        if (status == GR_SUCCESS)
            fmpz_mpoly_q_set_fmpq(res, c, MCTX(ctx));
        fmpq_clear(c);
        if (status == GR_SUCCESS)
            return status;
    }

    return GR_UNABLE;
}

static int
_randtest(fmpz_mpoly_q_t res, flint_rand_t state, gr_ctx_t ctx)
{
    /* random polynomial in the actual generators (unused slots must not appear) */
    gr_tower_struct * T = TOWER(ctx);
    slong i, k, nterms = 1 + n_randint(state, 1 + n_randint(state, 4));
    ulong * exp = flint_calloc(FLAT(ctx)->nvars, sizeof(ulong));
    fmpz_t c;

    /* sparse and of low degree: generic test code takes powers of these
       elements, and the number of terms grows quickly with the number
       of generators */
    fmpz_init(c);
    fmpz_mpoly_zero(fmpz_mpoly_q_numref(res), MCTX(ctx));
    for (i = 0; i < nterms; i++)
    {
        for (k = 1; k <= T->length; k++)
            exp[VAR(ctx, k)] = n_randint(state, 2) ? n_randint(state, 3) : 0;
        for (k = 1; k <= T->num_trans; k++)
            exp[GR_TOWER_FLAT_TVAR(FCORE(ctx), k)] = n_randint(state, 2) ? n_randint(state, 3) : 0;
        fmpz_randtest(c, state, 8);
        fmpz_mpoly_set_coeff_fmpz_ui(fmpz_mpoly_q_numref(res), c, exp, MCTX(ctx));
    }
    fmpz_mpoly_one(fmpz_mpoly_q_denref(res), MCTX(ctx));
    fmpz_clear(c);
    flint_free(exp);
    return _reduce(res, ctx);
}

int
gr_tower_field_flat_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx)
{
    return gr_tower_flat_get_acb(res, x, prec, FCORE(ctx));
}

/* The zero test of the numerator: the elements of this context are
   polynomials in a fixed context, so the tower and its flat layout
   must not change under them; the parts of the test which may change
   them (relations found between transcendental generators, which the
   lazy field records in place) run on a copy of the tower. */
static truth_t
_num_is_zero(fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    return _gr_tower_flat_num_is_zero_fixed(x, FCORE(ctx));
}

static truth_t
_is_zero(const fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    truth_t res;
    fmpz_mpoly_q_init(t, MCTX(ctx));
    fmpz_mpoly_q_set(t, x, MCTX(ctx));
    res = _num_is_zero(t, ctx);
    fmpz_mpoly_q_clear(t, MCTX(ctx));
    return res;
}

static int _neg(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_ctx_t ctx) { fmpz_mpoly_q_neg(res, x, MCTX(ctx)); return GR_SUCCESS; }
static int _add(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_ctx_t ctx) { fmpz_mpoly_q_add(res, x, y, MCTX(ctx)); return _reduce(res, ctx); }
static int _sub(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_ctx_t ctx) { fmpz_mpoly_q_sub(res, x, y, MCTX(ctx)); return _reduce(res, ctx); }
static int _mul(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_ctx_t ctx) { fmpz_mpoly_q_mul(res, x, y, MCTX(ctx)); return _reduce(res, ctx); }

static truth_t
_equal(const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    truth_t res;

    if (fmpz_mpoly_q_equal(x, y, MCTX(ctx)))
        return T_TRUE;

    fmpz_mpoly_q_init(t, MCTX(ctx));
    fmpz_mpoly_q_sub(t, x, y, MCTX(ctx));
    res = _num_is_zero(t, ctx);
    fmpz_mpoly_q_clear(t, MCTX(ctx));
    return res;
}

static truth_t
_is_one(const fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    fmpz_mpoly_q_t t;
    truth_t res;
    fmpz_t one;

    if (fmpz_mpoly_q_is_one(x, MCTX(ctx)))
        return T_TRUE;

    fmpz_init_set_ui(one, 1);
    fmpz_mpoly_q_init(t, MCTX(ctx));
    fmpz_mpoly_q_sub_fmpz(t, x, one, MCTX(ctx));
    res = _num_is_zero(t, ctx);
    fmpz_mpoly_q_clear(t, MCTX(ctx));
    fmpz_clear(one);
    return res;
}

static int
_inv(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, gr_ctx_t ctx)
{
    truth_t z = _is_zero(x, ctx);

    if (z == T_TRUE)
        return GR_DOMAIN;
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    fmpz_mpoly_q_inv(res, x, MCTX(ctx));
    return GR_SUCCESS;
}

static int
_div(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_q_t y, gr_ctx_t ctx)
{
    truth_t z = _is_zero(y, ctx);

    if (z == T_TRUE)
        return GR_DOMAIN;
    if (z == T_UNKNOWN)
        return GR_UNABLE;

    fmpz_mpoly_q_div(res, x, y, MCTX(ctx));
    return _reduce(res, ctx);
}

/* Field LU is fine here (division is a swap), but Berkowitz avoids
   denominator growth; keep the generic choice for now. */

int _gr_tower_flat_methods_initialized = 0;

gr_static_method_table _gr_tower_flat_methods;

gr_method_tab_input _gr_tower_flat_methods_input[] =
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
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) _ctx_is_canonical},
    {GR_METHOD_CTX_NGENS,               (gr_funcptr) _ctx_ngens},
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
    {GR_METHOD_EQUAL,           (gr_funcptr) _equal},
    {GR_METHOD_SET,             (gr_funcptr) _set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) _set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) _set_other},
    {GR_METHOD_NEG,             (gr_funcptr) _neg},
    {GR_METHOD_ADD,             (gr_funcptr) _add},
    {GR_METHOD_SUB,             (gr_funcptr) _sub},
    {GR_METHOD_MUL,             (gr_funcptr) _mul},
    {GR_METHOD_INV,             (gr_funcptr) _inv},
    {GR_METHOD_DIV,             (gr_funcptr) _div},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_tower_field_flat(gr_ctx_t ctx, gr_tower_t T)
{
    flat_ctx_struct * F;

    F = flint_calloc(1, sizeof(flat_ctx_struct));
    gr_tower_flat_init(&F->F, T, T->num_gens);
    _vars_update(F);

    ctx->which_ring = GR_CTX_GR_TOWER_FIELD_FLAT;
    ctx->sizeof_elem = sizeof(fmpz_mpoly_q_struct);
    ctx->size_limit = WORD_MAX;
    GR_CTX_DATA_AS_PTR(ctx) = F;

    ctx->methods = _gr_tower_flat_methods;

    if (!_gr_tower_flat_methods_initialized)
    {
        gr_method_tab_init(_gr_tower_flat_methods, _gr_tower_flat_methods_input);
        _gr_tower_flat_methods_initialized = 1;
    }
}

POP_OPTIONS
