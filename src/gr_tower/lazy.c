/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/* Lazy fields: contexts, elements, the locking layer and the method table (see lazy_impl.h for the layout of the implementation). */

/* (for recursive mutexes: pthread_mutexattr_settype) */
#define _GNU_SOURCE

#include "gr_tower/lazy_impl.h"

/* the depth of the lazy field operations in progress in this thread
   (dense forms are made and used only outside of them: the
   implementations read the flat data of their elements) */
FLINT_TLS_PREFIX slong _gr_tower_lazy_tls_depth = 0;

gr_tower_flat_struct *
_gr_tower_lazy_new_tower(gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_struct * T;
    gr_tower_flat_struct * F;

    if (L->num_towers == L->alloc_towers)
    {
        L->alloc_towers = FLINT_MAX(4, 2 * L->alloc_towers);
        L->towers = flint_realloc(L->towers, L->alloc_towers * sizeof(gr_tower_flat_struct *));
    }

    T = flint_malloc(sizeof(gr_tower_struct));
    gr_tower_init(T, L->base);
    T->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
    T->options = L->options;
    F = &T->flat;
    L->towers[L->num_towers++] = F;
    return F;
}

void
_gr_tower_lazy_ref(gr_tower_flat_struct * F)
{
    /* (pinned towers are never collected: not counted, so that the
       rational elements of the trivial tower need no lock) */
    if (!(F->gc & GR_TOWER_GC_PINNED))
        F->refs++;
}

void
_gr_tower_lazy_unref(gr_tower_lazy_ctx_struct * L, gr_tower_flat_struct * F)
{
    if (F->gc & GR_TOWER_GC_PINNED)
        return;
    F->refs--;
    if (F->refs == 0 && !(F->gc & GR_TOWER_GC_CANDIDATE))
    {
        if (L->gc_len == L->gc_alloc)
        {
            L->gc_alloc = FLINT_MAX(16, 2 * L->gc_alloc);
            L->gc_queue = flint_realloc(L->gc_queue, L->gc_alloc * sizeof(gr_tower_flat_struct *));
        }
        F->gc |= GR_TOWER_GC_CANDIDATE;
        L->gc_queue[L->gc_len++] = F;
    }
}

/* a tower of T's context containing the definition id (NULL if none) */
static gr_tower_flat_struct *
_gr_tower_lazy_find_live_def(gr_tower_lazy_ctx_struct * L, ulong def_id, slong * gid)
{
    slong i;
    for (i = L->num_towers - 1; i >= 0; i--)
    {
        gr_tower_flat_struct * F = L->towers[i];
        slong d;
        if (F->gc & GR_TOWER_GC_DEAD)
            continue;
        d = gr_tower_find_def_order(F->T, def_id);
        if (d >= 0)
        {
            *gid = GR_TOWER_GEN(F->T, d)->gid;
            return F;
        }
    }
    return NULL;
}

static void
_gr_tower_lazy_gc(gr_tower_lazy_ctx_struct * L)
{
    slong i, j, ndead = 0;

    for (i = 0; i < L->gc_len; i++)
    {
        gr_tower_flat_struct * F = L->gc_queue[i];
        F->gc &= ~GR_TOWER_GC_CANDIDATE;
        if (F->refs == 0 && !(F->gc & GR_TOWER_GC_PINNED) && F != L->trivial)
        {
            F->gc |= GR_TOWER_GC_DEAD;
            ndead++;
        }
    }
    L->gc_len = 0;

    if (ndead == 0)
        return;

    /* algebraic number cache: move or drop the entries */
    for (i = j = 0; i < L->num_qqbars; i++)
    {
        gr_tower_lazy_qqbar_entry_struct * e = L->qqbars + i;
        int keep = 1;
        if (e->F->gc & GR_TOWER_GC_DEAD)
        {
            slong d = gr_tower_gid_order(e->F->T, e->gid), gid;
            ulong id = (d >= 0) ? GR_TOWER_GEN(e->F->T, d)->def_id : 0;
            gr_tower_flat_struct * G = (id != 0) ? _gr_tower_lazy_find_live_def(L, id, &gid) : NULL;
            if (G != NULL)
                e->F = G, e->gid = gid;
            else
                keep = 0;
        }
        if (keep)
            L->qqbars[j++] = *e;
        else
            qqbar_clear(&e->x);
    }
    L->num_qqbars = j;

    /* aliases (expressions of definitions in the generators of a tower) */
    for (i = j = 0; i < L->num_aliases; i++)
    {
        gr_tower_lazy_alias_struct * a = L->aliases + i;
        if (a->F->gc & GR_TOWER_GC_DEAD)
            fmpz_mpoly_q_clear(&a->data, a->mctx);
        else
            L->aliases[j++] = *a;
    }
    L->num_aliases = j;

    for (i = j = 0; i < L->num_towers; i++)
    {
        gr_tower_flat_struct * F = L->towers[i];
        if (F->gc & GR_TOWER_GC_DEAD)
        {
            gr_tower_struct * T = F->T;
            gr_tower_clear(T);
            flint_free(T);
            L->gc_collected++;
        }
        else
            L->towers[j++] = F;
    }
    L->num_towers = j;
}

/* -------------------------------------------------------------------- */
/* context                                                               */
/* -------------------------------------------------------------------- */

/*
    Assigns a definition id to a generator, and a name unique in the
    context (a1, a2, ... for algebraic definitions, t1, t2, ... for
    transcendental ones, in order of creation; pi keeps its name), which
    follows the generator into every tower containing it.
*/
void
_gr_tower_lazy_new_def(gr_tower_gen_struct * g, gr_tower_t T, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    char buf[32];

    g->def_id = ++L->next_def_id;

    if (g->def_kind == GR_TOWER_PI)
        return;
    if (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param == 4)
    {
        _gr_tower_gen_rename(T, g, "i");
        return;
    }
    if (g->def_kind == GR_TOWER_ALGEBRAIC || g->def_kind == GR_TOWER_ROOT || g->def_kind == GR_TOWER_ROOT_OF_UNITY ||
        g->def_kind == GR_TOWER_TAN_PI)
        flint_sprintf(buf, "a%wd", ++L->num_alg_names);
    else
        flint_sprintf(buf, "t%wd", ++L->num_trans_names);
    _gr_tower_gen_rename(T, g, buf);
}

static int
_gr_tower_lazy_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int f = VIEW(ctx)->field_flags;
    return gr_stream_write(out,
        (f & GR_TOWER_LAZY_REAL) ?
            ((f & GR_TOWER_LAZY_ALGEBRAIC) ? "Real algebraic field (lazy towers)" : "Real field (lazy towers)") :
            ((f & GR_TOWER_LAZY_ALGEBRAIC) ? "Complex algebraic field (lazy towers)" : "Complex field (lazy towers)"));
}

void
gr_tower_lazy_ctx_set_gen_flags(gr_ctx_t ctx, int flags)
{
    /* (in terms of the options; a cyclotomic degree limit already set
       is kept when composite roots are requested) */
    slong * opt;
    _gr_tower_lazy_lock(ctx);
    opt = LAZY(ctx)->options;
    opt[GR_TOWER_OPT_SPLIT_IMAGINARY] = (flags & GR_TOWER_GENS_SPLIT_IMAGINARY) != 0;
    opt[GR_TOWER_OPT_COMPOSITE_RADICALS] = (flags & GR_TOWER_GENS_COMPOSITE_RADICALS) != 0;
    if (flags & GR_TOWER_GENS_COMPOSITE_ROOTS)
    {
        if (opt[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT] <= 0)
            opt[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT] = GR_TOWER_COMPOSITE_ROOTS_DEGREE;
    }
    else
        opt[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT] = 0;
    _gr_tower_lazy_unlock(ctx);
}

int
gr_tower_lazy_ctx_gen_flags(gr_ctx_t ctx)
{
    const slong * opt = LAZY(ctx)->options;
    return (opt[GR_TOWER_OPT_SPLIT_IMAGINARY] ? GR_TOWER_GENS_SPLIT_IMAGINARY : 0) |
           (opt[GR_TOWER_OPT_COMPOSITE_RADICALS] ? GR_TOWER_GENS_COMPOSITE_RADICALS : 0) |
           (opt[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT] > 0 ? GR_TOWER_GENS_COMPOSITE_ROOTS : 0);
}

int
gr_tower_lazy_ctx_set_option(gr_ctx_t ctx, slong option, slong value)
{
    if (option < 0 || option >= GR_TOWER_OPT_NUM_OPTIONS)
        flint_throw(FLINT_ERROR, "(%s): invalid option %wd\n", __func__, option);
    if (option == GR_TOWER_OPT_PRINT_DIGITS && value <= 0)
        value = GR_TOWER_PRINT_DIGITS_DEFAULT;
    if (!gr_tower_option_valid(option, value))
        return GR_DOMAIN;
    _gr_tower_lazy_lock(ctx);
    LAZY(ctx)->options[option] = value;
    _gr_tower_lazy_unlock(ctx);
    return GR_SUCCESS;
}

slong
gr_tower_lazy_ctx_get_option(gr_ctx_t ctx, slong option)
{
    if (option < 0 || option >= GR_TOWER_OPT_NUM_OPTIONS)
        flint_throw(FLINT_ERROR, "(%s): invalid option %wd\n", __func__, option);
    return LAZY(ctx)->options[option];
}

void
gr_tower_lazy_ctx_set_print(gr_ctx_t ctx, int flags, slong digits)
{
    LAZY(ctx)->options[GR_TOWER_OPT_PRINT_FLAGS] = flags;
    LAZY(ctx)->options[GR_TOWER_OPT_PRINT_DIGITS] = (digits <= 0) ? GR_TOWER_PRINT_DIGITS_DEFAULT : digits;
}

int gr_tower_lazy_ctx_print_flags(gr_ctx_t ctx) { return LAZY(ctx)->options[GR_TOWER_OPT_PRINT_FLAGS]; }
int gr_tower_lazy_ctx_field_flags(gr_ctx_t ctx) { return VIEW(ctx)->field_flags; }

/* the number of towers of the state (after a collection of unused ones) */
slong
gr_tower_lazy_ctx_num_towers(gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong n;
    _gr_tower_lazy_lock(ctx);
    if (L->gc_len > 0 && L->depth == 1)
        _gr_tower_lazy_gc(L);
    n = L->num_towers;
    _gr_tower_lazy_unlock(ctx);
    return n;
}

static void
_gr_tower_lazy_ctx_clear(gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong i;

    if (--L->refcount > 0)
        return;

    for (i = 0; i < L->num_aliases; i++)
        fmpz_mpoly_q_clear(&L->aliases[i].data, L->aliases[i].mctx);
    flint_free(L->aliases);

    for (i = 0; i < L->num_qqbars; i++)
        qqbar_clear(&L->qqbars[i].x);
    flint_free(L->qqbars);

    for (i = 0; i < L->num_roots; i++)
        fmpz_clear(&L->roots[i].p);
    flint_free(L->roots);

    for (i = 0; i < L->num_const_trans; i++)
    {
        fmpq_clear(&L->const_trans[i].x);
        fmpq_clear(&L->const_trans[i].y);
    }
    flint_free(L->const_trans);

    for (i = 0; i < L->num_towers; i++)
    {
        gr_tower_struct * T = L->towers[i]->T;
        gr_tower_clear(T);
        flint_free(T);
    }
    flint_free(L->towers);
    flint_free(L->gc_queue);

    #if FLINT_USES_PTHREAD
    pthread_mutex_destroy(&L->mutex);
#endif
    flint_free(L);
}

static truth_t _gr_tower_lazy_ctx_is_field(gr_ctx_t ctx) { return gr_ctx_is_field(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_finite_characteristic(gr_ctx_t ctx) { return gr_ctx_is_finite_characteristic(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_rational_vector_space(gr_ctx_t ctx) { return gr_ctx_is_rational_vector_space(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_real_vector_space(gr_ctx_t ctx) { return gr_ctx_is_real_vector_space(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_complex_vector_space(gr_ctx_t ctx) { return gr_ctx_is_complex_vector_space(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_exact(gr_ctx_t ctx) { return gr_ctx_is_exact(LAZY(ctx)->base); }
static truth_t _gr_tower_lazy_ctx_is_canonical(gr_ctx_t ctx) { return T_FALSE; }
#if FLINT_USES_PTHREAD
static truth_t _gr_tower_lazy_ctx_is_threadsafe(gr_ctx_t ctx) { return T_TRUE; }
#else
static truth_t _gr_tower_lazy_ctx_is_threadsafe(gr_ctx_t ctx) { return T_FALSE; }
#endif

void
_gr_tower_lazy_lock(gr_ctx_t ctx)
{
#if FLINT_USES_PTHREAD
    pthread_mutex_lock(&LAZY(ctx)->mutex);
#endif
    LAZY(ctx)->depth++;
    _gr_tower_lazy_tls_depth++;
}

static void _gr_tower_lazy_scratch_release(void);

void
_gr_tower_lazy_unlock(gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    if (L->depth == 1 && L->gc_len >= LAZY_GC_BATCH)
        _gr_tower_lazy_gc(L);
    LAZY(ctx)->depth--;
    _gr_tower_lazy_tls_depth--;
    if (_gr_tower_lazy_tls_depth == 0)
        _gr_tower_lazy_scratch_release();
#if FLINT_USES_PTHREAD
    pthread_mutex_unlock(&LAZY(ctx)->mutex);
#endif
}

/* whether the operation in progress is the outermost one (the
   restrictions of the real and algebraic subfields apply to it, not to
   the intermediate values of the implementations) */
int
_gr_tower_lazy_outermost(gr_ctx_t ctx)
{
    return LAZY(ctx)->depth <= 1;
}

static gr_ptr _gr_tower_lazy_ctx_base(gr_ctx_t ctx) { return LAZY(ctx)->base; }

/* -------------------------------------------------------------------- */
/* elements                                                              */
/* -------------------------------------------------------------------- */

void
_gr_tower_lazy_init(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    x->F = LAZY(ctx)->trivial;
    _gr_tower_lazy_ref(x->F);
    x->level = 0;
    x->structure_version = x->F->T->structure_version;
    x->reduced_version = 0;
    x->shallow = 0;
    x->repr = LAZY_REPR_RATIONAL;
    fmpq_init(&x->elem.q);
    LAZY_POISON_FLAT(x);
    LAZY_POISON_DENSE(x);
}

/* clears the data of x (not the tower reference) */
void
_gr_tower_lazy_clear_data(gr_tower_lazy_elem_t x)
{
    if (x->repr == LAZY_REPR_RATIONAL)
        fmpq_clear(&x->elem.q);
    else if (x->repr == LAZY_REPR_FLAT)
    {
        if (x->shallow != LAZY_SHALLOW_STALE)
            fmpz_mpoly_q_clear(&x->elem.flat.data, x->elem.flat.mctx);
    }
    else
        fmpq_poly_clear(&x->elem.dense.poly);
}

void
_gr_tower_lazy_clear(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    _gr_tower_lazy_clear_data(x);
    _gr_tower_lazy_unref(LAZY(ctx), x->F);
}

void
_gr_tower_lazy_swap(gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct t = *x;
    *x = *y;
    *y = t;
}

static void
_gr_tower_lazy_set_shallow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    *res = *x;
    if (res->shallow == LAZY_SHALLOW_NONE)
        res->shallow = LAZY_SHALLOW_COPY;
}

/* Replaces the content of res by a zero element of (F, level): rational
   for the trivial tower, flat otherwise. */
void
_gr_tower_lazy_fresh(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong level, gr_ctx_t ctx)
{
    gr_tower_flat_ensure(F);
    /* (res is an output, so it owns its data, also when it is a shallow
       copy: the gr idioms move elements with set_shallow, as in
       gr_mat_nonsingular_solve_tril_classical) */
    _gr_tower_lazy_clear_data(res);
    res->shallow = LAZY_SHALLOW_NONE;
    _gr_tower_lazy_ref(F);
    _gr_tower_lazy_unref(LAZY(ctx), res->F);
    res->F = F;
    res->level = level;
    res->structure_version = F->T->structure_version;
    res->reduced_version = 0;
    if (F == LAZY(ctx)->trivial)
    {
        res->repr = LAZY_REPR_RATIONAL;
        fmpq_init(&res->elem.q);
        LAZY_POISON_FLAT(res);
        LAZY_POISON_DENSE(res);
    }
    else
    {
        res->repr = LAZY_REPR_FLAT;
        res->elem.flat.mctx = F->mctx;
        fmpz_mpoly_q_init(&res->elem.flat.data, res->elem.flat.mctx);
        LAZY_POISON_DENSE(res);
    }
}

/* Makes the representation of x (which must not be read concurrently:
   a result, or a private temporary) flat, in the current context of
   its tower. */
static void
_gr_tower_lazy_make_flat(gr_tower_lazy_elem_t x)
{
    gr_tower_flat_struct * F = x->F;

    if (x->repr == LAZY_REPR_FLAT)
        return;

    gr_tower_flat_ensure(F);

    if (x->repr == LAZY_REPR_RATIONAL)
    {
        fmpq_t q;
        fmpq_init(q);
        fmpq_swap(q, &x->elem.q);
        fmpq_clear(&x->elem.q);
        LAZY_POISON_RATIONAL(x);
        x->repr = LAZY_REPR_FLAT;
        x->elem.flat.mctx = F->mctx;
        fmpz_mpoly_q_init(&x->elem.flat.data, x->elem.flat.mctx);
        fmpz_mpoly_q_set_fmpq(&x->elem.flat.data, q, x->elem.flat.mctx);
        fmpq_clear(q);
        x->level = 0;
    }
    else
    {
        fmpz_mpoly_q_t t;
        fmpz_mpoly_q_init(t, F->mctx);
        _gr_tower_lazy_dense_get_flat(t, x, F);
        if (!x->shallow)
            fmpq_poly_clear(&x->elem.dense.poly);
        x->repr = LAZY_REPR_FLAT;
        x->elem.flat.mctx = F->mctx;
        fmpz_mpoly_q_init(&x->elem.flat.data, x->elem.flat.mctx);
        fmpz_mpoly_q_swap(&x->elem.flat.data, t, x->elem.flat.mctx);
        fmpz_mpoly_q_clear(t, x->elem.flat.mctx);
        LAZY_POISON_DENSE(x);
        x->shallow = 0;
        x->reduced_version = F->ideal_version;
        _gr_tower_lazy_shrink(x);
    }
    x->structure_version = F->T->structure_version;
}

void
_gr_tower_lazy_update(gr_tower_lazy_elem_t x)
{
    _gr_tower_lazy_make_flat(x);

    gr_tower_flat_ensure(x->F);

    if (x->elem.flat.mctx != x->F->mctx)
    {
        if (x->shallow)
        {
            /* the data belongs to another element (which stays valid in
               its old context): the converted copy is owned by F */
            fmpz_mpoly_q_struct * t = _gr_tower_flat_stale_alloc(x->F);
            gr_tower_flat_convert(t, &x->elem.flat.data, x->elem.flat.mctx, x->F);
            x->elem.flat.mctx = x->F->mctx;
            x->elem.flat.data = *t;
            x->shallow = LAZY_SHALLOW_STALE;
        }
        else
        {
            fmpz_mpoly_q_t t;
            fmpz_mpoly_q_init(t, x->F->mctx);
            gr_tower_flat_convert(t, &x->elem.flat.data, x->elem.flat.mctx, x->F);
            fmpz_mpoly_q_clear(&x->elem.flat.data, x->elem.flat.mctx);
            x->elem.flat.mctx = x->F->mctx;
            fmpz_mpoly_q_init(&x->elem.flat.data, x->elem.flat.mctx);
            fmpz_mpoly_q_swap(&x->elem.flat.data, t, x->elem.flat.mctx);
            fmpz_mpoly_q_clear(t, x->elem.flat.mctx);
        }
    }

    /* the definition order may have changed (the level is a prefix length) */
    if (x->structure_version != x->F->T->structure_version)
        _gr_tower_lazy_shrink(x);
}

/*
    Temporary flat copies of operands in the rational or the dense
    representation (which are never changed in place while they may be
    read without the lock). They live in a per-thread pool and are
    released when the outermost locked operation of the thread returns.
*/
typedef struct
{
    gr_tower_lazy_elem_struct ** elems;    /* (allocated one by one: their addresses are stable) */
    slong num;                    /* in use */
    slong alloc;                  /* allocated (the first alloc entries of elems are valid) */
    slong cap;                    /* capacity of the elems array */
}

lazy_scratch_struct;

static FLINT_TLS_PREFIX lazy_scratch_struct _gr_tower_lazy_tls_scratch = { NULL, 0, 0, 0 };
static FLINT_TLS_PREFIX int _gr_tower_lazy_tls_scratch_registered = 0;

/* (registered with flint_register_cleanup_function, so that the pool of
   a thread is freed by flint_cleanup in that thread) */
static void
_gr_tower_lazy_scratch_cleanup(void)
{
    lazy_scratch_struct * S = &_gr_tower_lazy_tls_scratch;
    slong i;
    for (i = 0; i < S->num; i++)
        _gr_tower_lazy_clear_data(S->elems[i]);
    for (i = 0; i < S->alloc; i++)
        flint_free(S->elems[i]);
    flint_free(S->elems);
    S->elems = NULL;
    S->num = S->alloc = S->cap = 0;
    _gr_tower_lazy_tls_scratch_registered = 0;
}

static gr_tower_lazy_elem_struct *
_gr_tower_lazy_scratch_alloc(void)
{
    lazy_scratch_struct * S = &_gr_tower_lazy_tls_scratch;
    if (!_gr_tower_lazy_tls_scratch_registered)
    {
        flint_register_cleanup_function(_gr_tower_lazy_scratch_cleanup);
        _gr_tower_lazy_tls_scratch_registered = 1;
    }
    if (S->num == S->cap)
    {
        S->cap = FLINT_MAX(8, 2 * S->cap);
        S->elems = flint_realloc(S->elems, S->cap * sizeof(gr_tower_lazy_elem_struct *));
    }
    if (S->num == S->alloc)
    {
        S->elems[S->num] = flint_malloc(sizeof(gr_tower_lazy_elem_struct));
        S->alloc++;
    }
    return S->elems[S->num++];
}

/* (called when the outermost locked operation returns) */
static void
_gr_tower_lazy_scratch_release(void)
{
    lazy_scratch_struct * S = &_gr_tower_lazy_tls_scratch;
    slong i;
    for (i = 0; i < S->num; i++)
        _gr_tower_lazy_clear_data(S->elems[i]);
    S->num = 0;
    if (S->alloc > 64)
    {
        for (i = 0; i < S->alloc; i++)
            flint_free(S->elems[i]);
        flint_free(S->elems);
        S->elems = NULL;
        S->alloc = S->cap = 0;
    }
}

/* The flat view of an operand (under the lock): x itself, brought up to
   date, when it is flat; otherwise a temporary flat copy, valid until
   the outermost locked operation returns. The view is not to be written
   to, and not to be kept beyond the operation. */
gr_tower_lazy_elem_struct *
_gr_tower_lazy_flat_view(const gr_tower_lazy_elem_t x_in)
{
    gr_tower_lazy_elem_struct * x = (gr_tower_lazy_elem_struct *) x_in;
    gr_tower_lazy_elem_struct * v;

    if (x->repr == LAZY_REPR_FLAT)
    {
        _gr_tower_lazy_update(x);
        return x;
    }

    /* (a shallow struct: the pool entry does not own a tower reference
       and is cleared by _gr_tower_lazy_clear_data) */
    v = _gr_tower_lazy_scratch_alloc();
    v->F = x->F;
    v->level = x->level;
    v->structure_version = x->structure_version;
    v->reduced_version = x->reduced_version;
    v->shallow = 1;
    gr_tower_flat_ensure(v->F);
    v->repr = LAZY_REPR_FLAT;
    v->elem.flat.mctx = v->F->mctx;
    fmpz_mpoly_q_init(&v->elem.flat.data, v->elem.flat.mctx);
    if (x->repr == LAZY_REPR_RATIONAL)
    {
        fmpz_mpoly_q_set_fmpq(&v->elem.flat.data, &x->elem.q, v->elem.flat.mctx);
        v->level = 0;
        v->reduced_version = 0;
    }
    else
    {
        _gr_tower_lazy_dense_get_flat(&v->elem.flat.data, x, v->F);
        v->reduced_version = v->F->ideal_version;
    }
    v->structure_version = v->F->T->structure_version;
    _gr_tower_lazy_shrink(v);
    return v;
}

/* Lowers the level to the shortest prefix actually needed. */
void
_gr_tower_lazy_shrink(gr_tower_lazy_elem_t x)
{
    x->level = gr_tower_flat_level(&x->elem.flat.data, x->F);
    x->structure_version = x->F->T->structure_version;
}

/* Algebraic level (number of algebraic generators) of x's prefix. */
slong
_gr_tower_lazy_alg_level(const gr_tower_lazy_elem_t x)
{
    return gr_tower_prefix_length(x->F->T, x->level);
}

/* res = x, in the representation of x (a dense form is kept) */
int
_gr_tower_lazy_set(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (res == x)
        return GR_SUCCESS;

    if (x->repr == LAZY_REPR_RATIONAL)
    {
        _gr_tower_lazy_fresh(res, x->F, 0, ctx);
        fmpq_set(&res->elem.q, &x->elem.q);
    }
    else if (x->repr == LAZY_REPR_DENSE)
    {
        _gr_tower_lazy_fresh(res, x->F, x->level, ctx);
        fmpz_mpoly_q_clear(&res->elem.flat.data, res->elem.flat.mctx);
        LAZY_POISON_FLAT(res);
        res->repr = LAZY_REPR_DENSE;
        fmpq_poly_init(&res->elem.dense.poly);
        fmpq_poly_set(&res->elem.dense.poly, &x->elem.dense.poly);
        res->elem.dense.nf = x->elem.dense.nf;
        res->structure_version = x->structure_version;
        res->reduced_version = x->reduced_version;
    }
    else
    {
        _gr_tower_lazy_update((gr_tower_lazy_elem_struct *) x);
        if (x->F == LAZY(ctx)->trivial)
        {
            /* (a flat view of a rational element) */
            _gr_tower_lazy_fresh(res, x->F, 0, ctx);
            (void) fmpz_mpoly_q_get_fmpq(&res->elem.q, &x->elem.flat.data, x->elem.flat.mctx);
            return GR_SUCCESS;
        }
        _gr_tower_lazy_fresh(res, x->F, x->level, ctx);
        fmpz_mpoly_q_set(&res->elem.flat.data, &x->elem.flat.data, res->elem.flat.mctx);
        res->reduced_version = x->reduced_version;
    }
    return GR_SUCCESS;
}

int _gr_tower_lazy_zero(gr_tower_lazy_elem_t res, gr_ctx_t ctx) { _gr_tower_lazy_fresh(res, LAZY(ctx)->trivial, 0, ctx); return GR_SUCCESS; }
/* res made a rational element (of the trivial tower), its storage reused
   when it is one already */
static void
_gr_tower_lazy_fresh_rational(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F = LAZY(ctx)->trivial;
    if (res->F == F && res->repr == LAZY_REPR_RATIONAL && !res->shallow)
    {
        res->level = 0;
        res->structure_version = F->T->structure_version;
        res->reduced_version = 0;
    }
    else
        _gr_tower_lazy_fresh(res, F, 0, ctx);
}

/* x as a rational number (q shallow: valid while x is unchanged) */
static void
_gr_tower_lazy_rational_view(fmpq_t q, const gr_tower_lazy_elem_t x)
{
    *q = x->elem.q;
}

/* res (rational, see LAZY_IS_RATIONAL) = q, which is consumed */
static void
_gr_tower_lazy_rational_store(gr_tower_lazy_elem_t res, fmpq_t q)
{
    fmpq_swap(&res->elem.q, q);
    fmpq_clear(q);
}

/* a result which is rational (level 0) in another tower is moved to the
   trivial tower */
void
_gr_tower_lazy_normalize_rational(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    if (res->level == 0 && !res->shallow && res->repr == LAZY_REPR_FLAT &&
        fmpz_mpoly_q_is_fmpq(&res->elem.flat.data, res->elem.flat.mctx))
    {
        fmpq_t q;
        fmpq_init(q);
        (void) fmpz_mpoly_q_get_fmpq(q, &res->elem.flat.data, res->elem.flat.mctx);
        _gr_tower_lazy_fresh(res, L->trivial, 0, ctx);
        fmpq_swap(&res->elem.q, q);
        fmpq_clear(q);
    }
}

int _gr_tower_lazy_one(gr_tower_lazy_elem_t res, gr_ctx_t ctx) { _gr_tower_lazy_fresh_rational(res, ctx); fmpq_one(&res->elem.q); return GR_SUCCESS; }
int _gr_tower_lazy_set_si(gr_tower_lazy_elem_t res, slong c, gr_ctx_t ctx) { _gr_tower_lazy_fresh_rational(res, ctx); fmpq_set_si(&res->elem.q, c, 1); return GR_SUCCESS; }
int _gr_tower_lazy_set_fmpz(gr_tower_lazy_elem_t res, const fmpz_t c, gr_ctx_t ctx) { _gr_tower_lazy_fresh_rational(res, ctx); fmpq_set_fmpz(&res->elem.q, c); return GR_SUCCESS; }
int _gr_tower_lazy_set_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx) { _gr_tower_lazy_fresh_rational(res, ctx); fmpq_set(&res->elem.q, c); return GR_SUCCESS; }

static int
_gr_tower_lazy_set_ui(gr_tower_lazy_elem_t res, ulong c, gr_ctx_t ctx)
{
    fmpz_t t;
    fmpz_init_set_ui(t, c);
    _gr_tower_lazy_set_fmpz(res, t, ctx);
    fmpz_clear(t);
    return GR_SUCCESS;
}

/* Sets res to generator k of the tower of F. */
void
_gr_tower_lazy_set_gen(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong k, gr_ctx_t ctx)
{
    _gr_tower_lazy_fresh(res, F, k, ctx);
    fmpz_mpoly_gen(fmpz_mpoly_q_numref(&res->elem.flat.data), GR_TOWER_FLAT_VAR(F, k), res->elem.flat.mctx);
    fmpz_mpoly_one(fmpz_mpoly_q_denref(&res->elem.flat.data), res->elem.flat.mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(&res->elem.flat.data, F));
    res->reduced_version = F->ideal_version;
    _gr_tower_lazy_shrink(res);
}

/* Sets res to the generator with definition order d of F. */
void
_gr_tower_lazy_set_gen_d(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong d, gr_ctx_t ctx)
{
    _gr_tower_lazy_fresh(res, F, F->T->num_gens, ctx);
    fmpz_mpoly_q_gen(&res->elem.flat.data, GR_TOWER_FLAT_VAR_D(F, d), res->elem.flat.mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(&res->elem.flat.data, F));
    res->reduced_version = F->ideal_version;
    _gr_tower_lazy_shrink(res);
}

/* r mod 1 in [0, 1): exp(2 pi i r) depends only on that (keeps r canonical) */
void
_gr_tower_lazy_fmpq_frac_part(fmpq_t r)
{
    fmpz_fdiv_r(fmpq_numref(r), fmpq_numref(r), fmpq_denref(r));
}

/*
    The tower QQ(x) of the algebraic number x (hash-consed), and the gid
    of its generator. When a new tower is created, the generator is
    adjoined as a root of unity or as a root of a positive integer if
    kind is GR_TOWER_ROOT_OF_UNITY or GR_TOWER_ROOT (with the order n and
    the radicand p), or as tan(pi / n) if kind is GR_TOWER_TAN_PI, so that
    it carries that definition.
*/
int
_gr_tower_lazy_qqbar_tower(gr_tower_flat_struct ** F_out, slong * gid, const qqbar_t x, int kind, ulong n, const fmpz_t p, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_flat_struct * F;
    slong i;
    int status;

    for (i = 0; i < L->num_qqbars; i++)
    {
        if (qqbar_equal(&L->qqbars[i].x, x))
        {
            *F_out = L->qqbars[i].F;
            *gid = L->qqbars[i].gid;
            return GR_SUCCESS;
        }
    }

    F = _gr_tower_lazy_new_tower(ctx);
    if (kind == GR_TOWER_ROOT_OF_UNITY)
        status = gr_tower_adjoin_root_of_unity(F->T, n, NULL);
    else if (kind == GR_TOWER_ROOT)
        status = gr_tower_adjoin_root_fmpz(F->T, p, n, NULL);
    else
        status = gr_tower_adjoin_qqbar(F->T, x, NULL);
    if (status != GR_SUCCESS)
        return status;
    if (kind == GR_TOWER_TAN_PI)
    {
        /* (x = tan(pi / n)) */
        GR_TOWER_STEP(F->T, 0)->def_kind = GR_TOWER_TAN_PI;
        GR_TOWER_STEP(F->T, 0)->def_param = n;
    }
    _gr_tower_lazy_new_def(GR_TOWER_STEP(F->T, 0), F->T, ctx);

    if (L->num_qqbars == L->alloc_qqbars)
    {
        L->alloc_qqbars = FLINT_MAX(4, 2 * L->alloc_qqbars);
        L->qqbars = flint_realloc(L->qqbars, L->alloc_qqbars * sizeof(gr_tower_lazy_qqbar_entry_struct));
    }
    qqbar_init(&L->qqbars[L->num_qqbars].x);
    qqbar_set(&L->qqbars[L->num_qqbars].x, x);
    L->qqbars[L->num_qqbars].F = F;
    L->qqbars[L->num_qqbars].gid = GR_TOWER_STEP(F->T, 0)->gid;
    L->num_qqbars++;

    *F_out = F;
    *gid = GR_TOWER_STEP(F->T, 0)->gid;
    return GR_SUCCESS;
}

int
_gr_tower_lazy_set_qqbar(gr_tower_lazy_elem_t res, const qqbar_t x, gr_ctx_t ctx)
{
    gr_tower_flat_struct * F;
    slong gid;
    int status;

    if (qqbar_is_rational(x))
    {
        fmpq_t c;
        fmpq_init(c);
        qqbar_get_fmpq(c, x);
        status = _gr_tower_lazy_set_fmpq(res, c, ctx);
        fmpq_clear(c);
        return status;
    }

    {
        slong p;
        ulong q;
        if (qqbar_is_root_of_unity(&p, &q, x))
            return _gr_tower_lazy_root_of_unity(res, p, q, ctx);
    }

    status = _gr_tower_lazy_qqbar_tower(&F, &gid, x, GR_TOWER_ALGEBRAIC, 0, NULL, ctx);
    if (status == GR_SUCCESS)
        _gr_tower_lazy_set_gen_d(res, F, gr_tower_gid_order(F->T, gid), ctx);
    return status;
}

/* x of degree 2 as (-b +/- sqrt(b^2 - 4ac)) / (2a), in terms of the
   principal square root (the structured form, used for the roots of
   polynomials; conversions of qqbar elements keep the generic form,
   which merges more cheaply) */
static int
_gr_tower_lazy_set_qqbar_quadratic(gr_tower_lazy_elem_t res, const qqbar_t x, gr_ctx_t ctx)
{
    int status;

    const fmpz * c = QQBAR_COEFFS(x);
    fmpq_t d, u;
    qqbar_t y;
    gr_tower_lazy_elem_struct s;
    int sign;

    fmpq_init(d);
    fmpq_init(u);
    qqbar_init(y);
    _gr_tower_lazy_init(&s, ctx);

    fmpz_mul(fmpq_numref(d), c + 1, c + 1);
    fmpz_submul(fmpq_numref(d), c + 2, c + 0);
    fmpz_submul(fmpq_numref(d), c + 2, c + 0);
    fmpz_submul(fmpq_numref(d), c + 2, c + 0);
    fmpz_submul(fmpq_numref(d), c + 2, c + 0);

    /* the sign: x == (-b + sqrt(d)) / (2a)? Numerically: 2 a x + b is
       sqrt(d) or -sqrt(d), which differ (d != 0), so the enclosures
       separate at some precision. */
    {
        acb_t w, sq;
        slong prec;

        acb_init(w);
        acb_init(sq);
        sign = 0;
        for (prec = 64; prec <= 65536 && sign == 0; prec *= 2)
        {
            qqbar_get_acb(w, x, prec);
            acb_mul_fmpz(w, w, c + 2, prec);
            acb_mul_2exp_si(w, w, 1);
            acb_add_fmpz(w, w, c + 1, prec);
            acb_set_fmpq(sq, d, prec);
            acb_sqrt(sq, sq, prec);
            if (!acb_overlaps(w, sq))
                sign = -1;
            else
            {
                acb_neg(sq, sq);
                if (!acb_overlaps(w, sq))
                    sign = 1;
            }
        }
        acb_clear(w);
        acb_clear(sq);

        if (sign == 0)
        {
            /* (not reached in practice) exact comparison */
            qqbar_set_fmpq(y, d);
            qqbar_sqrt(y, y);
            qqbar_sub_fmpz(y, y, c + 1);
            qqbar_div_fmpz(y, y, c + 2);
            qqbar_div_ui(y, y, 2);
            sign = qqbar_equal(x, y) ? 1 : -1;
        }
    }

    status = _gr_tower_lazy_root_fmpq(&s, d, 2, ctx);
    if (status == GR_SUCCESS)
    {
        fmpz_set_si(fmpq_numref(u), sign);
        fmpz_mul_ui(fmpq_denref(u), c + 2, 2);
        fmpq_canonicalise(u);
        status |= gr_mul_fmpq(res, &s, u, ctx);
        fmpz_neg(fmpq_numref(u), c + 1);
        fmpz_mul_ui(fmpq_denref(u), c + 2, 2);
        fmpq_canonicalise(u);
        status |= gr_add_fmpq(res, res, u, ctx);
    }

    _gr_tower_lazy_clear(&s, ctx);
    qqbar_clear(y);
    fmpq_clear(u);
    fmpq_clear(d);
    return status;
}

/*
    x as a structured element where possible: a quadratic irrationality
    through a square root, a root of a binomial a x^n + b through the
    principal root of -b/a times a root of unity; otherwise the generic
    algebraic generator (used for the roots of polynomials and for the
    real forms of algebraic elements, not for conversions in general;
    see _gr_tower_lazy_set_qqbar_quadratic).
*/
int
_gr_tower_lazy_set_qqbar_structured(gr_tower_lazy_elem_t res, const qqbar_t x, gr_ctx_t ctx)
{
    slong d = qqbar_degree(x), i;
    const fmpz * c = QQBAR_COEFFS(x);
    int binomial = (d >= 2);

    if (d == 2)
        return _gr_tower_lazy_set_qqbar_quadratic(res, x, ctx);

    for (i = 1; i < d && binomial; i++)
        binomial = fmpz_is_zero(c + i);

    if (binomial)
    {
        fmpq_t q;
        qqbar_t r, u;
        slong p;
        ulong n;
        int status = GR_UNABLE;

        fmpq_init(q);
        qqbar_init(r);
        qqbar_init(u);

        /* r = principal d-th root of -c_0 / c_d, u = x / r */
        fmpz_neg(fmpq_numref(q), c + 0);
        fmpz_set(fmpq_denref(q), c + d);
        fmpq_canonicalise(q);
        qqbar_set_fmpq(r, q);
        qqbar_root_ui(r, r, d);
        qqbar_div(u, x, r);

        if (qqbar_is_root_of_unity(&p, &n, u))
        {
            gr_tower_lazy_elem_struct z;
            _gr_tower_lazy_init(&z, ctx);
            status = _gr_tower_lazy_root_fmpq(res, q, d, ctx);
            if (status == GR_SUCCESS && !(p == 0 || n == 1))
            {
                status = _gr_tower_lazy_root_of_unity(&z, p, n, ctx);
                if (status == GR_SUCCESS)
                    status = _gr_tower_lazy_mul(res, res, &z, ctx);
            }
            _gr_tower_lazy_clear(&z, ctx);
        }

        qqbar_clear(r);
        qqbar_clear(u);
        fmpq_clear(q);

        if (status == GR_SUCCESS)
            return status;
    }

    return _gr_tower_lazy_set_qqbar(res, x, ctx);
}

/* res = the e-th power of the generator with the given gid of F. */
static void
_gr_tower_lazy_set_gen_pow(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong gid, ulong e, gr_ctx_t ctx)
{
    _gr_tower_lazy_fresh(res, F, F->T->num_gens, ctx);
    fmpz_mpoly_q_gen(&res->elem.flat.data, GR_TOWER_FLAT_VAR_D(F, gr_tower_gid_order(F->T, gid)), res->elem.flat.mctx);
    if (e != 1)
        fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(&res->elem.flat.data), fmpz_mpoly_q_numref(&res->elem.flat.data), e, res->elem.flat.mctx);
    GR_MUST_SUCCEED(gr_tower_flat_reduce(&res->elem.flat.data, F));
    res->reduced_version = F->ideal_version;
    _gr_tower_lazy_shrink(res);
}

/*
    The canonical root of unity of order l^e (p = NULL) or the principal
    l^e-th root of the positive integer p, where l is prime: the
    generator of the largest order l^f (f >= e) of its kind, raised to
    the power l^(f-e). Each kind has an entry in the registry of the
    context, whose tower may have grown (and whose canonical generator
    may have been replaced by one of a larger order when towers were
    merged, see _structured_image in map.c); a root of a larger order
    than any present replaces the entry by a new tower.
*/
int
_gr_tower_lazy_prime_power_root(gr_tower_lazy_elem_t res, const fmpz_t p, ulong l, ulong e, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    gr_tower_lazy_root_entry_struct * ent = NULL;
    slong i;
    int status;

    for (i = 0; i < L->num_roots; i++)
    {
        if (L->roots[i].l == l && ((p == NULL) ? fmpz_is_zero(&L->roots[i].p) : fmpz_equal(&L->roots[i].p, p)))
        {
            ent = L->roots + i;
            break;
        }
    }

    if (ent != NULL)
    {
        /* the generator of the largest order l^f in the entry's tower */
        gr_tower_struct * T = ent->F->T;
        fmpz_t pu;
        ulong m, best_f = 0;
        slong d, best_gid = -1;
        int kind = (p == NULL) ? 1 : 2;

        fmpz_init(pu);
        for (d = 0; d < T->num_gens; d++)
        {
            if (_gr_tower_gen_const_root(GR_TOWER_GEN(T, d), pu, &m) == kind && (kind == 1 || fmpz_equal(p, pu)))
            {
                ulong f = 0;
                while (m % l == 0)
                {
                    m /= l;
                    f++;
                }
                if (m == 1 && f > best_f)
                {
                    best_f = f;
                    best_gid = GR_TOWER_GEN(T, d)->gid;
                }
            }
        }
        fmpz_clear(pu);

        if (best_gid >= 0 && best_f >= e)
        {
            ent->gid = best_gid;
            _gr_tower_lazy_set_gen_pow(res, ent->F, best_gid, n_pow(l, best_f - e), ctx);
            return GR_SUCCESS;
        }
    }

    /* a new canonical generator of order l^e */
    {
        gr_tower_flat_struct * F;
        qqbar_t x;
        slong gid;
        ulong q = n_pow(l, e);

        qqbar_init(x);
        if (p == NULL)
            qqbar_root_of_unity(x, 1, q);
        else
        {
            qqbar_set_fmpz(x, p);
            qqbar_root_ui(x, x, q);
        }

        status = _gr_tower_lazy_qqbar_tower(&F, &gid, x, (p == NULL) ? GR_TOWER_ROOT_OF_UNITY : GR_TOWER_ROOT, q, p, ctx);
        qqbar_clear(x);
        if (status != GR_SUCCESS)
            return status;

        if (ent == NULL)
        {
            if (L->num_roots == L->alloc_roots)
            {
                L->alloc_roots = FLINT_MAX(4, 2 * L->alloc_roots);
                L->roots = flint_realloc(L->roots, L->alloc_roots * sizeof(gr_tower_lazy_root_entry_struct));
            }
            ent = L->roots + L->num_roots;
            fmpz_init(&ent->p);
            L->num_roots++;
        }

        if (p != NULL)
            fmpz_set(&ent->p, p);
        else
            fmpz_zero(&ent->p);
        ent->l = l;
        ent->F = F;
        ent->gid = gid;
        F->gc |= GR_TOWER_GC_PINNED;

        _gr_tower_lazy_set_gen_pow(res, F, gid, 1, ctx);
        return GR_SUCCESS;
    }
}

/*
    Sets res to the root of unity exp(2 pi i / n) (p = NULL) or to the
    principal n-th root of the positive integer p, as a product of
    powers of the canonical roots of prime power orders q_i = l_i^{e_i}
    (n = prod q_i): with c_i = (n/q_i)^{-1} mod q_i, sum c_i / q_i =
    1/n + k for an integer k >= 0, so that root_n = p^{-k} prod
    root_{q_i}^{c_i}. The roots of prime power orders of different
    primes generate linearly disjoint fields, so this keeps the
    representation canonical (irreducible moduli) and sparse.
*/
int
_gr_tower_lazy_const_root(gr_tower_lazy_elem_t res, const fmpz_t p, ulong n, gr_ctx_t ctx)
{
    n_factor_t fac;
    gr_tower_lazy_elem_struct t;
    slong i;
    ulong ksum = 0;
    int status = GR_SUCCESS;

    if (n == 1)
    {
        if (p == NULL)
            return _gr_tower_lazy_one(res, ctx);
        return _gr_tower_lazy_set_fmpz(res, p, ctx);
    }

    n_factor_init(&fac);
    n_factor(&fac, n, 1);

    status = _gr_tower_lazy_one(res, ctx);
    _gr_tower_lazy_init(&t, ctx);

    for (i = 0; i < fac.num && status == GR_SUCCESS; i++)
    {
        ulong l = fac.p[i], e = fac.exp[i], q = n_pow(l, e), c;

        status = _gr_tower_lazy_prime_power_root(&t, p, l, e, ctx);
        if (status != GR_SUCCESS)
            break;

        c = (fac.num == 1) ? 1 : n_invmod((n / q) % q, q);
        ksum += c * (n / q);
        if (c != 1)
            status = gr_pow_ui(&t, &t, c, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_mul(res, res, &t, ctx);
    }

    if (status == GR_SUCCESS && p != NULL && fac.num > 1)
    {
        ulong k = (ksum - 1) / n;
        if (k > 0)
        {
            fmpz_t pk;
            fmpz_init(pk);
            fmpz_pow_ui(pk, p, k);
            status = gr_div_fmpz(res, res, pk, ctx);
            fmpz_clear(pk);
        }
    }

    _gr_tower_lazy_clear(&t, ctx);
    return status;
}

/*
    GR_TOWER_GENS_COMPOSITE_ROOTS: zeta_q^p (q >= 3, 0 < p < q) as a power
    of the root of unity of the largest order divisible by q in the tower
    of the canonical root of unity (whose generator may have been
    replaced by one of a larger order when towers were merged); otherwise
    a new canonical root of unity of the least common multiple of q and
    the orders present (towers holding the old one identify it as a power
    of the new one when merged, see _structured_image in map.c).
*/
static int
_gr_tower_lazy_composite_root_of_unity(gr_tower_lazy_elem_t res, ulong p, ulong q, gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    ulong M = q;
    slong d;

    if (L->cyclo_F != NULL)
    {
        gr_tower_struct * T = L->cyclo_F->T;
        slong best = -1;
        ulong bestN = 0;

        for (d = 0; d < T->num_gens; d++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
            if (g->kind == GR_TOWER_ALGEBRAIC && g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param > 0)
            {
                ulong N = g->def_param;
                /* (capped: 0 when beyond the order limit) */
                if (M != 0)
                {
                    ulong a = M / n_gcd(M, N);
                    M = (a > (ulong) L->options[GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT] / N) ? 0 : a * N;
                }
                if (N % q == 0 && N > bestN)
                {
                    bestN = N;
                    best = g->gid;
                }
            }
        }

        if (best >= 0)
        {
            _gr_tower_lazy_set_gen_pow(res, L->cyclo_F, best, p * (bestN / q), ctx);
            return GR_SUCCESS;
        }
    }

    /* (beyond the limits: the prime power decomposition) */
    if (M == 0 || M > (ulong) L->options[GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT] || n_euler_phi(M) > (ulong) L->options[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT])
        return GR_UNABLE;

    {
        gr_tower_flat_struct * F;
        qqbar_t x;
        slong gid;
        int status;

        qqbar_init(x);
        qqbar_root_of_unity(x, 1, M);
        status = _gr_tower_lazy_qqbar_tower(&F, &gid, x, GR_TOWER_ROOT_OF_UNITY, M, NULL, ctx);
        qqbar_clear(x);
        if (status != GR_SUCCESS)
            return status;
        F->gc |= GR_TOWER_GC_PINNED;
        L->cyclo_F = F;
        _gr_tower_lazy_set_gen_pow(res, F, gid, p * (M / q), ctx);
        return GR_SUCCESS;
    }
}

/* res = exp(2 pi i p / q) */
int
_gr_tower_lazy_root_of_unity(gr_tower_lazy_elem_t res, slong p, ulong q, gr_ctx_t ctx)
{
    ulong g, pp;
    int status;

    if (q == 0)
        return GR_DOMAIN;

    /* reduce p/q modulo 1 */
    if (p < 0)
        pp = q - ((ulong) (-p) % q);
    else
        pp = (ulong) p % q;
    if (pp == q)
        pp = 0;
    g = n_gcd(pp, q);
    if (g == 0)
        return _gr_tower_lazy_one(res, ctx);
    pp /= g;
    q /= g;

    if (q == 1)
        return _gr_tower_lazy_one(res, ctx);
    if (q == 2)
        return _gr_tower_lazy_set_si(res, -1, ctx);

    if (LAZY(ctx)->options[GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT] > 0)
    {
        status = _gr_tower_lazy_composite_root_of_unity(res, pp, q, ctx);
        if (status != GR_UNABLE)
            return status;
    }

    /* zeta_q^pp = prod zeta_{q_i}^{pp c_i} over the prime power factors
       q_i of q, with c_i = (q/q_i)^{-1} mod q_i; the powers are taken
       of the (univariate) generators before multiplying */
    {
        n_factor_t fac;
        gr_tower_lazy_elem_struct t;
        slong i;

        n_factor_init(&fac);
        n_factor(&fac, q, 1);

        status = _gr_tower_lazy_one(res, ctx);
        _gr_tower_lazy_init(&t, ctx);

        for (i = 0; i < fac.num && status == GR_SUCCESS; i++)
        {
            ulong l = fac.p[i], e = fac.exp[i], qi = n_pow(l, e), c;

            status = _gr_tower_lazy_prime_power_root(&t, NULL, l, e, ctx);
            if (status != GR_SUCCESS)
                break;

            c = (fac.num == 1) ? 1 : n_invmod((q / qi) % qi, qi);
            c = n_mulmod2_preinv(c, pp % qi, qi, n_preinvert_limb(qi));
            if (c == 0)
                continue;
            if (c != 1)
            {
                fmpz_mpoly_pow_ui(fmpz_mpoly_q_numref(&t.elem.flat.data), fmpz_mpoly_q_numref(&t.elem.flat.data), c, t.elem.flat.mctx);
                status = gr_tower_flat_reduce(&t.elem.flat.data, t.F);
                _gr_tower_lazy_shrink(&t);
            }
            if (status == GR_SUCCESS)
                status = _gr_tower_lazy_mul(res, res, &t, ctx);
        }

        _gr_tower_lazy_clear(&t, ctx);
    }
    return status;
}

static int
_gr_tower_lazy_set_other(gr_tower_lazy_elem_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
        return _gr_tower_lazy_set(res, x, ctx);

    if (x_ctx->which_ring == GR_CTX_REAL_ALGEBRAIC_QQBAR || x_ctx->which_ring == GR_CTX_COMPLEX_ALGEBRAIC_QQBAR)
        return _gr_tower_lazy_set_qqbar(res, x, ctx);

    if (x_ctx->which_ring == GR_CTX_GR_TOWER_LAZY && LAZY(x_ctx) == LAZY(ctx))
        return _gr_tower_lazy_set(res, x, ctx);   /* (a view of the same state) */

    if (x_ctx->which_ring == GR_CTX_GR_TOWER_LAZY)
        return _gr_tower_lazy_transfer(res, x, x_ctx, ctx);

    if (x_ctx->which_ring == GR_CTX_FEXPR)
    {
        /* symbolic expressions are evaluated with the operations of the
           field */
        fexpr_vec_t inputs;
        gr_vec_t outputs;
        int status;
        fexpr_vec_init(inputs, 0);
        gr_vec_init(outputs, 0, ctx);
        status = gr_generic_set_fexpr(res, inputs, outputs, x, ctx);
        fexpr_vec_clear(inputs);
        gr_vec_clear(outputs, ctx);
        return status;
    }

    {
        fmpq_t c;
        int status;
        fmpq_init(c);
        status = gr_get_fmpq(c, x, x_ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_set_fmpq(res, c, ctx);
        fmpq_clear(c);
        return status;
    }
}

/* Element of the nested field F_k of x's tower, k = algebraic level of x. */
int
_gr_tower_lazy_get_nested(gr_ptr res, gr_tower_lazy_elem_t x)
{
    x = _gr_tower_lazy_flat_view(x);
    return gr_tower_flat_get_nested_at(res, &x->elem.flat.data, _gr_tower_lazy_alg_level(x), x->F);
}

static int _gr_tower_lazy_get_qqbar_inner(qqbar_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* (reduces x in place: a shallow copy is replaced by a private copy) */
int
_gr_tower_lazy_get_qqbar(qqbar_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;
    if (x->shallow)
    {
        gr_tower_lazy_elem_struct t;
        _gr_tower_lazy_init(&t, ctx);
        status = _gr_tower_lazy_set(&t, x, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_get_qqbar_inner(res, &t, ctx);
        _gr_tower_lazy_clear(&t, ctx);
        return status;
    }
    return _gr_tower_lazy_get_qqbar_inner(res, x, ctx);
}

static int
_gr_tower_lazy_get_qqbar_inner(qqbar_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    gr_tower_struct * T = TOWER(x);
    gr_ctx_struct * Fk;
    gr_ptr t;
    slong k;
    int status;

    x = _gr_tower_lazy_flat_view(x);

    if (T->num_trans > 0)
    {
        /* the algebraic part of the tower, as a tower over QQ: possible
           when neither x nor the moduli involve the transcendental
           generators */
        gr_tower_flat_struct * F = x->F;
        gr_tower_t U;
        slong * c;
        int * used;
        slong i;
        fmpz_mpoly_q_t xu;
        int ok = 1;

        gr_tower_flat_ensure(F);
        GR_MUST_SUCCEED(gr_tower_flat_reduce(&x->elem.flat.data, F));

        used = flint_malloc(sizeof(int) * F->cap);
        fmpz_mpoly_q_used_vars(used, &x->elem.flat.data, F->mctx);
        for (i = 1; i <= T->num_trans; i++)
            if (used[GR_TOWER_FLAT_TVAR(F, i)])
                ok = 0;
        for (k = 1; k <= T->length && ok; k++)
        {
            fmpz_mpoly_used_vars(used, (fmpz_mpoly_struct *) _gr_tower_flat_ideal_elem(F, k), F->mctx);
            for (i = 1; i <= T->num_trans; i++)
                if (used[GR_TOWER_FLAT_TVAR(F, i)])
                    ok = 0;
        }
        flint_free(used);

        if (!ok)
            return GR_UNABLE;

        /* variable c[v] of U for variable v of F */
        gr_tower_init(U, T->consts);
        U->keep_retired = 0;   /* (internal: no nested elements across rebuilds) */
        U->options = T->options;
        gr_tower_flat_ensure(&U->flat);
        {
            slong cap = FLINT_MAX(T->length, 1);
            c = flint_malloc(sizeof(slong) * F->cap);
            for (i = 0; i < F->cap; i++)
                c[i] = 0;
            /* U will have T->length generators: definition order k - 1 for step k */
            gr_tower_flat_clear(&U->flat);
            gr_tower_flat_init(&U->flat, U, cap);
            for (k = 1; k <= T->length; k++)
                c[GR_TOWER_FLAT_VAR(F, k)] = cap - 1 - (k - 1);
        }

        status = GR_SUCCESS;
        for (k = 1; k <= T->length && status == GR_SUCCESS; k++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_STEP(T, k - 1);
            gr_ctx_struct * top = gr_tower_field(U);
            fmpz_mpoly_t mu;
            gr_poly_t m;
            slong v = c[GR_TOWER_FLAT_VAR(F, k)];
            fmpz_mpoly_univar_t uni;
            fmpz_mpoly_t coeff;

            fmpz_mpoly_init(mu, U->flat.mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(mu, _gr_tower_flat_ideal_elem(F, k), c, F->mctx, U->flat.mctx);

            /* the modulus, coefficient by coefficient in a_k */
            gr_poly_init(m, top);
            fmpz_mpoly_univar_init(uni, U->flat.mctx);
            fmpz_mpoly_init(coeff, U->flat.mctx);
            fmpz_mpoly_to_univar(uni, mu, v, U->flat.mctx);
            for (i = 0; i < fmpz_mpoly_univar_length(uni, U->flat.mctx) && status == GR_SUCCESS; i++)
            {
                slong e = fmpz_mpoly_univar_get_term_exp_si(uni, i, U->flat.mctx);
                gr_ptr ce;
                GR_TMP_INIT(ce, top);
                fmpz_mpoly_univar_get_term_coeff(coeff, uni, i, U->flat.mctx);
                status |= gr_tower_flat_poly_get_nested_at(ce, coeff, U->length, &U->flat);
                status |= gr_poly_set_coeff_scalar(m, e, ce, top);
                GR_TMP_CLEAR(ce, top);
            }
            status |= gr_poly_make_monic(m, m, top);
            if (status == GR_SUCCESS)
                status = gr_tower_adjoin_algebraic(U, m, &g->enclosure, g->status, g->name);

            gr_poly_clear(m, top);
            fmpz_mpoly_univar_clear(uni, U->flat.mctx);
            fmpz_mpoly_clear(coeff, U->flat.mctx);
            fmpz_mpoly_clear(mu, U->flat.mctx);
        }

        if (status == GR_SUCCESS)
        {
            gr_ctx_struct * top = gr_tower_field(U);
            gr_ptr u;
            GR_TMP_INIT(u, top);
            fmpz_mpoly_q_init(xu, U->flat.mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_numref(xu), fmpz_mpoly_q_numref(&x->elem.flat.data), c, F->mctx, U->flat.mctx);
            fmpz_mpoly_compose_fmpz_mpoly_gen(fmpz_mpoly_q_denref(xu), fmpz_mpoly_q_denref(&x->elem.flat.data), c, F->mctx, U->flat.mctx);
            status = gr_tower_flat_get_nested_at(u, xu, U->length, &U->flat);
            if (status == GR_SUCCESS)
                status = gr_tower_get_qqbar(res, u, U);
            fmpz_mpoly_q_clear(xu, U->flat.mctx);
            GR_TMP_CLEAR(u, top);
        }

        flint_free(c);
        gr_tower_clear(U);
        return status;
    }

    k = _gr_tower_lazy_alg_level(x);
    Fk = gr_tower_field_at(T, k);
    GR_TMP_INIT(t, Fk);
    status = _gr_tower_lazy_get_nested(t, x);
    if (status == GR_SUCCESS)
    {
        if (k == 0)
        {
            fmpq_t c;
            fmpq_init(c);
            status = gr_get_fmpq(c, t, T->base);
            if (status == GR_SUCCESS)
                qqbar_set_fmpq(res, c);
            fmpq_clear(c);
        }
        else
        {
            gr_ptr u;
            GR_TMP_INIT(u, gr_tower_field(T));
            status = gr_tower_promote(u, t, k, T->length, T);
            if (status == GR_SUCCESS)
                status = gr_tower_get_qqbar(res, u, T);
            GR_TMP_CLEAR(u, gr_tower_field(T));
        }
    }
    GR_TMP_CLEAR(t, Fk);
    return status;
}

static int
_gr_tower_lazy_randtest(gr_tower_lazy_elem_t res, flint_rand_t state, gr_ctx_t ctx)
{
    if (n_randint(state, 3) == 0)
    {
        fmpq_t c;
        int status;
        fmpq_init(c);
        fmpq_randtest(c, state, 6);
        status = _gr_tower_lazy_set_fmpq(res, c, ctx);
        fmpq_clear(c);
        return status;
    }
    else
    {
        qqbar_t x;
        int status;
        qqbar_init(x);
        qqbar_randtest(x, state, 1 + n_randint(state, 3), 6);
        status = _gr_tower_lazy_set_qqbar(res, x, ctx);
        qqbar_clear(x);
        return status;
    }
}

/* the definition orders of the marked generators, sorted by creation
   (definition id): a definition depends only on earlier ones */
slong *
_gr_tower_lazy_marked_by_creation(slong * count, const int * mark, gr_tower_t T)
{
    slong d, n = 0, i, j;
    slong * order = flint_malloc(sizeof(slong) * FLINT_MAX(T->num_gens, 1));
    for (d = 0; d < T->num_gens; d++)
        if (mark[d])
            order[n++] = d;
    /* insertion sort by def_id (unassigned ids last, in tower order) */
    for (i = 1; i < n; i++)
    {
        slong x = order[i];
        ulong key = GR_TOWER_GEN(T, x)->def_id;
        for (j = i - 1; j >= 0; j--)
        {
            ulong kj = GR_TOWER_GEN(T, order[j])->def_id;
            if (kj != 0 && (key == 0 || kj <= key))
                break;
            order[j + 1] = order[j];
        }
        order[j + 1] = x;
    }
    *count = n;
    return order;
}

/* pi and i print as themselves */
int
_gr_tower_lazy_gen_is_constant_symbol(const gr_tower_gen_struct * g)
{
    return g->def_kind == GR_TOWER_PI || (g->def_kind == GR_TOWER_ROOT_OF_UNITY && g->def_param == 4);
}

/* marks (in mark[], by definition order) the generators of T which a
   flat element in the context mctx involves */
void
_gr_tower_lazy_mark_used_gens(int * mark, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_struct * mctx, gr_tower_t T)
{
    slong cap, v;
    slong * var_gid;
    int * used;

    if (!_gr_tower_flat_find_layout(&cap, &var_gid, mctx, &T->flat))
        return;

    used = flint_malloc(sizeof(int) * cap);
    fmpz_mpoly_q_used_vars(used, x, mctx);
    for (v = 0; v < cap; v++)
    {
        if (used[v] && var_gid[v] >= 0)
        {
            slong d = gr_tower_gid_order(T, var_gid[v]);
            if (d >= 0)
                mark[d] = 1;
        }
    }
    flint_free(used);
}

/* the generators involved in x, with the generators their definitions
   involve (recursively), by definition order */

void
_gr_tower_lazy_involved_gens(int * mark, const gr_tower_lazy_elem_t x, gr_tower_t T)
{
    slong d;

    for (d = 0; d < T->num_gens; d++)
        mark[d] = 0;
    _gr_tower_lazy_mark_used_gens(mark, &x->elem.flat.data, x->elem.flat.mctx, T);
    _gr_tower_involved_gens_closure(mark, T);
}

/* adds to the marked generators those their definitions involve */
void
_gr_tower_involved_gens_closure(int * mark, gr_tower_t T)
{
    slong d;
    int progress = 1;

    while (progress)
    {
        progress = 0;
        for (d = 0; d < T->num_gens; d++)
        {
            const gr_tower_gen_struct * g = GR_TOWER_GEN(T, d);
            if (mark[d] == 1 && g->arg.mctx != NULL)
            {
                slong e;
                int * mark2 = flint_calloc(T->num_gens, sizeof(int));
                _gr_tower_lazy_mark_used_gens(mark2, &g->arg.data, g->arg.mctx, T);
                for (e = 0; e < T->num_gens; e++)
                    if (mark2[e] && !mark[e])
                        mark[e] = 1, progress = 1;
                flint_free(mark2);
            }
            if (mark[d] == 1 && g->kind == GR_TOWER_ALGEBRAIC)
            {
                /* a root of its minimal polynomial: the generators in
                   the coefficients (also for a generator with a
                   definition: the modulus may involve others, after a
                   refinement or a relation, sqrt(c) = an expression) */
                slong e, k = g->index;
                const gr_poly_struct * m = gr_tower_step_minpoly(T, k);
                int * mark2 = flint_calloc(T->num_gens, sizeof(int));
                fmpz_mpoly_q_t c;

                fmpz_mpoly_q_init(c, T->flat.mctx);
                for (e = 0; e < m->length; e++)
                {
                    if (gr_tower_flat_set_nested_at(c, gr_poly_coeff_srcptr(m, e, gr_tower_field_at(T, k - 1)), k - 1, &T->flat) == GR_SUCCESS)
                        _gr_tower_lazy_mark_used_gens(mark2, c, T->flat.mctx, T);
                    else
                    {
                        /* (not expected for a polynomial conversion;
                           conservatively, every generator the
                           coefficient may involve) */
                        slong e2;
                        for (e2 = 0; e2 < d; e2++)
                            mark2[e2] = 1;
                    }
                }
                fmpz_mpoly_q_clear(c, T->flat.mctx);

                for (e = 0; e < T->num_gens; e++)
                    if (mark2[e] && !mark[e])
                        mark[e] = 1, progress = 1;
                flint_free(mark2);
            }
            if (mark[d] == 1)
                mark[d] = 2;   /* processed */
        }
    }
}

/* -------------------------------------------------------------------- */
/* locking layer: every operation on the context (which may refine or   */
/* merge the towers of the registry, and update elements) holds the     */
/* recursive mutex of the context                                       */
/* -------------------------------------------------------------------- */

#define LOCK_UNARY(pub, name) int pub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(name(res, x, ctx)) }
#define LOCK_PRED(name) static truth_t name##_locked(const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED_T(truth_t, name(x, ctx)) }

int gr_tower_lazy_poly_mullow(gr_tower_lazy_elem_struct * res, const gr_tower_lazy_elem_struct * p1, slong len1, const gr_tower_lazy_elem_struct * p2, slong len2, slong n, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_poly_mullow(res, p1, len1, p2, len2, n, ctx)) }
int gr_tower_lazy_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_mat_mul(C, A, B, ctx)) }
int gr_tower_lazy_mat_det(gr_tower_lazy_elem_t res, const gr_mat_t A, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_mat_det(res, A, ctx)) }
int gr_tower_lazy_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_mat_nonsingular_solve(X, A, B, ctx)) }
static void _gr_tower_lazy_clear_locked(gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED_V(_gr_tower_lazy_clear(x, ctx)) }
int gr_tower_lazy_randtest(gr_tower_lazy_elem_t res, flint_rand_t state, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_randtest(res, state, ctx)) }
int gr_tower_lazy_write(gr_stream_t out, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_write(out, x, ctx)) }
/* -------------------------------------------------------------------- */
/* the real and algebraic subfields: operations whose result would     */
/* leave the subfield fail with GR_DOMAIN (all elements of a real       */
/* context are real, so the checks are sign tests; the special values   */
/* of the elementary functions at which they are algebraic are the only */
/* ones allowed in an algebraic context, by Lindemann-Weierstrass and   */
/* Gelfond-Schneider)                                                   */
/* -------------------------------------------------------------------- */

#define IS_REAL_CTX(ctx) ((VIEW(ctx)->field_flags & GR_TOWER_LAZY_REAL) && _gr_tower_lazy_outermost(ctx))
#define IS_ALG_CTX(ctx) ((VIEW(ctx)->field_flags & GR_TOWER_LAZY_ALGEBRAIC) && _gr_tower_lazy_outermost(ctx))

/* the sign of x in a real context: GR_SUCCESS with the sign, else the
   status (GR_UNABLE) of the sign test */
static int
_real_sign(int * sgn, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_real_sign_locked(sgn, x, ctx);
}

/* status of a result which must be real (set_other, roots) */
static int
_check_real(const gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    truth_t t = _gr_tower_lazy_is_real_exact(res, ctx);
    return (t == T_TRUE) ? GR_SUCCESS : ((t == T_FALSE) ? GR_DOMAIN : GR_UNABLE);
}

/* whether x depends only on algebraic generators (its definitions
   included) */
truth_t
_gr_tower_lazy_is_algebraic_repr(gr_tower_lazy_elem_t x)
{
    gr_tower_struct * T;
    int * mark;
    slong d;
    truth_t res = T_TRUE;

    x = _gr_tower_lazy_flat_view(x);
    T = TOWER(x);
    mark = flint_calloc(FLINT_MAX(T->num_gens, 1), sizeof(int));
    _gr_tower_lazy_involved_gens(mark, x, T);
    for (d = 0; d < T->num_gens; d++)
        if (mark[d] && GR_TOWER_GEN(T, d)->kind != GR_TOWER_ALGEBRAIC)
            res = T_UNKNOWN;
    flint_free(mark);
    return res;
}

truth_t
_gr_tower_lazy_is_algebraic_repr_locked(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    truth_t res;
    gr_tower_lazy_elem_struct t;
    _gr_tower_lazy_lock(ctx);
    _gr_tower_lazy_init(&t, ctx);
    if (_gr_tower_lazy_set(&t, x, ctx) == GR_SUCCESS)
        res = _gr_tower_lazy_is_algebraic_repr(&t);
    else
        res = T_UNKNOWN;
    _gr_tower_lazy_clear(&t, ctx);
    _gr_tower_lazy_unlock(ctx);
    return res;
}

/* status of a value entering the field of a view from outside the
   operations (parsing, conversions): GR_DOMAIN if it is not in the
   field, GR_UNABLE if that cannot be decided */
int
_gr_tower_lazy_check_member(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;

    if (VIEW(ctx)->field_flags & GR_TOWER_LAZY_REAL)
    {
        status = _check_real(res, ctx);
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_realify(res, ctx);
    }
    if (status == GR_SUCCESS && (VIEW(ctx)->field_flags & GR_TOWER_LAZY_ALGEBRAIC) &&
            _gr_tower_lazy_is_algebraic_repr(res) != T_TRUE)
        status = GR_UNABLE;
    return status;
}

int
gr_tower_lazy_pi(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = IS_ALG_CTX(ctx) ? GR_DOMAIN : _gr_tower_lazy_pi(res, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_i(gr_tower_lazy_elem_t res, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = IS_REAL_CTX(ctx) ? GR_DOMAIN : _gr_tower_lazy_i(res, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* a principal n-th root (n >= 2) of a real x is real iff x >= 0 */
int
gr_tower_lazy_root_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong n, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    if (IS_REAL_CTX(ctx) && n >= 2)
    {
        int sgn;
        status = _real_sign(&sgn, x, ctx);
        if (status == GR_SUCCESS && sgn < 0)
            status = GR_DOMAIN;
    }
    else
        status = GR_SUCCESS;
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_root_ui(res, x, n, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_sqrt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    if (IS_REAL_CTX(ctx))
    {
        int sgn;
        status = _real_sign(&sgn, x, ctx);
        if (status == GR_SUCCESS && sgn < 0)
            status = GR_DOMAIN;
    }
    else
        status = GR_SUCCESS;
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_sqrt(res, x, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_pow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (IS_REAL_CTX(ctx) || IS_ALG_CTX(ctx))
    {
        fmpq_t q;
        int st;
        fmpq_init(q);
        st = _gr_tower_lazy_get_fmpq(q, y, ctx);
        if (st == GR_SUCCESS)
        {
            /* a rational power: real for x >= 0 or an integer exponent */
            if (IS_REAL_CTX(ctx) && !fmpz_is_one(fmpq_denref(q)))
            {
                int sgn;
                status = _real_sign(&sgn, x, ctx);
                if (status == GR_SUCCESS && sgn < 0)
                    status = GR_DOMAIN;
            }
        }
        else if (IS_ALG_CTX(ctx))
        {
            /* x^y with y irrational: algebraic only for x = 0 or x = 1
               (when y is not known to be irrational, the answer is
               GR_UNABLE rather than GR_DOMAIN) */
            if (_gr_tower_lazy_is_one(x, ctx) == T_TRUE)
            {
                status = _gr_tower_lazy_one(res, ctx);
                fmpq_clear(q);
                _gr_tower_lazy_unlock(ctx);
                return status;
            }
            if (_gr_tower_lazy_is_zero(x, ctx) == T_TRUE)
                status = GR_SUCCESS;
            else
                status = (st == GR_DOMAIN) ? GR_DOMAIN : GR_UNABLE;
        }
        else
        {
            /* x^y = exp(y log x): real for x > 0 */
            int sgn;
            status = _real_sign(&sgn, x, ctx);
            if (status == GR_SUCCESS && sgn < 0)
                status = GR_DOMAIN;
        }
        fmpq_clear(q);
    }
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_pow(res, x, y, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_exp(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    if (IS_ALG_CTX(ctx))
    {
        truth_t z = _gr_tower_lazy_is_zero(x, ctx);
        status = (z == T_TRUE) ? _gr_tower_lazy_one(res, ctx) : ((z == T_FALSE) ? GR_DOMAIN : GR_UNABLE);
    }
    else
        status = _gr_tower_lazy_exp(res, x, ctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    _gr_tower_lazy_lock(ctx);
    if (IS_ALG_CTX(ctx))
    {
        truth_t o = _gr_tower_lazy_is_one(x, ctx);
        status = (o == T_TRUE) ? _gr_tower_lazy_zero(res, ctx) : ((o == T_FALSE) ? GR_DOMAIN : GR_UNABLE);
    }
    else
    {
        if (IS_REAL_CTX(ctx))
        {
            int sgn;
            status = _real_sign(&sgn, x, ctx);
            if (status == GR_SUCCESS && sgn <= 0)
                status = GR_DOMAIN;
        }
        if (status == GR_SUCCESS)
            status = _gr_tower_lazy_log(res, x, ctx);
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

int
gr_tower_lazy_set_other(gr_tower_lazy_elem_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = _gr_tower_lazy_set_other(res, x, x_ctx, ctx);
    if (status == GR_SUCCESS && _gr_tower_lazy_outermost(ctx))
        status = _gr_tower_lazy_check_member(res, ctx);
    if (status != GR_SUCCESS)
        GR_MUST_SUCCEED(_gr_tower_lazy_zero(res, ctx));
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/* in a real context, the real roots only */
int
gr_tower_lazy_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    int status;
    _gr_tower_lazy_lock(ctx);
    status = _gr_tower_lazy_poly_roots(roots, mult, poly, flags, ctx);
    if (status == GR_SUCCESS && IS_REAL_CTX(ctx))
    {
        slong i, j = 0, n = roots->length;
        for (i = 0; i < n && status == GR_SUCCESS; i++)
        {
            truth_t t = _gr_tower_lazy_is_real_exact(gr_vec_entry_ptr(roots, i, ctx), ctx);
            if (t == T_UNKNOWN)
                status = GR_UNABLE;
            else if (t == T_TRUE)
            {
                if (j != i)
                {
                    gr_swap(gr_vec_entry_ptr(roots, j, ctx), gr_vec_entry_ptr(roots, i, ctx), ctx);
                    fmpz_swap(mult->entries + j, mult->entries + i);
                }
                j++;
            }
        }
        if (status == GR_SUCCESS)
        {
            gr_vec_set_length(roots, j, ctx);
            fmpz_vec_set_length(mult, j);
            for (i = 0; i < j && status == GR_SUCCESS; i++)
                status = _gr_tower_lazy_realify(gr_vec_entry_ptr(roots, i, ctx), ctx);
        }
    }
    _gr_tower_lazy_unlock(ctx);
    return status;
}

/*
    Factorization over the field: linear factors x - r for the roots
    (complex fields), and in the real fields linear factors for the real
    roots and quadratic factors (x - r)(x - conj(r)) for the pairs of
    nonreal roots.
*/
int
gr_tower_lazy_poly_factor(gr_poly_t c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx)
{
    gr_ctx_t pctx;
    gr_vec_t roots;
    fmpz_vec_t rmult;
    slong i, total = 0;
    int status;

    if (poly->length == 0)
        return GR_DOMAIN;

    _gr_tower_lazy_lock(ctx);
    gr_ctx_init_gr_poly(pctx, ctx);
    gr_vec_init(roots, 0, ctx);
    fmpz_vec_init(rmult, 0);
    gr_vec_set_length(fac, 0, pctx);
    fmpz_vec_set_length(mult, 0);

    status = gr_poly_set_scalar(c, gr_poly_coeff_srcptr(poly, poly->length - 1, ctx), ctx);
    if (status == GR_SUCCESS)
        status = _gr_tower_lazy_poly_roots(roots, rmult, poly, 0, ctx);

    for (i = 0; i < roots->length && status == GR_SUCCESS; i++)
    {
        gr_ptr r = gr_vec_entry_ptr(roots, i, ctx);
        gr_poly_t h;
        truth_t real = T_TRUE;

        if (IS_REAL_CTX(ctx))
            real = _gr_tower_lazy_is_real_exact(r, ctx);
        if (real == T_UNKNOWN)
        {
            status = GR_UNABLE;
            break;
        }

        gr_poly_init(h, ctx);

        if (real == T_TRUE)
        {
            status |= gr_poly_set_coeff_si(h, 1, 1, ctx);
            status |= _gr_tower_lazy_neg(gr_poly_coeff_ptr(h, 0, ctx), r, ctx);
            total += fmpz_get_si(rmult->entries + i);
        }
        else
        {
            /* one factor for the root with positive imaginary part */
            acb_t z;
            slong prec;
            int sign = 0;

            acb_init(z);
            for (prec = 64; prec <= 65536 && sign == 0 && status == GR_SUCCESS; prec *= 2)
            {
                status = _gr_tower_lazy_get_acb_impl(z, r, prec, ctx);
                if (arb_is_positive(acb_imagref(z)))
                    sign = 1;
                else if (arb_is_negative(acb_imagref(z)))
                    sign = -1;
            }
            acb_clear(z);
            if (sign == 0)
                status = GR_UNABLE;

            if (status == GR_SUCCESS && sign > 0)
            {
                gr_tower_lazy_elem_struct rc, t;
                _gr_tower_lazy_init(&rc, ctx);
                _gr_tower_lazy_init(&t, ctx);
                status |= _gr_tower_lazy_conj(&rc, r, ctx);
                status |= gr_poly_set_coeff_si(h, 2, 1, ctx);
                status |= _gr_tower_lazy_add(&t, r, &rc, ctx);
                status |= _gr_tower_lazy_neg(&t, &t, ctx);
                status |= gr_poly_set_coeff_scalar(h, 1, &t, ctx);
                status |= _gr_tower_lazy_mul(&t, r, &rc, ctx);
                status |= gr_poly_set_coeff_scalar(h, 0, &t, ctx);
                _gr_tower_lazy_clear(&rc, ctx);
                _gr_tower_lazy_clear(&t, ctx);
                total += 2 * fmpz_get_si(rmult->entries + i);
            }
        }

        if (status == GR_SUCCESS && h->length > 0 && IS_REAL_CTX(ctx))
        {
            slong j;
            for (j = 0; j < h->length && status == GR_SUCCESS; j++)
                status = _gr_tower_lazy_realify(gr_poly_coeff_ptr(h, j, ctx), ctx);
        }

        if (status == GR_SUCCESS && h->length > 0)
        {
            gr_vec_set_length(fac, fac->length + 1, pctx);
            gr_swap(gr_vec_entry_ptr(fac, fac->length - 1, pctx), h, pctx);
            fmpz_vec_append(mult, rmult->entries + i);
        }
        gr_poly_clear(h, ctx);
    }

    if (status == GR_SUCCESS && total != poly->length - 1)
        status = GR_UNABLE;

    fmpz_vec_clear(rmult);
    gr_vec_clear(roots, ctx);
    gr_ctx_clear(pctx);
    _gr_tower_lazy_unlock(ctx);
    return status;
}

LOCK_PRED(_gr_tower_lazy_is_zero)
LOCK_PRED(_gr_tower_lazy_is_one)
static truth_t _gr_tower_lazy_equal_locked(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx) { LOCKED_T(truth_t, _gr_tower_lazy_equal(x, y, ctx)) }
static int _gr_tower_lazy_set_si_locked(gr_tower_lazy_elem_t res, slong c, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_set_si(res, c, ctx)) }
static int _gr_tower_lazy_set_ui_locked(gr_tower_lazy_elem_t res, ulong c, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_set_ui(res, c, ctx)) }
static int _gr_tower_lazy_set_fmpz_locked(gr_tower_lazy_elem_t res, const fmpz_t c, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_set_fmpz(res, c, ctx)) }
static int _gr_tower_lazy_set_fmpq_locked(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_set_fmpq(res, c, ctx)) }
static int _gr_tower_lazy_get_fmpz_locked(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_fmpz(res, x, ctx)) }
int gr_tower_lazy_get_si(slong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_si(res, x, ctx)) }
int gr_tower_lazy_get_ui(ulong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_ui(res, x, ctx)) }
static int _gr_tower_lazy_get_fmpq_locked(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_get_fmpq(res, x, ctx)) }

/* -------------------------------------------------------------------- */
/* rational fast paths (no towers, no lock; see LAZY_IS_RATIONAL)        */
/* -------------------------------------------------------------------- */

#define RAT(x) LAZY_IS_RATIONAL(x, LAZY(ctx))
#define LAZY_TRUTH(b) ((b) ? T_TRUE : T_FALSE)

void gr_tower_lazy_init(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* (the trivial tower is pinned: no reference count) */
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    x->F = L->trivial;
    x->level = 0;
    x->structure_version = x->F->T->structure_version;
    x->reduced_version = 0;
    x->shallow = 0;
    x->repr = LAZY_REPR_RATIONAL;
    fmpq_init(&x->elem.q);
    LAZY_POISON_FLAT(x);
    LAZY_POISON_DENSE(x);
}

void gr_tower_lazy_clear(gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (x->F == LAZY(ctx)->trivial)
        _gr_tower_lazy_clear_data(x);
    else
        _gr_tower_lazy_clear_locked(x, ctx);
}

/* (pure moves of the structs) */
void gr_tower_lazy_swap(gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx) { _gr_tower_lazy_swap(x, y, ctx); }
void gr_tower_lazy_set_shallow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx) { _gr_tower_lazy_set_shallow(res, x, ctx); }

#define RAT_BINARY(pub, name, fmpq_op, dense_op, nonzero_y) \
int pub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx) \
{ \
    if (RAT(x) && RAT(y) && RAT(res) && (!(nonzero_y) || !fmpq_is_zero(&y->elem.q))) \
    { \
        fmpq_op(&res->elem.q, &x->elem.q, &y->elem.q); \
        return GR_SUCCESS; \
    } \
    if ((LAZY_HAS_DENSE(x) || LAZY_HAS_DENSE(y)) && _gr_tower_lazy_dense_binary(res, x, y, dense_op, ctx)) \
        return GR_SUCCESS; \
    LOCKED_DENSE(name(res, x, y, ctx), res, x, y) \
}

RAT_BINARY(gr_tower_lazy_add, _gr_tower_lazy_add, fmpq_add, DENSE_ADD, 0)
RAT_BINARY(gr_tower_lazy_sub, _gr_tower_lazy_sub, fmpq_sub, DENSE_SUB, 0)
RAT_BINARY(gr_tower_lazy_mul, _gr_tower_lazy_mul, fmpq_mul, DENSE_MUL, 0)
RAT_BINARY(gr_tower_lazy_div, _gr_tower_lazy_div, fmpq_div, DENSE_DIV, 1)

int gr_tower_lazy_neg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x) && RAT(res))
    {
        fmpq_t a, r;
        _gr_tower_lazy_rational_view(a, x);
        fmpq_init(r);
        fmpq_neg(r, a);
        _gr_tower_lazy_rational_store(res, r);
        return GR_SUCCESS;
    }
    if (LAZY_HAS_DENSE(x))
    {
        fmpq_t m;
        int ok;
        fmpq_init(m);
        fmpz_set_si(fmpq_numref(m), -1);
        ok = _gr_tower_lazy_dense_scalar(res, x, m, DENSE_MUL, ctx);
        fmpq_clear(m);
        if (ok)
            return GR_SUCCESS;
    }
    LOCKED_DENSE(_gr_tower_lazy_neg(res, x, ctx), res, x, NULL)
}

int gr_tower_lazy_inv(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x) && RAT(res) && !fmpq_is_zero(&x->elem.q))
    {
        fmpq_t a, r;
        _gr_tower_lazy_rational_view(a, x);
        fmpq_init(r);
        fmpq_inv(r, a);
        _gr_tower_lazy_rational_store(res, r);
        return GR_SUCCESS;
    }
    if (LAZY_HAS_DENSE(x) && _gr_tower_lazy_tls_depth == 0 && !res->shallow && res->F == x->F)
    {
        gr_tower_lazy_dense_view_struct X, Y;
        fmpq_t one;
        fmpq_init(one);
        fmpq_one(one);
        _gr_tower_lazy_dense_view_dense(&Y, x);
        _gr_tower_lazy_dense_view_fmpq(&X, one);
        if (_gr_tower_lazy_dense_op(res, &X, &Y, DENSE_DIV, x->elem.dense.nf))
            return GR_SUCCESS;
    }
    LOCKED_DENSE(_gr_tower_lazy_inv(res, x, ctx), res, x, NULL)
}

int gr_tower_lazy_set(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (res == x)
        return GR_SUCCESS;
    if (RAT(x) && RAT(res))
    {
        fmpq_t a, r;
        _gr_tower_lazy_rational_view(a, x);
        fmpq_init(r);
        fmpq_set(r, a);
        _gr_tower_lazy_rational_store(res, r);
        return GR_SUCCESS;
    }
    if (LAZY_HAS_DENSE(x))
    {
        fmpq_t one;
        int ok;
        fmpq_init(one);
        fmpq_one(one);
        ok = _gr_tower_lazy_dense_scalar(res, x, one, DENSE_MUL, ctx);
        fmpq_clear(one);
        if (ok)
            return GR_SUCCESS;
    }
    LOCKED_DENSE(_gr_tower_lazy_set(res, x, ctx), res, x, NULL)
}

/* with a scalar c of type T, as the rational number q: r = a op q */
#define RAT_SCALAR(pub, name, T, body, setq, dense_op, nonzero_c) \
int pub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, T c, gr_ctx_t ctx) \
{ \
    if (RAT(x) && RAT(res) && !(nonzero_c)) \
    { \
        fmpq_t a, r; \
        _gr_tower_lazy_rational_view(a, x); \
        fmpq_init(r); \
        body; \
        _gr_tower_lazy_rational_store(res, r); \
        return GR_SUCCESS; \
    } \
    if (LAZY_HAS_DENSE(x) && !(nonzero_c)) \
    { \
        fmpq_t q; \
        int ok; \
        fmpq_init(q); \
        setq; \
        ok = _gr_tower_lazy_dense_scalar(res, x, q, dense_op, ctx); \
        fmpq_clear(q); \
        if (ok) \
            return GR_SUCCESS; \
    } \
    LOCKED_DENSE(name(res, x, c, ctx), res, x, NULL) \
}

RAT_SCALAR(gr_tower_lazy_add_si, _gr_tower_lazy_add_si, slong, fmpq_add_si(r, a, c), fmpz_set_si(fmpq_numref(q), c), DENSE_ADD, 0)
RAT_SCALAR(gr_tower_lazy_add_ui, _gr_tower_lazy_add_ui, ulong, fmpq_add_ui(r, a, c), fmpz_set_ui(fmpq_numref(q), c), DENSE_ADD, 0)
RAT_SCALAR(gr_tower_lazy_add_fmpz, _gr_tower_lazy_add_fmpz, const fmpz_t, fmpq_add_fmpz(r, a, c), fmpz_set(fmpq_numref(q), c), DENSE_ADD, 0)
RAT_SCALAR(gr_tower_lazy_add_fmpq, _gr_tower_lazy_add_fmpq, const fmpq_t, fmpq_add(r, a, c), fmpq_set(q, c), DENSE_ADD, 0)
RAT_SCALAR(gr_tower_lazy_sub_si, _gr_tower_lazy_sub_si, slong, fmpq_sub_si(r, a, c), fmpz_set_si(fmpq_numref(q), c), DENSE_SUB, 0)
RAT_SCALAR(gr_tower_lazy_sub_ui, _gr_tower_lazy_sub_ui, ulong, fmpq_sub_ui(r, a, c), fmpz_set_ui(fmpq_numref(q), c), DENSE_SUB, 0)
RAT_SCALAR(gr_tower_lazy_sub_fmpz, _gr_tower_lazy_sub_fmpz, const fmpz_t, fmpq_sub_fmpz(r, a, c), fmpz_set(fmpq_numref(q), c), DENSE_SUB, 0)
RAT_SCALAR(gr_tower_lazy_sub_fmpq, _gr_tower_lazy_sub_fmpq, const fmpq_t, fmpq_sub(r, a, c), fmpq_set(q, c), DENSE_SUB, 0)
RAT_SCALAR(gr_tower_lazy_mul_si, _gr_tower_lazy_mul_si, slong, fmpq_mul_si(r, a, c), fmpz_set_si(fmpq_numref(q), c), DENSE_MUL, 0)
RAT_SCALAR(gr_tower_lazy_mul_ui, _gr_tower_lazy_mul_ui, ulong, fmpq_mul_ui(r, a, c), fmpz_set_ui(fmpq_numref(q), c), DENSE_MUL, 0)
RAT_SCALAR(gr_tower_lazy_mul_fmpz, _gr_tower_lazy_mul_fmpz, const fmpz_t, fmpq_mul_fmpz(r, a, c), fmpz_set(fmpq_numref(q), c), DENSE_MUL, 0)
RAT_SCALAR(gr_tower_lazy_mul_fmpq, _gr_tower_lazy_mul_fmpq, const fmpq_t, fmpq_mul(r, a, c), fmpq_set(q, c), DENSE_MUL, 0)
RAT_SCALAR(gr_tower_lazy_div_si, _gr_tower_lazy_div_si, slong, { fmpz_t t; fmpz_init_set_si(t, c); fmpq_div_fmpz(r, a, t); fmpz_clear(t); }, fmpz_set_si(fmpq_numref(q), c), DENSE_DIV, c == 0)
RAT_SCALAR(gr_tower_lazy_div_ui, _gr_tower_lazy_div_ui, ulong, { fmpz_t t; fmpz_init_set_ui(t, c); fmpq_div_fmpz(r, a, t); fmpz_clear(t); }, fmpz_set_ui(fmpq_numref(q), c), DENSE_DIV, c == 0)
RAT_SCALAR(gr_tower_lazy_div_fmpz, _gr_tower_lazy_div_fmpz, const fmpz_t, fmpq_div_fmpz(r, a, c), fmpz_set(fmpq_numref(q), c), DENSE_DIV, fmpz_is_zero(c))
RAT_SCALAR(gr_tower_lazy_div_fmpq, _gr_tower_lazy_div_fmpq, const fmpq_t, fmpq_div(r, a, c), fmpq_set(q, c), DENSE_DIV, fmpq_is_zero(c))

#define RAT_SET(pub, name, T, body) \
int pub(gr_tower_lazy_elem_t res, T c, gr_ctx_t ctx) \
{ \
    if (RAT(res)) \
    { \
        fmpq_t r; \
        fmpq_init(r); \
        body; \
        _gr_tower_lazy_rational_store(res, r); \
        return GR_SUCCESS; \
    } \
    return name##_locked(res, c, ctx); \
}

RAT_SET(gr_tower_lazy_set_si, _gr_tower_lazy_set_si, slong, fmpq_set_si(r, c, 1))
RAT_SET(gr_tower_lazy_set_ui, _gr_tower_lazy_set_ui, ulong, fmpq_set_ui(r, c, 1))
RAT_SET(gr_tower_lazy_set_fmpz, _gr_tower_lazy_set_fmpz, const fmpz_t, fmpz_set(fmpq_numref(r), c))
RAT_SET(gr_tower_lazy_set_fmpq, _gr_tower_lazy_set_fmpq, const fmpq_t, fmpq_set(r, c))

int gr_tower_lazy_zero(gr_tower_lazy_elem_t res, gr_ctx_t ctx) { return gr_tower_lazy_set_si(res, 0, ctx); }
int gr_tower_lazy_one(gr_tower_lazy_elem_t res, gr_ctx_t ctx) { return gr_tower_lazy_set_si(res, 1, ctx); }

/* (the dense forms are canonical: structural comparisons) */
static int
_dense_view_equal(const gr_tower_lazy_dense_view_struct * X, const gr_tower_lazy_dense_view_struct * Y)
{
    return fmpq_poly_equal(X->p, Y->p);
}

truth_t gr_tower_lazy_is_zero(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x))
        return LAZY_TRUTH(fmpq_is_zero(&x->elem.q));
    if (LAZY_HAS_DENSE(x))
        return LAZY_TRUTH(x->elem.dense.poly.length == 0);
    return _gr_tower_lazy_is_zero_locked(x, ctx);
}

truth_t gr_tower_lazy_is_one(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x))
    {
        fmpq_t a;
        _gr_tower_lazy_rational_view(a, x);
        return LAZY_TRUTH(fmpq_is_one(a));
    }
    if (LAZY_HAS_DENSE(x))
    {
        return LAZY_TRUTH(fmpq_poly_is_one(&x->elem.dense.poly));
    }
    return _gr_tower_lazy_is_one_locked(x, ctx);
}

truth_t gr_tower_lazy_equal(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx)
{
    if (RAT(x) && RAT(y))
    {
        fmpq_t a, b;
        _gr_tower_lazy_rational_view(a, x);
        _gr_tower_lazy_rational_view(b, y);
        return LAZY_TRUTH(fmpq_equal(a, b));
    }
    if (LAZY_HAS_DENSE(x) || LAZY_HAS_DENSE(y))
    {
        const gr_tower_lazy_elem_struct * D = LAZY_HAS_DENSE(x) ? x : y;
        gr_tower_lazy_dense_view_struct X, Y;
        if (_gr_tower_lazy_dense_view(&X, x, D->elem.dense.nf, D->F, ctx) && _gr_tower_lazy_dense_view(&Y, y, D->elem.dense.nf, D->F, ctx))
            return LAZY_TRUTH(_dense_view_equal(&X, &Y));
    }
    return _gr_tower_lazy_equal_locked(x, y, ctx);
}

int gr_tower_lazy_get_fmpq(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x))
    {
        fmpq_t a;
        _gr_tower_lazy_rational_view(a, x);
        fmpq_set(res, a);
        return GR_SUCCESS;
    }
    if (LAZY_HAS_DENSE(x))
    {
        if (x->elem.dense.poly.length >= 2)
            return GR_DOMAIN;
        fmpq_poly_get_coeff_fmpq(res, &x->elem.dense.poly, 0);
        return GR_SUCCESS;
    }
    return _gr_tower_lazy_get_fmpq_locked(res, x, ctx);
}

int gr_tower_lazy_get_fmpz(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    if (RAT(x))
    {
        fmpq_t a;
        _gr_tower_lazy_rational_view(a, x);
        if (!fmpz_is_one(fmpq_denref(a)))
            return GR_DOMAIN;
        fmpz_set(res, fmpq_numref(a));
        return GR_SUCCESS;
    }
    return _gr_tower_lazy_get_fmpz_locked(res, x, ctx);
}

#undef RAT

LOCK_UNARY(gr_tower_lazy_conj, _gr_tower_lazy_conj)
LOCK_UNARY(gr_tower_lazy_abs, _gr_tower_lazy_abs)
int gr_tower_lazy_gens(gr_vec_t vec, gr_ctx_t ctx) { LOCKED(_gr_tower_lazy_gens(vec, ctx)) }

/*
    Reads "root(m(name), approximation)": the polynomial m in the
    variable name over the field, in the terminals given, and the root
    nearest to the approximation "re", "re + im*i", "re - im*i" or "im*i".
*/
/* the unit in the last place of a decimal string such as "-0.00990099",
   "1.5e-10", "3.14 - 2.7*i" (the largest among its parts): 10^(exponent -
   fraction digits) */
/* -------------------------------------------------------------------- */
/* method table                                                          */
/* -------------------------------------------------------------------- */

int _gr_tower_lazy_methods_initialized = 0;

gr_static_method_table _gr_tower_lazy_methods;

gr_method_tab_input _gr_tower_lazy_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_tower_lazy_ctx_write},
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_tower_lazy_ctx_clear},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) _gr_tower_lazy_ctx_is_field},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) _gr_tower_lazy_ctx_is_field},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN, (gr_funcptr) _gr_tower_lazy_ctx_is_field},
    {GR_METHOD_CTX_IS_RATIONAL_VECTOR_SPACE, (gr_funcptr) _gr_tower_lazy_ctx_is_rational_vector_space},
    {GR_METHOD_CTX_IS_REAL_VECTOR_SPACE, (gr_funcptr) _gr_tower_lazy_ctx_is_real_vector_space},
    {GR_METHOD_CTX_IS_COMPLEX_VECTOR_SPACE, (gr_funcptr) _gr_tower_lazy_ctx_is_complex_vector_space},
    {GR_METHOD_CTX_IS_THREADSAFE,       (gr_funcptr) _gr_tower_lazy_ctx_is_threadsafe},
    {GR_METHOD_CTX_IS_FINITE,           (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC, (gr_funcptr) _gr_tower_lazy_ctx_is_finite_characteristic},
    {GR_METHOD_CTX_IS_EXACT,            (gr_funcptr) _gr_tower_lazy_ctx_is_exact},
    {GR_METHOD_CTX_IS_CANONICAL,        (gr_funcptr) _gr_tower_lazy_ctx_is_canonical},
    {GR_METHOD_CTX_BASE,                (gr_funcptr) _gr_tower_lazy_ctx_base},

    {GR_METHOD_INIT,            (gr_funcptr) gr_tower_lazy_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) gr_tower_lazy_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) gr_tower_lazy_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) gr_tower_lazy_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) gr_tower_lazy_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) gr_tower_lazy_write},
    {GR_METHOD_ZERO,            (gr_funcptr) gr_tower_lazy_zero},
    {GR_METHOD_ONE,             (gr_funcptr) gr_tower_lazy_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) gr_tower_lazy_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) gr_tower_lazy_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) gr_tower_lazy_equal},
    {GR_METHOD_SET,             (gr_funcptr) gr_tower_lazy_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) gr_tower_lazy_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) gr_tower_lazy_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) gr_tower_lazy_set_fmpz},
    {GR_METHOD_SET_FMPQ,        (gr_funcptr) gr_tower_lazy_set_fmpq},
    {GR_METHOD_SET_OTHER,       (gr_funcptr) gr_tower_lazy_set_other},
    {GR_METHOD_NEG,             (gr_funcptr) gr_tower_lazy_neg},
    {GR_METHOD_ADD,             (gr_funcptr) gr_tower_lazy_add},
    {GR_METHOD_SUB,             (gr_funcptr) gr_tower_lazy_sub},
    {GR_METHOD_MUL,             (gr_funcptr) gr_tower_lazy_mul},
    {GR_METHOD_INV,             (gr_funcptr) gr_tower_lazy_inv},
    {GR_METHOD_ADD_SI, (gr_funcptr) gr_tower_lazy_add_si},
    {GR_METHOD_ADD_UI, (gr_funcptr) gr_tower_lazy_add_ui},
    {GR_METHOD_ADD_FMPZ, (gr_funcptr) gr_tower_lazy_add_fmpz},
    {GR_METHOD_ADD_FMPQ, (gr_funcptr) gr_tower_lazy_add_fmpq},
    {GR_METHOD_SUB_SI, (gr_funcptr) gr_tower_lazy_sub_si},
    {GR_METHOD_SUB_UI, (gr_funcptr) gr_tower_lazy_sub_ui},
    {GR_METHOD_SUB_FMPZ, (gr_funcptr) gr_tower_lazy_sub_fmpz},
    {GR_METHOD_SUB_FMPQ, (gr_funcptr) gr_tower_lazy_sub_fmpq},
    {GR_METHOD_MUL_SI, (gr_funcptr) gr_tower_lazy_mul_si},
    {GR_METHOD_MUL_UI, (gr_funcptr) gr_tower_lazy_mul_ui},
    {GR_METHOD_MUL_FMPZ, (gr_funcptr) gr_tower_lazy_mul_fmpz},
    {GR_METHOD_MUL_FMPQ, (gr_funcptr) gr_tower_lazy_mul_fmpq},
    {GR_METHOD_DIV_SI, (gr_funcptr) gr_tower_lazy_div_si},
    {GR_METHOD_DIV_UI, (gr_funcptr) gr_tower_lazy_div_ui},
    {GR_METHOD_DIV_FMPZ, (gr_funcptr) gr_tower_lazy_div_fmpz},
    {GR_METHOD_DIV_FMPQ, (gr_funcptr) gr_tower_lazy_div_fmpq},
    {GR_METHOD_DIV,             (gr_funcptr) gr_tower_lazy_div},
    {GR_METHOD_SQRT,            (gr_funcptr) gr_tower_lazy_sqrt},
    {GR_METHOD_POW,             (gr_funcptr) gr_tower_lazy_pow},
    {GR_METHOD_SET_STR,         (gr_funcptr) gr_tower_lazy_set_str},
    {GR_METHOD_ROOT_UI,         (gr_funcptr) gr_tower_lazy_root_ui},
    {GR_METHOD_GET_SI,          (gr_funcptr) gr_tower_lazy_get_si},
    {GR_METHOD_GET_UI,          (gr_funcptr) gr_tower_lazy_get_ui},
    {GR_METHOD_GET_FMPZ,        (gr_funcptr) gr_tower_lazy_get_fmpz},
    {GR_METHOD_GET_FMPQ,        (gr_funcptr) gr_tower_lazy_get_fmpq},
    {GR_METHOD_PI,              (gr_funcptr) gr_tower_lazy_pi},
    {GR_METHOD_I,               (gr_funcptr) gr_tower_lazy_i},
    {GR_METHOD_CONJ,            (gr_funcptr) gr_tower_lazy_conj},
    {GR_METHOD_GET_D,           (gr_funcptr) gr_tower_lazy_get_d},
    {GR_METHOD_RE,              (gr_funcptr) gr_tower_lazy_re},
    {GR_METHOD_IM,              (gr_funcptr) gr_tower_lazy_im},
    {GR_METHOD_ARG,             (gr_funcptr) gr_tower_lazy_arg},
    {GR_METHOD_SGN,             (gr_funcptr) gr_tower_lazy_sgn},
    {GR_METHOD_CSGN,            (gr_funcptr) gr_tower_lazy_csgn},
    {GR_METHOD_IS_REAL,         (gr_funcptr) gr_tower_lazy_is_real},
    {GR_METHOD_CMP,             (gr_funcptr) gr_tower_lazy_cmp},
    {GR_METHOD_CMPABS,          (gr_funcptr) gr_tower_lazy_cmpabs},
    {GR_METHOD_FLOOR,           (gr_funcptr) gr_tower_lazy_floor},
    {GR_METHOD_CEIL,            (gr_funcptr) gr_tower_lazy_ceil},
    {GR_METHOD_TRUNC,           (gr_funcptr) gr_tower_lazy_trunc},
    {GR_METHOD_NINT,            (gr_funcptr) gr_tower_lazy_nint},
    {GR_METHOD_SIN,             (gr_funcptr) gr_tower_lazy_sin},
    {GR_METHOD_COS,             (gr_funcptr) gr_tower_lazy_cos},
    {GR_METHOD_TAN,             (gr_funcptr) gr_tower_lazy_tan},
    {GR_METHOD_SINH,            (gr_funcptr) gr_tower_lazy_sinh},
    {GR_METHOD_COSH,            (gr_funcptr) gr_tower_lazy_cosh},
    {GR_METHOD_TANH,            (gr_funcptr) gr_tower_lazy_tanh},
    {GR_METHOD_ASIN,            (gr_funcptr) gr_tower_lazy_asin},
    {GR_METHOD_ACOS,            (gr_funcptr) gr_tower_lazy_acos},
    {GR_METHOD_ATAN,            (gr_funcptr) gr_tower_lazy_atan},
    {GR_METHOD_GAMMA,           (gr_funcptr) gr_tower_lazy_gamma},
    {GR_METHOD_RGAMMA,          (gr_funcptr) gr_tower_lazy_rgamma},
    {GR_METHOD_BETA,            (gr_funcptr) gr_tower_lazy_beta},
    {GR_METHOD_DIGAMMA,         (gr_funcptr) gr_tower_lazy_digamma},
    {GR_METHOD_POLYGAMMA,       (gr_funcptr) gr_tower_lazy_polygamma},
    {GR_METHOD_ERF,             (gr_funcptr) gr_tower_lazy_erf},
    {GR_METHOD_ERFC,            (gr_funcptr) gr_tower_lazy_erfc},
    {GR_METHOD_ERFI,            (gr_funcptr) gr_tower_lazy_erfi},
    {GR_METHOD_LAMBERTW,        (gr_funcptr) gr_tower_lazy_lambertw},
    {GR_METHOD_LAMBERTW_FMPZ,   (gr_funcptr) gr_tower_lazy_lambertw_fmpz},
    {GR_METHOD_ZETA,            (gr_funcptr) gr_tower_lazy_zeta},
    {GR_METHOD_HURWITZ_ZETA,    (gr_funcptr) gr_tower_lazy_hurwitz_zeta},
    {GR_METHOD_POLYLOG,         (gr_funcptr) gr_tower_lazy_polylog},
    {GR_METHOD_DILOG,           (gr_funcptr) gr_tower_lazy_dilog},
    {GR_METHOD_ELLIPTIC_K,      (gr_funcptr) gr_tower_lazy_elliptic_k},
    {GR_METHOD_ELLIPTIC_E,      (gr_funcptr) gr_tower_lazy_elliptic_e},
    {GR_METHOD_EULER,           (gr_funcptr) gr_tower_lazy_euler},
    {GR_METHOD_CATALAN,         (gr_funcptr) gr_tower_lazy_catalan},
    {GR_METHOD_ASINH,           (gr_funcptr) gr_tower_lazy_asinh},
    {GR_METHOD_ACOSH,           (gr_funcptr) gr_tower_lazy_acosh},
    {GR_METHOD_ATANH,           (gr_funcptr) gr_tower_lazy_atanh},
    {GR_METHOD_ABS,             (gr_funcptr) gr_tower_lazy_abs},
    {GR_METHOD_EXP,             (gr_funcptr) gr_tower_lazy_exp},
    {GR_METHOD_LOG,             (gr_funcptr) gr_tower_lazy_log},
    {GR_METHOD_MAT_DET,         (gr_funcptr) gr_tower_lazy_mat_det},
    {GR_METHOD_POLY_ROOTS,      (gr_funcptr) gr_tower_lazy_poly_roots},
    {GR_METHOD_POLY_FACTOR,     (gr_funcptr) gr_tower_lazy_poly_factor},
    {GR_METHOD_GENS,            (gr_funcptr) gr_tower_lazy_gens},
    /* polynomial gcds, xgcds and resultants by the subresultant PRS: the
       Euclidean remainder sequence over a function field with algebraic
       generators has exponential coefficient growth */
    {GR_METHOD_POLY_GCD,        (gr_funcptr) _gr_poly_gcd_subresultant},
    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) gr_tower_lazy_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) gr_tower_lazy_mat_mul},
    {GR_METHOD_MAT_NONSINGULAR_SOLVE, (gr_funcptr) gr_tower_lazy_mat_nonsingular_solve},
    {GR_METHOD_POLY_XGCD,       (gr_funcptr) _gr_poly_xgcd_subresultant},
    {GR_METHOD_POLY_RESULTANT,  (gr_funcptr) _gr_poly_resultant_subresultant},
    {GR_METHOD_GET_FEXPR,       (gr_funcptr) gr_tower_lazy_get_fexpr},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_tower_lazy(gr_ctx_t ctx, gr_ctx_t base, int flags)
{
    gr_tower_lazy_ctx_struct * L;

    L = flint_calloc(1, sizeof(gr_tower_lazy_ctx_struct));
    L->base = base;
    L->merge_flags = flags & GR_TOWER_MERGE_EXPRESS;
    if (flags & GR_TOWER_LAZY_REAL)
        L->merge_flags |= GR_TOWER_MERGE_REAL_FIRST;
    L->refcount = 1;
    {
        slong k;
        for (k = 0; k < GR_TOWER_OPT_NUM_OPTIONS; k++)
            L->options[k] = gr_tower_default_options[k];
    }
    L->depth = 0;
    L->num_alg_names = L->num_trans_names = 0;
    L->conj_depth = 0;
#if FLINT_USES_PTHREAD
    {
        pthread_mutexattr_t attr;
        pthread_mutexattr_init(&attr);
        pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_RECURSIVE);
        pthread_mutex_init(&L->mutex, &attr);
        pthread_mutexattr_destroy(&attr);
    }
#endif

    ctx->which_ring = GR_CTX_GR_TOWER_LAZY;
    ctx->sizeof_elem = sizeof(gr_tower_lazy_elem_struct);
    ctx->size_limit = WORD_MAX;
    VIEW(ctx)->L = L;
    VIEW(ctx)->field_flags = flags & (GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC);

    L->trivial = _gr_tower_lazy_new_tower(ctx);
    L->trivial->gc |= GR_TOWER_GC_PINNED;

    ctx->methods = _gr_tower_lazy_methods;

    if (!_gr_tower_lazy_methods_initialized)
    {
        gr_method_tab_init(_gr_tower_lazy_methods, _gr_tower_lazy_methods_input);
        _gr_tower_lazy_methods_initialized = 1;
    }
}

/* a context sharing the state (and the elements) of parent, restricted
   further by field_flags */
void
gr_ctx_init_tower_lazy_view(gr_ctx_t ctx, gr_ctx_t parent, int field_flags)
{
    gr_tower_lazy_ctx_struct * L = LAZY(parent);

    _gr_tower_lazy_lock(parent);
    L->refcount++;
    /* (real views prefer real generators when towers are merged) */
    if (field_flags & GR_TOWER_LAZY_REAL)
        L->merge_flags |= GR_TOWER_MERGE_REAL_FIRST;
    _gr_tower_lazy_unlock(parent);

    ctx->which_ring = GR_CTX_GR_TOWER_LAZY;
    ctx->sizeof_elem = sizeof(gr_tower_lazy_elem_struct);
    ctx->size_limit = WORD_MAX;
    VIEW(ctx)->L = L;
    VIEW(ctx)->field_flags = VIEW(parent)->field_flags | (field_flags & (GR_TOWER_LAZY_REAL | GR_TOWER_LAZY_ALGEBRAIC));
    ctx->methods = _gr_tower_lazy_methods;
}

/* whether the two contexts share their elements */
int
gr_tower_lazy_ctx_same_state(gr_ctx_t ctx1, gr_ctx_t ctx2)
{
    return ctx1->which_ring == GR_CTX_GR_TOWER_LAZY && ctx2->which_ring == GR_CTX_GR_TOWER_LAZY && LAZY(ctx1) == LAZY(ctx2);
}

int
_gr_tower_lazy_get_qqbar_impl(qqbar_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    return _gr_tower_lazy_get_qqbar(res, (gr_tower_lazy_elem_struct *) x, ctx);
}


int
_gr_tower_lazy_get_acb_impl(acb_t res, const gr_tower_lazy_elem_t x, slong prec, gr_ctx_t ctx)
{
    gr_tower_lazy_elem_struct * e = (gr_tower_lazy_elem_struct *) x;
    e = _gr_tower_lazy_flat_view(e);
    return gr_tower_flat_get_acb(res, &e->elem.flat.data, prec, e->F);
}

gr_tower_struct *
_gr_tower_lazy_get_tower_impl(slong * level, const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    const gr_tower_lazy_elem_struct * e = x;
    if (level != NULL)
        *level = e->level;
    return e->F->T;
}

const fmpz_mpoly_q_struct *
_gr_tower_lazy_get_data_impl(const gr_tower_lazy_elem_t x, gr_ctx_t ctx)
{
    /* The caller owns x, so it is converted to the flat representation
       in place; the pointer stays valid until the next operation on x. */
    gr_tower_lazy_elem_struct * e = (gr_tower_lazy_elem_struct *) x;
    _gr_tower_lazy_make_flat(e);
    _gr_tower_lazy_update(e);
    return &e->elem.flat.data;
}

/* number of base-field coefficients of a nested element */
static slong
_nested_size(gr_srcptr x, gr_ctx_t ctx)
{
    if (ctx->which_ring == GR_CTX_GR_POLY_QUOTIENT)
    {
        const gr_poly_struct * p = (const gr_poly_struct *) x;
        gr_ctx_struct * cctx = gr_poly_quotient_ctx_base(ctx);
        slong j, s = 0;
        for (j = 0; j < p->length; j++)
            s += _nested_size(gr_poly_coeff_ptr((gr_poly_struct *) p, j, cctx), cctx);
        return s;
    }
    return 1;
}

void
_gr_tower_lazy_ctx_stats_impl(gr_ctx_t ctx)
{
    gr_tower_lazy_ctx_struct * L = LAZY(ctx);
    slong i, maxdeg = 0, totlen = 0;
    for (i = 0; i < L->num_towers; i++)
    {
        maxdeg = FLINT_MAX(maxdeg, gr_tower_degree(L->towers[i]->T));
        totlen += gr_tower_length(L->towers[i]->T);
    }
    flint_printf("towers: %wd (max degree %wd, total steps %wd), aliases: %wd, defs: %wu, collected: %wd\n", L->num_towers, maxdeg, totlen, L->num_aliases, L->next_def_id, L->gc_collected);

    /* sizes of the moduli of the largest towers (VERBOSE >= 2) */
    if (L->options[GR_TOWER_OPT_VERBOSE] >= 2)
    {
        for (i = 0; i < L->num_towers; i++)
        {
            gr_tower_flat_struct * F = L->towers[i];
            gr_tower_struct * T = F->T;
            slong k;

            if (gr_tower_length(T) < 2)
                continue;

            flint_printf("tower %wd: %wd algebraic steps, flat ideal of length %wd in %wd variables\n",
                i, gr_tower_length(T), F->ideal_len, F->mctx ? (slong) fmpz_mpoly_ctx_nvars(F->mctx) : 0);

            for (k = 0; k < F->ideal_len; k++)
            {
                slong terms, bits, nested_terms;
                fmpz_mpoly_struct * f = F->ideal[k];
                terms = (f == NULL) ? 0 : f->length;
                bits = terms ? _fmpz_vec_max_bits(f->coeffs, terms) : 0;
                if (bits < 0) bits = -bits;
                nested_terms = 0;
                {
                    /* the nested modulus: number of base coefficients */
                    const gr_poly_struct * m = gr_poly_quotient_ctx_modulus(GR_TOWER_STEP_CTX(T, k));
                    gr_ctx_struct * cctx = gr_poly_quotient_ctx_base(GR_TOWER_STEP_CTX(T, k));
                    slong j;
                    for (j = 0; j < m->length; j++)
                        nested_terms += _nested_size(gr_poly_coeff_ptr((gr_poly_struct *) m, j, cctx), cctx);
                }
                if (f == NULL)
                    flint_printf("  modulus %wd: flat conversion deferred, %wd nested base coefficients\n", k + 1, nested_terms);
                else
                    flint_printf("  modulus %wd: %wd flat terms, max %wd coefficient bits (~%wd KB), %wd nested base coefficients\n",
                        k + 1, terms, bits, (terms * (bits / 8 + 8 * (1 + (fmpz_mpoly_ctx_nvars(F->mctx) * f->bits) / FLINT_BITS))) / 1024, nested_terms);
            }
        }
    }
}
