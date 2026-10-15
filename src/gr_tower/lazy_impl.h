/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Lazy field of algebraic numbers: each element carries a pointer to a
    tower (shared between elements) and a level in that tower, and is
    stored in one of three representations (see gr_tower_lazy_elem_struct below):
    a rational number, the flat representation (fmpz_mpoly_q in the
    generators of the tower), or a dense form in a number field of one
    generator. Towers are created and extended on demand; arithmetic
    between elements of different towers goes through a tower containing
    the definitions of both, found by definition ids or constructed.

    Towers only grow (new steps appended) or are refined in place, so an
    element referring to (tower, level) never becomes invalid. When a
    tower grows beyond the capacity of its polynomial context, elements
    created with the old context are converted on use.

    Every operation of the lazy field holds the mutex of its shared
    state (gr_tower_lazy_ctx_struct.mutex; the dense and rational fast paths
    read elements without it, see the comments on the forms).

    This is the internal header of the implementation (the files
    src/gr_tower/lazy*.c): the structures, the macros and the functions
    shared between the files, which are

    lazy.c          contexts, elements, the locking layer, method table
    lazy_arith.c    arithmetic, comparisons, inverses
    lazy_merge.c    common towers: merging, primitive elements
    lazy_radical.c  radicals
    lazy_explog.c   exponentials, logarithms, transcendental generators
    lazy_conj.c     conversions, conjugation, transfer between contexts
    lazy_roots.c    roots of polynomials
    lazy_real.c     real forms
    lazy_dense.c    dense products and dense forms
    lazy_io.c       printing, expressions, parsing
*/

#ifndef GR_TOWER_LAZY_IMPL_H
#define GR_TOWER_LAZY_IMPL_H

#include <string.h>
#include <math.h>
#include "fmpz_mpoly_factor.h"
#include "fmpq.h"
#include "fmpzi.h"
#include "fmpz_factor.h"
#include "fmpz_poly.h"
#include "fmpz_poly_factor.h"
#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "fmpz_mpoly.h"
#include "fmpz_mpoly_q.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_special.h"
#include "fexpr.h"
#include "fexpr_builtin.h"
#include "gr_generic.h"
#include "gr_vec.h"
#include "gr_poly.h"
#include "gr_mat.h"
#include "fmpq_mat.h"
#include "fmpq_vec.h"
#include "acb_poly.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "arb.h"
#include "mpoly.h"
#include "gr_tower.h"
#include "gr_tower_lazy.h"
#include <stdlib.h>
#include "gr_tower/impl.h"
#if FLINT_USES_PTHREAD
#include <pthread.h>
#endif



/*
    An element of a lazy field: a value of the tower F (via F->T) in one
    of three representations (repr):

    - LAZY_REPR_RATIONAL: a rational number q; F is the trivial tower
      (no generators), level 0.
    - LAZY_REPR_FLAT: the flat representation data (fmpz_mpoly_q) in the
      polynomial context mctx of F (the current one, or a superseded one
      of F: elements are brought up to date on use), with the prefix
      length level.
    - LAZY_REPR_DENSE: an element of a number field of one generator a
      (level 1 of F, with a a proven algebraic generator of degree d
      with a monic integral modulus, the descriptor nf of F):
      (c[0] + c[1] a + ... + c[d-1] a^(d-1)) / den with den > 0 and the
      content of c coprime to it.

    Elements in the rational and the dense representations are read
    without the lock by the fast paths of the method table, so they are
    never changed in place while they may be operands: the locked
    implementations read every operand through _gr_tower_lazy_flat_view (a
    temporary flat copy for those two representations, the element
    itself for a flat one, brought up to date in place), and write only
    their result (its representation is whatever they produce: flat by
    default, rational for a constant, dense for a result of the dense
    arithmetic).
*/
/* (the element type gr_tower_lazy_elem_struct is public, read-only: gr_tower_lazy.h) */
#define LAZY_REPR_FLAT GR_TOWER_LAZY_REPR_FLAT
#define LAZY_REPR_RATIONAL GR_TOWER_LAZY_REPR_RATIONAL
#define LAZY_REPR_DENSE GR_TOWER_LAZY_REPR_DENSE

#ifdef LAZY_UNION_DEBUG
/* (debugging the representations: the fields of the representations an
   element is not in are poisoned, so that a read of them crashes) */
#define LAZY_POISON_FLAT(x) do { (x)->elem.flat.mctx = (fmpz_mpoly_ctx_struct *) 8; fmpz_mpoly_q_numref(&(x)->elem.flat.data)->coeffs = (fmpz *) 8; fmpz_mpoly_q_numref(&(x)->elem.flat.data)->exps = (ulong *) 8; fmpz_mpoly_q_numref(&(x)->elem.flat.data)->length = 1; fmpz_mpoly_q_denref(&(x)->elem.flat.data)->coeffs = (fmpz *) 8; fmpz_mpoly_q_denref(&(x)->elem.flat.data)->exps = (ulong *) 8; fmpz_mpoly_q_denref(&(x)->elem.flat.data)->length = 1; } while (0)
#define LAZY_POISON_DENSE(x) do { (x)->elem.dense.poly.coeffs = (fmpz *) 8; (x)->elem.dense.poly.alloc = 0; (x)->elem.dense.poly.length = 1; (x)->elem.dense.nf = (const _gr_tower_dense_field_struct *) 8; } while (0)
#define LAZY_POISON_RATIONAL(x) do { fmpq_numref(&(x)->elem.q)[0] = 0; fmpq_denref(&(x)->elem.q)[0] = 0; } while (0)
#else
#define LAZY_POISON_FLAT(x) do { } while (0)
#define LAZY_POISON_DENSE(x) do { } while (0)
#define LAZY_POISON_RATIONAL(x) do { } while (0)
#endif

/* (the representation of an operand may be changed by the lock-holding
   thread -- the flat representation into the dense one, see
   _gr_tower_lazy_try_dense -- while other threads read it without the lock: the
   data is written before the representation is published) */
#if defined(__GNUC__)
#define LAZY_REPR(x) __atomic_load_n(&(x)->repr, __ATOMIC_ACQUIRE)
#define LAZY_REPR_SET(x, v) __atomic_store_n(&(x)->repr, (v), __ATOMIC_RELEASE)
#else
#define LAZY_REPR(x) ((x)->repr)
#define LAZY_REPR_SET(x, v) ((x)->repr = (v))
#endif
#define LAZY_HAS_DENSE(x) (LAZY_REPR(x) == LAZY_REPR_DENSE)
#define LAZY_IS_FLAT(x) (LAZY_REPR(x) == LAZY_REPR_FLAT)

/* the depth of the lazy field operations in progress in this thread
   (dense forms are made and used only outside of them: the
   implementations read the flat data of their elements) */
extern FLINT_TLS_PREFIX slong _gr_tower_lazy_tls_depth;

/* A definition found expressible in a tower (rather than adjoined). */
typedef struct
{
    gr_tower_flat_struct * F;
    ulong def_id;
    fmpz_mpoly_ctx_struct * mctx;
    fmpz_mpoly_q_struct data;
}
gr_tower_lazy_alias_struct;

typedef struct
{
    qqbar_struct x;
    gr_tower_flat_struct * F;
    slong gid;
}
gr_tower_lazy_qqbar_entry_struct;

struct gr_tower_lazy_ctx_struct_tag;

/* A definition found to be a rational function of later definitions
   (lazy_merge.c): wherever elements involving the generator def_id and
   one of the generators trigger[0], trigger[1] (0: none) meet, the
   generator def_id becomes the value in their tower (a linear modulus),
   the generators of the value being moved before it. */
typedef struct
{
    ulong def_id;
    ulong trigger[2];
    gr_tower_lazy_elem_struct value;
    int have_value;     /* (0: computed by fn(value, data, ctx) on first use) */
    int (*fn)(gr_ptr, void *, gr_ctx_t);
    void * data;
    void (*data_clear)(void *, struct gr_tower_lazy_ctx_struct_tag *);
    /* (adds delta to the reference counts of the towers of the elements
       the data holds, by _gr_tower_lazy_ref_adjust: NULL if none) */
    void (*data_refs)(void *, slong delta);
}
gr_tower_lazy_rebase_struct;

/* The canonical root of unity of order q = l^e for a prime l (p = 0),
   or the principal q-th root of the positive integer p: the generator
   of a tower of its own (avoiding the construction of the algebraic
   number for the cache of algebraic numbers at every request). */
typedef struct
{
    fmpz p;
    ulong q;
    gr_tower_flat_struct * F;
    slong gid;
}
gr_tower_lazy_root_entry_struct;

/* A transcendental generator with a constant argument (exp(c), log(c)),
   or pi, hash-consed. */
typedef struct
{
    int kind;
    slong param;                    /* parameter of a special function (def_param) */
    fmpq x;
    fmpq y;                         /* imaginary part of the argument (log of a Gaussian rational) */
    slong nxs;                      /* functions of several arguments: the other (rational) arguments */
    fmpq * xs;
    gr_tower_flat_struct * F;
    slong gid;                      /* the generator in F */
    ulong def_id;                   /* its definition id (the same generator in other towers) */
}
gr_tower_lazy_const_trans_entry_struct;

typedef struct gr_tower_lazy_ctx_struct_tag
{
    gr_ctx_struct * base;
    gr_tower_flat_struct * trivial;
    gr_tower_lazy_const_trans_entry_struct * const_trans;
    slong num_const_trans;
    slong alloc_const_trans;
    gr_tower_flat_struct ** towers;
    slong num_towers;
    slong alloc_towers;
    gr_tower_flat_struct ** gc_queue;   /* towers whose element count dropped to zero */
    slong gc_len;
    slong gc_alloc;
    slong gc_collected;                 /* statistics */
    gr_tower_lazy_alias_struct * aliases;
    slong num_aliases;
    slong alloc_aliases;
    gr_tower_lazy_rebase_struct * rebases;
    slong num_rebases;
    slong alloc_rebases;
    ulong rebase_serial; /* (changed whenever the records change: the key of the caches of _rebase_pending) */
    int rebasing;        /* (a rebase in progress: the merges inside it do not start another) */
    gr_tower_lazy_qqbar_entry_struct * qqbars;
    slong num_qqbars;
    slong alloc_qqbars;
    gr_tower_lazy_root_entry_struct * roots;
    slong num_roots;
    slong alloc_roots;
    gr_tower_flat_struct * cyclo_F;   /* GR_TOWER_GENS_COMPOSITE_ROOTS: the tower of the canonical root of unity */
    ulong next_def_id;
    slong num_alg_names;   /* global generator names a1, a2, ... */
    slong num_trans_names; /* t1, t2, ... */
    int merge_flags;
    slong refcount;      /* the contexts sharing this state (views) */
    slong options[GR_TOWER_OPT_NUM_OPTIONS];   /* GR_TOWER_OPT_*, shared by the towers of the context */
    slong depth;         /* nesting of operations in progress (the subfield restrictions apply to the outermost) */
    int conj_depth;      /* nesting of conjugations in progress */
    int hyp_anchored;    /* nesting of the linking steps (Landen, quadratic transformations) in progress */
    struct gr_tower_lazy_theta_point_struct * theta_pts;   /* the points of the theta functions of z (lazy_modular.c) */
    slong num_theta_pts;
    slong alloc_theta_pts;
    struct gr_tower_lazy_mod_point_struct * mod_pts;   /* the values at the reduced points tau0 (lazy_modular.c) */
    slong num_mod_pts;
    slong alloc_mod_pts;
    ulong cache_clock;
#if FLINT_USES_PTHREAD
    pthread_mutex_t mutex;   /* recursive: every operation on the context holds it */
#endif
}
gr_tower_lazy_ctx_struct;

/* a point of the theta functions of z and the method of its evaluation
   (lazy_modular.c) */
typedef struct gr_tower_lazy_theta_point_struct
{
    gr_tower_lazy_elem_struct tau0;
    gr_tower_lazy_elem_struct z0;
    int kind;
    slong i, j, l;
    slong n, dv;
    int s, s2;
    int have_vals;
    gr_tower_lazy_elem_struct vals[4];
    ulong stamp;
    int pinned;     /* the values are in use: not evicted */
}
gr_tower_lazy_theta_point_struct;

/* the values at a reduced point tau0 (lazy_modular.c): for the ratios
   (0, 1) theta_2, theta_3, theta_4, lambda, E_2 */
typedef struct gr_tower_lazy_mod_point_struct
{
    gr_tower_lazy_elem_struct tau0;
    int have[2][5];
    gr_tower_lazy_elem_struct v[2][5];
    int anchor;                         /* 0: not chosen, 1: the anchor tau1, gam; 2: a new generator to be linked to it, 3: linked; -1: none */
    ulong root;                         /* the definition of the generator lambda the values come from (0: none, or not known) */
    ulong stamp;                        /* the last use (the values of the least recently used points are dropped) */
    gr_tower_lazy_elem_struct tau1;
    fmpz gam[4];
}
gr_tower_lazy_mod_point_struct;

/*
    The context data: the shared state (towers, registry, names, lock)
    and the field restriction of this context. Several contexts may share
    one state (gr_ctx_init_tower_lazy_view): the real and algebraic
    fields are views of a complex field, whose elements they share.
*/
typedef struct
{
    gr_tower_lazy_ctx_struct * L;
    int field_flags;     /* GR_TOWER_LAZY_REAL, GR_TOWER_LAZY_ALGEBRAIC */
}
gr_tower_lazy_view_struct;

#define VIEW(ctx) ((gr_tower_lazy_view_struct *) ((ctx)->data))
#define LAZY(ctx) (VIEW(ctx)->L)

/* towers continuing with at least this many generators beyond the
   prefix an element lives in are not extended for that element: new
   transcendental generators defined from it go into a copy of the
   prefix (see _gr_tower_lazy_trans_gen and _gr_tower_lazy_prefix_copy; a tower of only
   the generators the element involves was tried and is slower, since
   the relations with the skipped generators are then rediscovered) */
#define LAZY_FORK_SUFFIX 4

/* the same for the roots of an element (_gr_tower_lazy_root_ui), with a
   longer suffix: roots are found among the later generators of the
   tower, which a copy of the prefix would have to rediscover */
#define LAZY_ROOT_FORK_SUFFIX 12

#define TOWER(x) ((x)->F->T)

/* values of gr_tower_lazy_elem_struct.shallow */
#define LAZY_SHALLOW_NONE 0         /* the element owns its data */
#define LAZY_SHALLOW_COPY 1         /* a bitwise copy (gr set_shallow): the data may be another element's */
#define LAZY_SHALLOW_STALE 2        /* the data is owned by the tower (F->stale): the converted data of a
                                       shallow copy, never freed through the element */


/*
    Garbage collection of towers: every (non-shallow) element counts as a
    reference to its tower; a tower whose count drops to zero is queued,
    and the queue is processed at the end of the outermost operation on
    the context, when no internal pointers to towers are live. Towers
    referenced by the registries of roots of unity, radicals of integers
    and constant transcendentals (canonical definitions) are pinned; the
    entries of the algebraic number cache and the aliases of a collected
    tower are moved to another tower containing the definition, or
    dropped.
*/
#define LAZY_GC_BATCH 32






/* ---- types shared between the files ---- */

/* an operand of a dense operation as a polynomial in the generator:
   p points to the dense element's poly, or to tmp (a shallow view of a
   rational number) */
typedef struct
{
    const fmpq_poly_struct * p;
    fmpq_poly_struct tmp;
    fmpz c0;
}
gr_tower_lazy_dense_view_struct;

/* ---- macros shared between the files ---- */

/*
    Rational elements: those of the trivial tower (no generators), in the
    rational representation. Arithmetic between them needs neither the
    towers nor the lock (the trivial tower is pinned and never changes,
    and a rational element is never changed in place while it may be an
    operand), see the fast paths in lazy.c.
*/
#define LAZY_IS_RATIONAL(x, L) (LAZY_REPR(x) == LAZY_REPR_RATIONAL && !(x)->shallow)
#define LOCKED(call) int status; _gr_tower_lazy_lock(ctx); status = call; _gr_tower_lazy_unlock(ctx); return status;
#define LOCKED_T(T, call) T status; _gr_tower_lazy_lock(ctx); status = call; _gr_tower_lazy_unlock(ctx); return status;
#define LOCKED_V(call) _gr_tower_lazy_lock(ctx); call; _gr_tower_lazy_unlock(ctx);
#define DENSE_ADD 0
#define DENSE_SUB 1
#define DENSE_MUL 2
#define DENSE_DIV 3
/* the operation under the lock; then the dense forms of the result and
   the operands where they apply */
#define LOCKED_DENSE(call, a, b, c) \
    { \
        int status; \
        _gr_tower_lazy_lock(ctx); \
        status = call; \
        if (status == GR_SUCCESS && _gr_tower_lazy_tls_depth == 1) \
            _gr_tower_lazy_dense_after((gr_tower_lazy_elem_struct *) (a), (const gr_tower_lazy_elem_struct *) (b), (const gr_tower_lazy_elem_struct *) (c), ctx); \
        _gr_tower_lazy_unlock(ctx); \
        return status; \
    }

/* ---- variables shared between the files ---- */


/* ---- functions shared between the files ---- */

/* lazy.c */
gr_tower_flat_struct *
_gr_tower_lazy_new_tower(gr_ctx_t ctx);
void
_gr_tower_lazy_ref(gr_tower_flat_struct * F);
void
_gr_tower_lazy_unref(gr_tower_lazy_ctx_struct * L, gr_tower_flat_struct * F);
void
_gr_tower_lazy_ref_adjust(gr_tower_flat_struct * F, slong delta);
void
_gr_tower_lazy_new_def(gr_tower_gen_struct * g, gr_tower_t T, gr_ctx_t ctx);
ulong
_gr_tower_lazy_find_def(char ** name, int kind, slong param, const gr_tower_lazy_elem_struct * args, slong n,
    const gr_tower_struct * skip, gr_ctx_t ctx);
void
_gr_tower_lazy_set_def(gr_tower_gen_struct * g, gr_tower_t T, ulong def_id, const char * name, gr_ctx_t ctx);
void
_gr_tower_lazy_init(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void
_gr_tower_lazy_clear_data(gr_tower_lazy_elem_t x);
void
_gr_tower_lazy_clear(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void
_gr_tower_lazy_swap(gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx);
void
_gr_tower_lazy_fresh(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong level, gr_ctx_t ctx);
void
_gr_tower_lazy_update(gr_tower_lazy_elem_t x);
gr_tower_lazy_elem_struct *
_gr_tower_lazy_flat_view(const gr_tower_lazy_elem_t x_in);
void
_gr_tower_lazy_shrink(gr_tower_lazy_elem_t x);
slong
_gr_tower_lazy_alg_level(const gr_tower_lazy_elem_t x);
int
_gr_tower_lazy_set(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int _gr_tower_lazy_zero(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
void
_gr_tower_lazy_normalize_rational(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int _gr_tower_lazy_one(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int _gr_tower_lazy_set_si(gr_tower_lazy_elem_t res, slong c, gr_ctx_t ctx);
int _gr_tower_lazy_set_fmpz(gr_tower_lazy_elem_t res, const fmpz_t c, gr_ctx_t ctx);
int _gr_tower_lazy_set_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, gr_ctx_t ctx);
void
_gr_tower_lazy_set_gen(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong k, gr_ctx_t ctx);
void
_gr_tower_lazy_set_gen_d(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, slong d, gr_ctx_t ctx);
void
_gr_tower_lazy_fmpq_frac_part(fmpq_t r);
int
_gr_tower_lazy_qqbar_tower(gr_tower_flat_struct ** F_out, slong * gid, const qqbar_t x, int kind, ulong n, const fmpz_t p, gr_ctx_t ctx);
int
_gr_tower_lazy_set_qqbar(gr_tower_lazy_elem_t res, const qqbar_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_set_qqbar_structured(gr_tower_lazy_elem_t res, const qqbar_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_prime_power_root(gr_tower_lazy_elem_t res, const fmpz_t p, ulong l, ulong e, gr_ctx_t ctx);
int
_gr_tower_lazy_const_root(gr_tower_lazy_elem_t res, const fmpz_t p, ulong n, gr_ctx_t ctx);
int
_gr_tower_lazy_root_of_unity(gr_tower_lazy_elem_t res, slong p, ulong q, gr_ctx_t ctx);
int
_gr_tower_lazy_get_nested(gr_ptr res, gr_tower_lazy_elem_t x);
int
_gr_tower_lazy_get_qqbar(qqbar_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx);
slong *
_gr_tower_lazy_marked_by_creation(slong * count, const int * mark, gr_tower_t T);
int
_gr_tower_lazy_gen_is_constant_symbol(const gr_tower_gen_struct * g);
void
_gr_tower_lazy_mark_used_gens(int * mark, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_struct * mctx, gr_tower_t T);
void
_gr_tower_lazy_involved_gens(int * mark, const gr_tower_lazy_elem_t x, gr_tower_t T);
truth_t
_gr_tower_lazy_is_algebraic_repr(gr_tower_lazy_elem_t x);
int
_gr_tower_lazy_check_member(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int
_gr_tower_lazy_get_qqbar_impl(qqbar_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_get_acb_impl(acb_t res, const gr_tower_lazy_elem_t x, slong prec, gr_ctx_t ctx);
gr_tower_struct *
_gr_tower_lazy_get_tower_impl(slong * level, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
const fmpz_mpoly_q_struct *
_gr_tower_lazy_get_data_impl(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void
_gr_tower_lazy_ctx_stats_impl(gr_ctx_t ctx);

/* lazy_io.c */
int
_gr_tower_lazy_write(gr_stream_t out, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_gens(gr_vec_t vec, gr_ctx_t ctx);

/* lazy_merge.c */
int
_gr_tower_lazy_map_element(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_flat_struct * S, slong p, gr_tower_flat_struct * U, gr_ctx_t ctx);
void
_gr_tower_lazy_assign_def_ids(gr_tower_flat_struct * F, gr_ctx_t ctx);
gr_tower_flat_struct *
_gr_tower_lazy_common_tower(gr_tower_flat_struct * Fx, slong lx, gr_tower_flat_struct * Fy, slong ly, gr_ctx_t ctx);
int
_gr_tower_lazy_prim_coords(fmpq * v, const fmpz_mpoly_q_t x, const slong * var, const slong * C, slong m, slong D, const fmpz_mpoly_ctx_t mctx);
int
_gr_tower_lazy_common(gr_tower_flat_struct ** U, fmpz_mpoly_q_t xx, fmpz_mpoly_q_t yy,
    gr_tower_lazy_elem_t x, gr_tower_lazy_elem_t y, gr_ctx_t ctx);

/* lazy_arith.c */
void
_gr_tower_lazy_install(gr_tower_lazy_elem_t res, gr_tower_flat_struct * U, fmpz_mpoly_q_t r, gr_ctx_t ctx);
void
_gr_tower_lazy_prefix_copy(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* rebasing (lazy_merge.c): records that the generator of definition
   def_id equals value (a rational function of generators defined later),
   to be applied in the towers containing it and the generator trigger
   or trigger2 (0: none), when elements of a tower involving the
   definition and a trigger meet (_gr_tower_lazy_rebase_check) */
void _gr_tower_lazy_rebase_add(ulong def_id, ulong trigger, ulong trigger2, const gr_tower_lazy_elem_t value, gr_ctx_t ctx);
/* (the same with the value computed by fn(res, data, ctx) when first
   needed: a tower contains the definition and a trigger; data is
   cleared by data_clear(data, L), and the record ignored if fn fails) */
void _gr_tower_lazy_rebase_add_lazy(ulong def_id, ulong trigger, int (*fn)(gr_ptr, void *, gr_ctx_t), void * data, void (*data_clear)(void *, struct gr_tower_lazy_ctx_struct_tag *), void (*data_refs)(void *, slong), gr_ctx_t ctx);
void _gr_tower_lazy_rebase_clear(gr_tower_lazy_rebase_struct * r, gr_ctx_t ctx);
void _gr_tower_lazy_rebase_tower(gr_tower_flat_struct * F, gr_ctx_t ctx);
int _gr_tower_lazy_rebase_check(gr_tower_flat_struct * F, const fmpz_mpoly_q_struct * x, slong n, gr_ctx_t ctx);
ulong _gr_tower_lazy_gen_def_of(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_mul(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int
_gr_tower_lazy_neg(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t
_gr_tower_lazy_is_zero(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
truth_t
_gr_tower_lazy_equal(const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
truth_t
_gr_tower_lazy_is_one(const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_inv(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_div(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);

/* lazy_radical.c */
int
_gr_tower_lazy_root_fmpq(gr_tower_lazy_elem_t res, const fmpq_t c, ulong n, gr_ctx_t ctx);
int
_gr_tower_lazy_root_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, ulong n, gr_ctx_t ctx);

/* lazy_explog.c */
int
_gr_tower_lazy_pi(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int
_gr_tower_lazy_log_root_of_unity(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int * status, gr_ctx_t ctx);
int
_gr_tower_lazy_trans_gen(gr_tower_lazy_elem_t res, gr_tower_lazy_elem_t x, int kind, gr_ctx_t ctx);
int
_gr_tower_lazy_exp_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, int kind, gr_ctx_t ctx);
int
_gr_tower_lazy_exp(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_log(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* lazy_conj.c */
int
_gr_tower_lazy_get_fmpq(fmpq_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_get_fmpz(fmpz_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_get_si(slong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_get_ui(ulong * res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_i(gr_tower_lazy_elem_t res, gr_ctx_t ctx);
int
_gr_tower_lazy_gen_is_i(gr_tower_struct * T, slong d);
void
_gr_tower_lazy_set_flat(gr_tower_lazy_elem_t res, gr_tower_flat_struct * F, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_ctx_t ctx);
int
_gr_tower_lazy_conj(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_transfer(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x_in, gr_ctx_t x_ctx, gr_ctx_t ctx);
int
_gr_tower_lazy_abs(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);
int
_gr_tower_lazy_pow(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y_in, gr_ctx_t ctx);
int
_gr_tower_lazy_sqrt(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* lazy_roots.c */
int
_gr_tower_lazy_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t poly, int flags, gr_ctx_t ctx);
/* (flag of _gr_tower_lazy_poly_roots: the real roots suffice, as for the
   real fields; other roots may still be returned, the caller filters) */
#define LAZY_ROOTS_REAL_ONLY (1 << 30)
int
_gr_tower_lazy_poly_root_near(gr_tower_lazy_elem_t res, const gr_poly_t f, const acb_t ref, int pm, gr_ctx_t ctx);

/* lazy_dense.c */
int
_gr_tower_lazy_poly_mullow(gr_tower_lazy_elem_struct * res, const gr_tower_lazy_elem_struct * poly1, slong len1,
    const gr_tower_lazy_elem_struct * poly2, slong len2, slong n, gr_ctx_t ctx);
int
_gr_tower_lazy_mat_mul(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
int
_gr_tower_lazy_mat_det(gr_tower_lazy_elem_t res, const gr_mat_t A, gr_ctx_t ctx);
int
_gr_tower_lazy_mat_nonsingular_solve(gr_mat_t X, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx);
void
_gr_tower_lazy_try_dense(gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void
_gr_tower_lazy_dense_get_flat(fmpz_mpoly_q_t res, const gr_tower_lazy_elem_t x, gr_tower_flat_struct * F);
void
_gr_tower_lazy_dense_view_dense(gr_tower_lazy_dense_view_struct * V, const gr_tower_lazy_elem_t x);
void _gr_tower_lazy_dense_view_fmpq(gr_tower_lazy_dense_view_struct * V, const fmpq_t q);
int _gr_tower_lazy_get_fmpq_poly_impl(fmpq_poly_t res, fmpz_poly_t modulus, gr_tower_lazy_elem_t x, gr_ctx_t ctx);
void _gr_tower_lazy_dense_after(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int
_gr_tower_lazy_dense_view(gr_tower_lazy_dense_view_struct * V, const gr_tower_lazy_elem_t x, const _gr_tower_dense_field_struct * D, const gr_tower_flat_struct * F, gr_ctx_t ctx);
int
_gr_tower_lazy_dense_op(gr_tower_lazy_elem_t res, const gr_tower_lazy_dense_view_struct * X, const gr_tower_lazy_dense_view_struct * Y, int op, const _gr_tower_dense_field_struct * D);
int
_gr_tower_lazy_dense_binary(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, int op, gr_ctx_t ctx);
int
_gr_tower_lazy_dense_scalar(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, int op, gr_ctx_t ctx);
int
_gr_tower_lazy_dense_candidate(const gr_tower_lazy_elem_t x);

/* lazy_real.c */
int
_gr_tower_lazy_tan_pi_gen(gr_tower_lazy_elem_t res, ulong M, qqbar_t z, gr_ctx_t ctx);
int
_gr_tower_lazy_realify(gr_tower_lazy_elem_t x, gr_ctx_t ctx);

/* lazy_arith.c (defined by macros) */
int _gr_tower_lazy_add(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int _gr_tower_lazy_sub(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const gr_tower_lazy_elem_t y, gr_ctx_t ctx);
int _gr_tower_lazy_add_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int _gr_tower_lazy_add_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int _gr_tower_lazy_add_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int _gr_tower_lazy_add_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int _gr_tower_lazy_sub_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int _gr_tower_lazy_sub_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int _gr_tower_lazy_sub_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int _gr_tower_lazy_sub_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int _gr_tower_lazy_mul_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int _gr_tower_lazy_mul_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int _gr_tower_lazy_mul_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int _gr_tower_lazy_mul_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);
int _gr_tower_lazy_div_fmpq(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpq_t c, gr_ctx_t ctx);
int _gr_tower_lazy_div_fmpz(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, const fmpz_t c, gr_ctx_t ctx);
int _gr_tower_lazy_div_si(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, slong c, gr_ctx_t ctx);
int _gr_tower_lazy_div_ui(gr_tower_lazy_elem_t res, const gr_tower_lazy_elem_t x, ulong c, gr_ctx_t ctx);

#endif
