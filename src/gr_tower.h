/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef GR_TOWER_H
#define GR_TOWER_H

#ifdef GR_TOWER_INLINES_C
#define GR_TOWER_INLINE
#else
#define GR_TOWER_INLINE static inline
#endif

#include "acb_types.h"
#include "fmpz_mpoly_q.h"
#include "qqbar.h"
#include "gr.h"
#include "gr_poly.h"
#include "gr_mat.h"

#ifdef __cplusplus
extern "C" {
#endif

/*
    A tower F_0 subset F_1 subset ... subset F_n of fields, where F_0 is a
    base field (typically QQ) and each F_k = F_{k-1}[x]/(m_k) is generated
    by an algebraic element a_k with monic minimal polynomial m_k over F_{k-1}.
    Each generator is embedded in C by a numerical enclosure, which fixes
    the choice of conjugate.

    The polynomial m_k need not be known to be irreducible: the ring
    F_{k-1}[x]/(m_k) is treated as a field (dynamic evaluation), and
    whenever an inversion exposes a factorization of m_k, the tower is
    refined by replacing m_k with the factor vanishing at the enclosure.
    Refinement does not invalidate existing elements.
*/

#define GR_TOWER_ALGEBRAIC 0

/* kinds of transcendental generators */
#define GR_TOWER_EXP 1
#define GR_TOWER_LOG 2
#define GR_TOWER_PI 3
#define GR_TOWER_FREE 4
#define GR_TOWER_ROOT 5               /* def_kind only: principal def_param-th root of arg */
#define GR_TOWER_ROOT_OF_UNITY 6      /* def_kind only: exp(2 pi i / def_param) */
#define GR_TOWER_TAN 7                /* tan(arg), arg real */
#define GR_TOWER_ATAN 8               /* atan(arg), arg real */

/* special functions (see gr_tower_adjoin_special) */
#define GR_TOWER_GAMMA 9              /* Gamma(arg) */
#define GR_TOWER_ERF 10               /* erf(arg) (def_param = 0), erfi(arg) (def_param = 1) */
#define GR_TOWER_LAMBERTW 11          /* W_k(arg), k = def_param */
#define GR_TOWER_POLYGAMMA 12         /* psi^(m)(arg), m = def_param >= 0 */
#define GR_TOWER_POLYLOG 13           /* Li_s(arg), s = def_param */
#define GR_TOWER_ZETA 14              /* zeta(arg) */
#define GR_TOWER_ELLIPTIC_K 15        /* K(arg), complete elliptic integral of the first kind (parameter m = arg) */
#define GR_TOWER_ELLIPTIC_E 16        /* E(arg), of the second kind */
#define GR_TOWER_CONSTANT 17          /* a named constant (def_param = GR_TOWER_CONST_*), no argument */
#define GR_TOWER_TAN_PI 18            /* def_kind only: tan(pi / def_param), a real algebraic number (the tangent
                                         normal form of the real trigonometric constants) */

#define GR_TOWER_CONST_EULER 1        /* Euler's constant gamma */
#define GR_TOWER_CONST_CATALAN 2      /* Catalan's constant G */

#define GR_TOWER_KIND_IS_SPECIAL(kind) ((kind) >= GR_TOWER_GAMMA && (kind) <= GR_TOWER_CONSTANT)

/* the transcendental kinds with an argument */
#define GR_TOWER_KIND_HAS_ARG(kind) ((kind) == GR_TOWER_EXP || (kind) == GR_TOWER_LOG || (kind) == GR_TOWER_TAN || (kind) == GR_TOWER_ATAN || \
                                     (GR_TOWER_KIND_IS_SPECIAL(kind) && (kind) != GR_TOWER_CONSTANT))

/* Status of a generator. For an algebraic generator: the modulus m_k is
   known to be irreducible over F_{k-1} (PROVEN) or only assumed to be
   (DYNAMIC). For a transcendental generator: PROVEN (alias INDEPENDENT)
   means that it is known to be algebraically independent of all the
   generators preceding it in the definition order (a proof of
   transcendence alone, as of e by Lindemann, does not give this status
   when a transcendental generator precedes it); SCHANUEL and CONJECTURAL
   mean that the independence is conjectural. */
#define GR_TOWER_STATUS_PROVEN 0      /* m_k is known to be irreducible over F_{k-1} / t is independent of the earlier generators */
#define GR_TOWER_STATUS_INDEPENDENT GR_TOWER_STATUS_PROVEN
#define GR_TOWER_STATUS_DYNAMIC 1     /* m_k is only assumed irreducible */
#define GR_TOWER_STATUS_SCHANUEL 2    /* independence conjectural: no relation found at the certification effort */
#define GR_TOWER_STATUS_CONJECTURAL 3 /* special function value: algebraic independence conjectural (no relation
                                         following from the known identities), beyond Schanuel's conjecture */

/*
    A generator of a tower. Generators are stored in definition order;
    an algebraic generator a_k has an index k in the nested chain and a
    transcendental generator t_j an index j in the base field. A
    transcendental generator may become algebraic when a relation is
    discovered (kind changes, def_kind and the argument are kept).
*/
/*
    A flat element of a tower together with the polynomial context it is
    stored in (the flat machinery replaces its context when the tower
    grows beyond its capacity or is reordered; old contexts stay alive,
    and elements in them are converted on use, gr_tower_flat_convert).
*/
typedef struct
{
    fmpz_mpoly_q_struct data;
    fmpz_mpoly_ctx_struct * mctx;
}
gr_tower_flat_elem_struct;

typedef gr_tower_flat_elem_struct gr_tower_flat_elem_t[1];

typedef struct
{
    int kind;                     /* current kind: GR_TOWER_ALGEBRAIC or a transcendental kind */
    int def_kind;                 /* kind of the definition */
    slong def_param;              /* parameter of the definition (the n of an n-th root) */
    int status;
    slong index;                  /* k (algebraic) or j (transcendental), 1-based */
    gr_ctx_struct * ctx;          /* algebraic: F_k = F_{k-1}[x]/(m_k), owned */
    gr_tower_flat_elem_struct arg;    /* EXP/LOG/...: the argument, a flat element of the tower (arg.mctx NULL if none) */
    fmpz_poly_struct * origin;    /* an integer polynomial the generator is a root of (NULL if unknown); conjugate
                                     roots of one polynomial are recognized by it when towers are merged */
    acb_struct enclosure;
    slong enclosure_prec;         /* precision at which the enclosure was last refined */
    char * name;
    ulong def_id;                 /* identifies the generator across towers (0 = unassigned) */
    slong def_order;              /* position in the definition order */
    slong gid;                    /* identifies the generator within the tower (stable under reorderings) */
    ulong proof_version;          /* 1 + the moduli version at the last modular proof attempt (0 = none) */
    int real;                     /* the generator is asserted real (tan(u), atan(u) for the real u the
                                     caller promised): its enclosures get an exactly zero imaginary part */
}
gr_tower_gen_struct;

/*
    Flat machinery of a tower: multivariate rational functions (fmpz_mpoly_q)
    in all generators, reduced modulo the triangular set with denominators
    cleared. The generator with definition order d corresponds to variable
    cap - 1 - d, so that later generators are more significant in lex order
    (which makes the triangular set a Groebner basis, each modulus involving
    only earlier generators) and the tower can grow up to the capacity
    without renumbering. When the capacity is exceeded, a new polynomial
    context is created; elements in an old context are converted with
    gr_tower_flat_convert (old contexts are kept alive).
*/
typedef struct
{
    struct gr_tower_struct_tag * T;
    fmpz_mpoly_ctx_struct * mctx;
    slong cap;
    slong * var_gid;              /* gid of the generator represented by each variable (-1 if none) */
    fmpz_mpoly_struct ** ideal;
    slong ideal_len;
    ulong ideal_version;
    ulong ideal_moduli_version;   /* moduli version of the tower for which the existing ideal elements are valid */
    int ideal_lc_const;           /* all leading coefficients of the ideal are integers */
    fmpz_poly_struct ** ideal_univar;  /* ideal[k] as a monic univariate integer polynomial in its main variable (NULL if it is not one) */
    fmpz_mpoly_ctx_struct ** old_mctx;
    slong ** old_var_gid;
    slong num_old_mctx;
    ulong layout_version;         /* structure version of the tower for which var_gid was made */
    fmpz_mpoly_q_struct ** stale; /* elements owned by F (converted data of shallow copies), freed with F */
    fmpz_mpoly_ctx_struct ** stale_mctx;
    slong num_stale;
    fmpz_mpoly_q_struct * inv_cache;   /* the last dense inverse: an element and its inverse (NULL: none), for repeated divisions by the same element */
    fmpz_mpoly_ctx_struct * inv_cache_mctx;
    ulong inv_cache_version[2];   /* ideal_version and layout_version of inv_cache */
    void * modp;                  /* primes modulo which the moduli split, with their roots (dense inverses; NULL: none) */
    slong refs;                   /* lazy fields: elements living in the tower */
    int gc;                       /* lazy fields: GR_TOWER_GC_* flags */
    void * dense_fields;          /* lazy fields: descriptors of the dense forms of elements (_gr_tower_dense_field_struct), freed with F */
    ulong primitive_tried[2];     /* lazy fields: version + 1 and length of the tower when a primitive element was last sought (0: never) */
}
gr_tower_flat_struct;

/* flags of towers of lazy fields (gr_tower_flat_struct.gc) */
#define GR_TOWER_GC_PINNED 1      /* referenced by a registry: never collected */
#define GR_TOWER_GC_CANDIDATE 2   /* queued for the next collection */
#define GR_TOWER_GC_DEAD 4        /* being collected */

typedef gr_tower_flat_struct gr_tower_flat_t[1];

typedef struct gr_tower_struct_tag
{
    gr_ctx_struct * consts;       /* the constant field (QQ), not owned */
    gr_ctx_struct * base;         /* F_0: consts if there are no transcendental generators,
                                     otherwise an owned rational function field context */
    slong num_gens;
    slong alloc_gens;
    gr_tower_gen_struct * gens;   /* in definition order */
    slong length;                 /* number of algebraic generators */
    slong * alg;                  /* alg[k-1] = definition order of a_k */
    slong num_trans;
    slong * trans;                /* trans[j-1] = definition order of t_j */
    slong next_gid;
    ulong version;                /* incremented on every refinement or rebuild */
    ulong moduli_version;         /* incremented when the modulus of an existing step changes (not when steps are appended) */
    ulong structure_version;      /* incremented when the definition order changes */
    gr_tower_flat_struct flat;
    const slong * options;        /* GR_TOWER_OPT_* (not owned: the defaults, or a lazy context's table) */
    int caps;                     /* GR_TOWER_CAP_*: what the constant field supports, set at initialization */
    int keep_retired;             /* keep the nested contexts superseded by rebuilds until the tower is
                                     cleared (towers made by gr_tower_init; the towers of lazy fields and
                                     internal temporaries never hold nested elements across a rebuild) */
    gr_ctx_struct ** retired;     /* those contexts, in the order in which they can be cleared */
    slong num_retired;
}
gr_tower_struct;

/*
    Capabilities of the constant field K of a tower, from which the
    algorithms select what applies (currently K = QQ, which has all of
    them; gr_tower_init rejects other fields). A future constant field
    would set the bits it supports: a finite field has places and
    factoring but no complex embedding; the p-adic numbers have an
    embedding (Hensel refinement of roots) but no involution or order.
*/
#define GR_TOWER_CAP_RATIONAL 1       /* K is the rational field (integer arithmetic, Gauss sums, cyclotomic rules) */
#define GR_TOWER_CAP_FACTOR 2         /* univariate polynomials over K can be factored */
#define GR_TOWER_CAP_PLACES 4         /* K has finite places (modular proofs of irreducibility) */
#define GR_TOWER_CAP_EMBEDDING 8      /* the elements have numerical enclosures (acb) refined by Newton's method */
#define GR_TOWER_CAP_INVOLUTION 16    /* the embedding has a conjugation (real parts, realness) */
#define GR_TOWER_CAP_ORDERED 32       /* the fixed field of the involution is ordered (signs, comparisons, principal branches) */
#define GR_TOWER_HAS_CAP(T, cap) (((T)->caps & (cap)) != 0)
/* the base field K(t_1, ..., t_r) is K itself: no transcendental generators */
#define GR_TOWER_BASE_IS_CONSTS(T) ((T)->num_trans == 0)

typedef gr_tower_struct gr_tower_t[1];

/*
    Tuning options. A tower reads them from the table its options pointer
    refers to: gr_tower_default_options, or the table of the lazy context
    that created it (gr_tower_lazy_ctx_set_option), shared by its towers.
*/
enum
{
    GR_TOWER_OPT_VERBOSE,                   /* 1: print the relations found between generators; 2: also statistics */
    GR_TOWER_OPT_PRINT_FLAGS,               /* GR_TOWER_PRINT_* (lazy fields) */
    GR_TOWER_OPT_PRINT_DIGITS,              /* digits of numerical values printed */
    GR_TOWER_OPT_PREC_LIMIT,                /* bits: numerical separation from zero before the exact tests */
    GR_TOWER_OPT_CERTIFY_PREC_LIMIT,        /* bits: certification of roots, signs and real parts, and the relation searches */
    GR_TOWER_OPT_NUMERIC_PREC_LIMIT,        /* bits: numerical evaluation pushed this far when the exact methods are undecided */
    GR_TOWER_OPT_SMOOTH_LIMIT,              /* number of primes for the trial division of radicands and of Gaussian integers in logarithms */
    GR_TOWER_OPT_EXPRESS_DEGREE_LIMIT,      /* lattice-based expression searches in fields up to this degree */
    GR_TOWER_OPT_EXPRESS_PREC,              /* bits (plus 4 per degree of the field) for the lattice-based expression searches */
    GR_TOWER_OPT_TRAGER_DEGREE_LIMIT,       /* Trager's method for steps up to this degree */
    GR_TOWER_OPT_FACTOR_DEGREE_LIMIT,       /* factoring over towers up to this degree */
    GR_TOWER_OPT_ROOTS_FACTOR_DEGREE_LIMIT, /* lazy fields: roots of polynomials by factoring up to this degree */
    GR_TOWER_OPT_MODULAR_TRIES,             /* primes tried by modular irreducibility proofs */
    GR_TOWER_OPT_MODULAR_TERMS_LIMIT,       /* terms of a modulus beyond which modular proofs are not attempted */
    GR_TOWER_OPT_NO_ROOTS_TRIES,            /* places tried for modular proofs of no roots */
    GR_TOWER_OPT_RATIONALIZE_LIMIT,         /* norm degree for rationalizing denominators */
    GR_TOWER_OPT_RELATION_COST_LIMIT,       /* multiplicative relations costlier to verify are not verified */
    GR_TOWER_OPT_SQRT_BUDGET,               /* recursive square root searches per call */
    GR_TOWER_OPT_GAUSS_SUM_LIMIT,           /* square roots as Gauss sums up to this order */
    GR_TOWER_OPT_DENSE_LIMIT,               /* Kronecker array length for dense products */
    GR_TOWER_OPT_INV_DENSE_DEGREE_LIMIT,    /* dense inverses in fields up to this degree */
    GR_TOWER_OPT_INV_LINEAR_DEGREE_LIMIT,   /* lazy fields: inverses by a linear system up to this degree */
    GR_TOWER_OPT_DEFERRED_IDEAL_SIZE,       /* reductions by moduli of more than this many terms are deferred */
    GR_TOWER_OPT_GROW_THRESHOLD,            /* towers grown in place beyond this many generators */
    GR_TOWER_OPT_ROOT_OF_UNITY_ORDER_LIMIT, /* exp(r pi i) recognized as a root of unity up to this order */
    GR_TOWER_OPT_CYCLOTOMIC_ORDER_LIMIT,    /* roots of unity as generators up to this order */
    GR_TOWER_OPT_REALIFY_ORDER_LIMIT,       /* real forms of elements with roots of unity up to this order */
    GR_TOWER_OPT_TRIG_PI_LIMIT,             /* tangent normal form for denominators up to this size */
    GR_TOWER_OPT_TRIG_ALGEBRAIC_LIMIT,      /* trigonometric values at r pi as algebraic numbers for denominators up to this size */
    GR_TOWER_OPT_TRIG_FORM,                 /* GR_TOWER_TRIG_*: trigonometric constants in the complex field */
    GR_TOWER_OPT_SPLIT_IMAGINARY,           /* sqrt(-A) as i sqrt(A) rather than one generator */
    GR_TOWER_OPT_COMPOSITE_RADICALS,        /* sqrt(A), A squarefree, as one generator */
    GR_TOWER_OPT_CYCLOTOMIC_DEGREE_LIMIT,   /* one root of unity zeta_N for phi(N) up to this (0: prime powers) */
    GR_TOWER_OPT_SAME_ORIGIN_DEGREE_LIMIT,  /* merges: conjugate roots of one polynomial recognized up to this degree */
    GR_TOWER_OPT_ANNIHILATING_FLAT_LIMIT,   /* annihilating polynomials by flat linear algebra up to this degree */
    GR_TOWER_OPT_GAMMA_LATTICE_LIMIT,       /* normal form of Gamma at rationals up to this denominator */
    GR_TOWER_OPT_GAMMA_LINE_LIMIT,          /* Gamma relations on rational lines for multipliers up to this size */
    GR_TOWER_OPT_HURWITZ_LATTICE_LIMIT,     /* normal form of Hurwitz zeta values up to this level */
    GR_TOWER_OPT_HURWITZ_WEIGHT_LIMIT,      /* Hurwitz zeta relations up to this weight */
    GR_TOWER_OPT_SPECIAL_RELATION_LEVEL_LIMIT, /* Gamma and Hurwitz relations in the zero test up to this level */
    GR_TOWER_OPT_GAUSS_DIGAMMA_LIMIT,       /* Gauss's digamma theorem up to this denominator */
    GR_TOWER_OPT_DENSE_FORM_DEGREE_LIMIT,   /* lazy elements of number fields of one generator up to this degree in dense form (0: never) */
    GR_TOWER_OPT_DENSE_FORM_SPARSITY,       /* ... and, if longer than 256, only with at least 1/this of their coefficients nonzero (0: always) */
    GR_TOWER_OPT_PRIMITIVE_DEGREE_LIMIT,    /* lazy fields: number fields of several steps up to this degree merged into Q(theta) (0: never) */
    GR_TOWER_OPT_SPLIT_DEGREE_LIMIT,        /* lazy fields: zero tests in towers of algebraic numbers beyond this degree split independent parts (0: never) */
    GR_TOWER_OPT_MINPOLY_DEGREE_LIMIT,      /* lazy fields: equality of split parts by the minimal polynomial of one of degree up to this */
    GR_TOWER_OPT_INV_DENSE_ALG,             /* dense inverses in fields of several generators: 0 automatic, 1 modular only, 2 linear algebra only (testing) */
    GR_TOWER_OPT_NUM_OPTIONS
};

#define GR_TOWER_TRIG_EXPONENTIAL 0
#define GR_TOWER_TRIG_TANGENT 1

/* element printing in lazy fields (GR_TOWER_OPT_PRINT_FLAGS; gr_tower_lazy.h) */
#define GR_TOWER_PRINT_NUMERIC 1      /* a numerical value */
#define GR_TOWER_PRINT_SYMBOLIC 2     /* the expression in the generators */
#define GR_TOWER_PRINT_DEFS 4         /* the definitions of the generators involved */
#define GR_TOWER_PRINT_DIGITS_DEFAULT 6

FLINT_DLL extern const slong gr_tower_default_options[GR_TOWER_OPT_NUM_OPTIONS];

/* names, defaults and validation of the options; a table for
   gr_tower_set_options is made with gr_tower_options_init (a copy of
   the defaults) and changed with gr_tower_options_set (GR_DOMAIN for a
   value out of range) */
const char * gr_tower_option_name(slong option);
slong gr_tower_option_find(const char * name);
slong gr_tower_option_default(slong option);
int gr_tower_option_valid(slong option, slong value);
void gr_tower_options_init(slong * options);
int gr_tower_options_set(slong * options, slong option, slong value);

#define GR_TOWER_OPTION(T, k) ((T)->options[k])

#define GR_TOWER_GEN(T, d) ((T)->gens + (d))
#define GR_TOWER_STEP(T, k) ((T)->gens + (T)->alg[k])          /* k 0-based */
#define GR_TOWER_STEP_CTX(T, k) (GR_TOWER_STEP(T, k)->ctx)
#define GR_TOWER_STEP_ENCLOSURE(T, k) (&(GR_TOWER_STEP(T, k)->enclosure))
#define GR_TOWER_TRANS(T, j) ((T)->gens + (T)->trans[j])       /* j 0-based */

/* Generators in definition order. A prefix of length p of a tower is the
   set of its first p generators in this order; it is itself a tower. */
GR_TOWER_INLINE slong gr_tower_num_gens(const gr_tower_t T) { return T->num_gens; }

/* Code of the generator with definition order d: k >= 1 for the algebraic
   generator a_k, -j for the transcendental generator t_j. */
GR_TOWER_INLINE slong gr_tower_order_code(const gr_tower_t T, slong d)
{
    const gr_tower_gen_struct * g = T->gens + d;
    return (g->kind == GR_TOWER_ALGEBRAIC) ? g->index : -g->index;
}

/* Number of algebraic, respectively transcendental, generators among the
   first p generators. */
slong gr_tower_prefix_length(const gr_tower_t T, slong p);
slong gr_tower_prefix_num_trans(const gr_tower_t T, slong p);

/* Definition order of the generator with the given gid, or -1. */
slong gr_tower_gid_order(const gr_tower_t T, slong gid);

#define GR_TOWER_FLAT_VAR_D(F, d) ((F)->cap - 1 - (d))
#define GR_TOWER_FLAT_VAR(F, k) GR_TOWER_FLAT_VAR_D(F, (F)->T->alg[(k) - 1])
#define GR_TOWER_FLAT_TVAR(F, j) GR_TOWER_FLAT_VAR_D(F, (F)->T->trans[(j) - 1])
#define GR_TOWER_FLAT(T) (&(T)->flat)

/* Memory management and basic access */

void gr_tower_init(gr_tower_t T, gr_ctx_t base);
void gr_tower_set_options(gr_tower_t T, const slong * options);
void gr_tower_clear(gr_tower_t T);

/* Heap-allocated towers, a string representation and non-inline
   accessors (for language bindings). */
gr_tower_struct * gr_tower_heap_init(gr_ctx_t base);
void gr_tower_heap_clear(gr_tower_struct * T);
int gr_tower_get_str(char ** s, const gr_tower_t T);
gr_ctx_struct * gr_tower_field_ptr(const gr_tower_t T);
slong gr_tower_length_si(const gr_tower_t T);
slong gr_tower_num_gens_si(const gr_tower_t T);
slong gr_tower_num_trans_si(const gr_tower_t T);
ulong gr_tower_version(const gr_tower_t T);
const char * gr_tower_gen_name(const gr_tower_t T, slong d);
int gr_tower_gen_status(const gr_tower_t T, slong d);
int gr_tower_gen_kind(const gr_tower_t T, slong d);
int gr_tower_gen_get(gr_ptr res, const gr_tower_t T, slong d);
void gr_tower_set(gr_tower_t res, const gr_tower_t T);

/* Sets res to a copy of the first p generators (in definition order) of T. */
void gr_tower_set_prefix(gr_tower_t res, const gr_tower_t T, slong p);
void gr_tower_set_subset(gr_tower_t res, const gr_tower_t T, const int * mark);

/* Definition ids identify generators across towers (gr_tower_absorb
   identifies generators with equal nonzero ids); the lazy fields assign
   them, and a user of the tower API may assign its own. */
void gr_tower_gen_set_def_id(gr_tower_t T, slong d, ulong def_id);
/* Definition order of the generator of T with the given definition id, or -1. */
slong gr_tower_find_def_order(const gr_tower_t T, ulong def_id);

/* Index (1-based) of the algebraic generator of T with the given definition id, or 0. */
slong gr_tower_find_def(const gr_tower_t T, ulong def_id);

GR_TOWER_INLINE slong gr_tower_length(const gr_tower_t T) { return T->length; }
GR_TOWER_INLINE gr_ctx_struct * gr_tower_base(const gr_tower_t T) { return T->base; }

/* Context of F_k, 0 <= k <= n */
GR_TOWER_INLINE gr_ctx_struct * gr_tower_field_at(const gr_tower_t T, slong k)
{
    return (k == 0) ? T->base : GR_TOWER_STEP(T, k - 1)->ctx;
}

/* Context of the top field F_n */
GR_TOWER_INLINE gr_ctx_struct * gr_tower_field(const gr_tower_t T)
{
    return gr_tower_field_at(T, T->length);
}

/* Degree of m_k over F_{k-1}, 1 <= k <= n */
slong gr_tower_step_degree(const gr_tower_t T, slong k);
slong gr_tower_degree_at(const gr_tower_t T, slong k);   /* [F_k : F_0], saturating at WORD_MAX */

/* Degree [F_n : F_0] */
slong gr_tower_degree(const gr_tower_t T);

const gr_poly_struct * gr_tower_step_minpoly(const gr_tower_t T, slong k);

int gr_tower_write(gr_stream_t out, const gr_tower_t T);
void gr_tower_print(const gr_tower_t T);

/* Adjoining generators */

/*
    Adjoins the root of the monic polynomial m over the current top field
    nearest to the midpoint of z: the root is certified (by interval
    Newton iteration) to be the only root of m in a disk around the
    midpoint of z, whose radius is that of z inflated by up to |z|/4 when
    z is tighter than the roots can be separated. z may thus be a
    rigorous enclosure or an approximation to a few digits. Returns
    GR_UNABLE if no such disk is found (several roots near z, or
    insufficient precision), GR_DOMAIN if m is not monic of positive
    degree.
*/
int gr_tower_adjoin_algebraic(gr_tower_t T, const gr_poly_t m, const acb_t z, int status, const char * name);

/* Modular irreducibility proofs (towers over QQ): see modular.c */
int gr_tower_prove_step_modular(gr_tower_t T, slong k, slong tries);
int gr_tower_prove_modular(gr_tower_t T, slong tries);
int gr_tower_binomial_irreducible_modular(gr_tower_t T, const fmpz_mpoly_q_t a, const fmpz_mpoly_ctx_t actx, ulong n, slong tries);
int gr_tower_poly_no_roots_modular(const gr_poly_t g, gr_tower_t T, slong tries);
int gr_tower_prove_step_trager(gr_tower_t T, slong k, slong degree_limit);

/* Factorization of polynomials over F_k (Trager's method; the steps up
   to k must be provable, the base field QQ or QQ(t_1, ..., t_r)) */
int gr_tower_poly_factor(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T);
int gr_tower_poly_factor_limit(gr_ptr c, gr_vec_t fac, fmpz_vec_t mult, const gr_poly_t f, slong k, slong degree_limit, gr_tower_t T);
int gr_tower_poly_roots(gr_vec_t roots, fmpz_vec_t mult, const gr_poly_t f, slong k, gr_tower_t T);

/* Adjoins the algebraic number x (assuming the base field is QQ). */
int gr_tower_adjoin_qqbar(gr_tower_t T, const qqbar_t x, const char * name);

/* Adjoins the principal n-th root of the element x of the top field. */
int gr_tower_adjoin_root_ui(gr_tower_t T, gr_srcptr x, ulong n, const char * name);

/* Adjoins the root of unity exp(2 pi i / n), recorded as such
   (GR_TOWER_ROOT_OF_UNITY) so that roots of unity of orders dividing n
   can be identified with its powers. */
int gr_tower_adjoin_root_of_unity(gr_tower_t T, ulong n, const char * name);

/* Adjoins the principal n-th root of the positive integer p, recorded
   as such (GR_TOWER_ROOT with a constant argument). */
int gr_tower_adjoin_root_fmpz(gr_tower_t T, const fmpz_t p, ulong n, const char * name);

GR_TOWER_INLINE int gr_tower_adjoin_sqrt(gr_tower_t T, gr_srcptr x, const char * name)
{
    return gr_tower_adjoin_root_ui(T, x, 2, name);
}

/*
    Adjoins a transcendental generator: pi, exp(u) or log(u). The argument
    u is an element of the top field of the tower (nested representation)
    or, in the _flat variants, a flat element in the current or a
    superseded context of the tower's flat machinery. The base field of the
    tower changes, so the nested contexts are rebuilt: nested elements held
    by the caller (u among them) remain elements of the superseded
    contexts, which a tower made by gr_tower_init keeps until it is
    cleared, so that they can still be used with those contexts and
    cleared, but they are not elements of the new top field; flat
    elements remain valid (they may need conversion). The generator is initially of status
    GR_TOWER_STATUS_SCHANUEL (or PROVEN for pi). Adjoining exp(0), log(0)
    or log(1) (whose values are not transcendental) fails with GR_DOMAIN;
    if it cannot be decided whether u is 0 or 1, GR_UNABLE.
*/
int gr_tower_adjoin_pi(gr_tower_t T, const char * name);
/* Adjoins a free transcendental generator (a formal variable, GR_TOWER_FREE):
   algebraically independent of all generators by definition, with no
   numerical value; the top field is then a rational function field. */
int gr_tower_adjoin_free(gr_tower_t T, const char * name);
int gr_tower_adjoin_exp(gr_tower_t T, gr_srcptr u, const char * name);
int gr_tower_adjoin_log(gr_tower_t T, gr_srcptr u, const char * name);
int gr_tower_adjoin_exp_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);
int gr_tower_adjoin_log_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);

/* The real trigonometric generators tan(u) and atan(u), for a real u
   (which the caller must ensure; u = 0 is rejected with GR_DOMAIN, as is
   u known to be non-real numerically). The relations between them, pi
   and the exponentials and logarithms are found by Richardson's
   algorithm (on the angles 2u of tan(u), 2 atan(u) and pi). */
int gr_tower_adjoin_tan(gr_tower_t T, gr_srcptr u, const char * name);
int gr_tower_adjoin_atan(gr_tower_t T, gr_srcptr u, const char * name);
int gr_tower_adjoin_tan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);
int gr_tower_adjoin_atan_flat(gr_tower_t T, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);

/*
    Adjoins the value of a special function (kind GR_TOWER_GAMMA, ...,
    GR_TOWER_CONSTANT, with the parameter param where the kind has one)
    at u, an element of the top field (nested representation) or a flat
    element, as a transcendental generator of status
    GR_TOWER_STATUS_CONJECTURAL (GR_TOWER_STATUS_SCHANUEL for Lambert W).
    The caller is responsible for canonical arguments: values which are
    algebraic or expressible through other generators (Gamma(1), erf(0),
    ...) must not be adjoined, since the generator would contradict its
    status. Returns GR_DOMAIN at a pole or for u = 0 where the value is
    trivial (erf, Lambert W, polylogarithms), GR_UNABLE if the value
    cannot be evaluated numerically.
*/
int gr_tower_adjoin_special(gr_tower_t T, int kind, slong param, gr_srcptr u, const char * name);
int gr_tower_adjoin_special_flat(gr_tower_t T, int kind, slong param, const fmpz_mpoly_q_t u, const fmpz_mpoly_ctx_t mctx, const char * name);

/* Enclosure of a special function value (kind and param as above). */
int _gr_tower_special_eval(acb_t res, int kind, slong param, const acb_t u, slong prec);

/* Writes the definition of a special function value, with the argument
   given as a string (ignored for constants). */
int _gr_tower_special_write(gr_stream_t out, int kind, slong param, const char * arg);

/* Whether the special function is real at the real point u (u with an
   exactly zero imaginary part), certainly. */
int _gr_tower_special_real_at(int kind, slong param, const acb_t u, slong prec);

/* Index (1-based) of the transcendental generator with the given definition id, or 0. */
slong gr_tower_find_trans_def(const gr_tower_t T, ulong def_id);

/* Numerical evaluation */

/* Enclosure of the transcendental generator t_j (1 <= j <= num_trans). */
int gr_tower_trans_get_acb(acb_t res, gr_tower_t T, slong j, slong prec);

/* Enclosure of the generator a_k (1 <= k <= n), refined to precision prec. */
int gr_tower_step_get_acb(acb_t res, gr_tower_t T, slong k, slong prec);

/* Enclosure of an element x of F_k. */
int gr_tower_get_acb_at(acb_t res, gr_srcptr x, slong k, slong prec, gr_tower_t T);

GR_TOWER_INLINE int gr_tower_get_acb(acb_t res, gr_srcptr x, slong prec, gr_tower_t T)
{
    return gr_tower_get_acb_at(res, x, T->length, prec, T);
}

/* Dynamic evaluation */

/*
    Processes the zero divisors recorded by the quotient ring contexts,
    replacing each affected minimal polynomial by the factor vanishing at
    the corresponding enclosure. Returns nonzero if the tower changed.
*/
int gr_tower_refine(gr_tower_t T);

/* Refines step k of the tower using the factor g of m_k or m_k / g,
   whichever vanishes at the enclosure. */
int gr_tower_refine_step(gr_tower_t T, slong k, const gr_poly_t g);

/* Complete predicates and inversion in the top field, or in F_k. */
truth_t gr_tower_is_zero(gr_srcptr x, gr_tower_t T);
truth_t gr_tower_equal(gr_srcptr x, gr_srcptr y, gr_tower_t T);
int gr_tower_inv(gr_ptr res, gr_srcptr x, gr_tower_t T);
int gr_tower_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_tower_t T);
truth_t gr_tower_is_zero_at(gr_srcptr x, slong k, gr_tower_t T);
truth_t gr_tower_equal_at(gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T);
int gr_tower_inv_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T);
int gr_tower_div_at(gr_ptr res, gr_srcptr x, gr_srcptr y, slong k, gr_tower_t T);

/* Absolute representation */

/* Coordinates of an element of F_k in the monomial basis over F_0. */
int gr_tower_get_coeffs_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T);
int gr_tower_set_coeffs_at(gr_ptr res, gr_srcptr coeffs, slong k, gr_tower_t T);

/* Matrix of multiplication by x on F_n as an F_0-vector space. */
int gr_tower_multiplication_matrix(gr_mat_t res, gr_srcptr x, gr_tower_t T);

/* Characteristic polynomial of x over F_0 (annihilating polynomial of degree [F_n:F_0]). */
int gr_tower_charpoly(gr_poly_t res, gr_srcptr x, gr_tower_t T);

/* Minimal polynomial and qqbar of an element (base field QQ). */
int gr_tower_get_fmpz_poly_minpoly(fmpz_poly_t res, gr_srcptr x, gr_tower_t T);
int gr_tower_get_qqbar(qqbar_t res, gr_srcptr x, gr_tower_t T);

/* Promotes an element of F_j to an element of F_k (j <= k). */
int gr_tower_promote(gr_ptr res, gr_srcptr x, slong j, slong k, gr_tower_t T);

/* Expressing numbers in a tower */

/*
    Given a polynomial q over the top field and an enclosure z isolating
    one root of q, attempts to find an element res of the top field which
    is that root. Returns GR_SUCCESS if found (and verified exactly),
    GR_UNABLE if no representation was found within the precision limits.
    The tower may be refined as a side effect of exact verification.
*/
int gr_tower_express(gr_ptr res, const gr_poly_t q, const acb_t z, gr_tower_t T);

/* As above with an explicit limit on the working precision (bits);
   prec_limit <= 0 selects the default. */
int gr_tower_express_limit(gr_ptr res, const gr_poly_t q, const acb_t z, slong prec_limit, gr_tower_t T);

/* The principal square root of x in F_k / F_n, if it lies there: GR_DOMAIN
   if x is not a square in that field, GR_UNABLE if the tower has a step of
   degree other than 2 below k (the search is exact and needs no lattice
   reduction; it recurses through quadratic steps). */
int gr_tower_sqrt_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_t T);
int gr_tower_sqrt(gr_ptr res, gr_srcptr x, gr_tower_t T);

/* Expresses the algebraic number x in the tower (base field QQ).
   Returns GR_DOMAIN if x provably does not belong to the field. */
int gr_tower_express_qqbar(gr_ptr res, const qqbar_t x, gr_tower_t T);

/* Maps between towers */

/*
    A homomorphism from a source tower S to a target tower T over the
    same constant field, given by the images of the generators of S as
    flat elements of T. The map is defined on the first `length`
    generators of S in definition order. Images live in the context
    `mctx` of T's flat machinery (possibly a superseded one; call
    gr_tower_map_sync to bring them to the current context). Since flat
    elements survive extensions of the target, a map remains valid when
    the target grows.
*/
typedef struct
{
    gr_tower_struct * source;
    gr_tower_struct * target;
    slong length;                 /* number of generators (in definition order) of the source on which the map is defined */
    slong alloc;
    fmpz_mpoly_q_struct * images; /* images[d] = image of the generator with def_order d */
    fmpz_mpoly_ctx_struct * mctx;
}
gr_tower_map_struct;

typedef gr_tower_map_struct gr_tower_map_t[1];

/* Initializes a map defined on no generators (length 0); the images are
   filled in by gr_tower_map_set_inclusion or gr_tower_absorb, or set
   with gr_tower_map_set_image. */
void gr_tower_map_init(gr_tower_map_t map, gr_tower_t source, gr_tower_t target);
void gr_tower_map_clear(gr_tower_map_t map);
void gr_tower_map_set(gr_tower_map_t res, const gr_tower_map_t map);
int gr_tower_map_fit_length(gr_tower_map_t map, slong len);

/* Brings the images to the current flat context of the target. */
void gr_tower_map_sync(gr_tower_map_t map);

/* Image of the generator with definition order d (0 <= d < length),
   as a flat element in the context map->mctx. */
GR_TOWER_INLINE fmpz_mpoly_q_struct * gr_tower_map_image(gr_tower_map_t map, slong d)
{
    return map->images + d;
}

/* Sets the image of the generator with definition order d (extending
   length to d + 1 if needed) to the element x of F_k of the target
   (nested representation), or to a flat element in the given context. */
int gr_tower_map_set_image(gr_tower_map_t map, slong d, gr_srcptr x, slong k);
int gr_tower_map_set_image_flat(gr_tower_map_t map, slong d, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t mctx);

/* Sets the map to the inclusion of a tower into a tower which extends it
   (the generators of the source must be, in definition order, the first
   generators of the target). */
int gr_tower_map_set_inclusion(gr_tower_map_t map);

/* Applies the map to an element x of F_k of the source (using only
   generators covered by the map), respectively of the top field of the
   source. The result is an element of the top field of the target. */
int gr_tower_map_apply_at(gr_ptr res, gr_srcptr x, slong k, gr_tower_map_t map);
int gr_tower_map_apply(gr_ptr res, gr_srcptr x, gr_tower_map_t map);

/* Applies the map to a flat element x of the source (in the current or a
   superseded context of the source's flat machinery). The result is a
   flat element of the target, in the context map->mctx after
   gr_tower_map_sync (which this function calls first). */
int gr_tower_map_apply_flat(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, gr_tower_map_t map);

/* Applies the map to a polynomial over the top field of the source. */
int gr_tower_map_apply_poly(gr_poly_t res, const gr_poly_t f, gr_tower_map_t map);

/*
    Extends the tower U (which must currently be the target of map)
    so that it contains the tower B, and sets map to a homomorphism
    B -> U. Each generator of B is either expressed in U (if flags
    contains GR_TOWER_MERGE_EXPRESS and the search succeeds) or adjoined
    to U. Generators of B whose minimal polynomial has become linear are
    always substituted. Transcendental generators of B are identified with
    generators of U having the same definition id, and adjoined otherwise.
*/
#define GR_TOWER_MERGE_EXPRESS 1
#define GR_TOWER_MERGE_REAL_FIRST 8   /* a real algebraic generator expressible only through nonreal ones is inserted before them instead (real views) */
/* (the flags of the lazy fields, gr_tower_lazy.h, share this space:
   gr_ctx_init_tower_lazy takes both kinds in one argument) */
#define GR_TOWER_LAZY_REAL 2          /* the real subfield: operations with non-real results fail with GR_DOMAIN */
#define GR_TOWER_LAZY_ALGEBRAIC 4     /* the algebraic subfield: exp, log, pi, ... only at their algebraic special values */

/* Precision limit for expression searches during merges (they are
   heuristic fast paths; missing a representation is harmless): the
   option GR_TOWER_OPT_EXPRESS_PREC plus 4 bits per degree. */
#define GR_TOWER_MERGE_EXPRESS_PREC(T, D) (GR_TOWER_OPTION(T, GR_TOWER_OPT_EXPRESS_PREC) + 4 * (D))

int gr_tower_absorb(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, int flags);

/* As above, absorbing only the first p generators of B in definition
   order (all if p < 0). */
int gr_tower_absorb_prefix(gr_tower_t U, gr_tower_map_t map, gr_tower_t B, slong p, int flags);

/* Builds U = A + B with maps A -> U and B -> U. mapA and mapB are
   initialized by this function. */
int gr_tower_merge(gr_tower_t U, gr_tower_map_t mapA, gr_tower_map_t mapB, gr_tower_t A, gr_tower_t B, int flags);

/* Builds a tower U isomorphic to T in which steps of degree 1 have been
   removed, with the isomorphism map: T -> U. map is initialized by
   this function. */
int gr_tower_eliminate(gr_tower_t U, gr_tower_map_t map, gr_tower_t T);

/* A gr context for the top field of a tower with complete
   zero testing and dynamic refinement. */
void gr_ctx_init_tower_field(gr_ctx_t ctx, gr_tower_t T);
int gr_tower_field_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx);

/* Flat representation */

void gr_tower_flat_init(gr_tower_flat_t F, gr_tower_t T, slong cap);
void gr_tower_flat_clear(gr_tower_flat_t F);
int gr_tower_flat_ensure(gr_tower_flat_t F);
void gr_tower_flat_convert(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t old_mctx, gr_tower_flat_t F);
int gr_tower_flat_rationalize(fmpz_mpoly_q_t x, gr_tower_flat_t F);
int gr_tower_flat_has_alg_var(const fmpz_mpoly_t p, gr_tower_flat_t F);
slong gr_tower_flat_max_step(const fmpz_mpoly_t p, gr_tower_flat_t F);


int gr_tower_flat_reduce(fmpz_mpoly_q_t x, gr_tower_flat_t F);
int gr_tower_flat_set_nested_at(fmpz_mpoly_q_t res, gr_srcptr x, slong k, gr_tower_flat_t F);
int gr_tower_flat_poly_get_nested_at(gr_ptr res, const fmpz_mpoly_t f, slong k, gr_tower_flat_t F);
int gr_tower_flat_get_nested_at(gr_ptr res, const fmpz_mpoly_q_t x, slong k, gr_tower_flat_t F);
int gr_tower_flat_get_acb(acb_t res, const fmpz_mpoly_q_t x, slong prec, gr_tower_flat_t F);
/* Length of the shortest prefix (in definition order) of the tower
   containing all generators occurring in x (0 for a constant). */
slong gr_tower_flat_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* Highest algebraic generator occurring in x (0 if none). */
slong gr_tower_flat_alg_level(const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* Substitutes the flat elements imgs[i] (elements of F's tower) for the
   variables of the context x_mctx in x (imgs[i] may be NULL for a
   variable not occurring in x). The result is reduced in F. */
int gr_tower_flat_compose(fmpz_mpoly_q_t res, const fmpz_mpoly_q_t x, const fmpz_mpoly_ctx_t x_mctx, fmpz_mpoly_q_struct ** imgs, gr_tower_flat_t F);
truth_t gr_tower_flat_num_is_zero(fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* Highest definition order of a generator occurring in x (-1 if none). */
slong gr_tower_flat_def_order(const fmpz_mpoly_q_t x, gr_tower_flat_t F);

/* The flat machinery of the tower itself (initialized on first use). */
gr_tower_flat_struct * gr_tower_flat(gr_tower_t T);

/* A gr context for the top field of a (fixed) tower using the flat
   representation. Inversion is free; denominators are not canonical. */
void gr_ctx_init_tower_field_flat(gr_ctx_t ctx, gr_tower_t T);
int gr_tower_flat_set_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
int gr_tower_flat_get_nested(gr_ptr res, gr_srcptr x, gr_ctx_t ctx);
int gr_tower_field_flat_get_acb(acb_t res, gr_srcptr x, slong prec, gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
