/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CA_H
#define CA_H

#ifdef CA_INLINES_C
#define CA_INLINE
#else
#define CA_INLINE static inline
#endif

#include "flint.h"
#include "fmpz_mpoly_q.h"
#include "nf_elem.h"
#include "qqbar.h"
#include "fmpz_mpoly.h"
#include "fexpr.h"

#ifdef __cplusplus
extern "C" {
#endif

#define CA_INFO(ctx, args) do { if (ctx->options[CA_OPT_VERBOSE]) \
    do { flint_printf("%s:%d:  ", __FILE__, __LINE__); flint_printf args; } \
    while (0); } while (0);

/* Number object *************************************************************/

typedef union
{
    fmpq q;                           /* rational number */
    nf_elem_struct nf;                /* algebraic number field element */
    fmpz_mpoly_q_struct * mpoly_q;    /* generic field element */
}
ca_elem_struct;

typedef struct
{
    ulong field;
    ca_elem_struct elem;
}
ca_struct;


typedef ca_struct ca_t[1];
typedef ca_struct * ca_ptr;
typedef const ca_struct * ca_srcptr;

#define CA_FMPQ(x)         (&((x)->elem.q))
#define CA_MPOLY_Q(x)      (&(((x)->elem.mpoly_q)[0]))
#define CA_NF_ELEM(x)      (&((x)->elem.nf))
#define CA_FMPQ_NUMREF(x)  (fmpq_numref(CA_FMPQ(x)))
#define CA_FMPQ_DENREF(x)  (fmpq_denref(CA_FMPQ(x)))

#define CA_FIELD(x, ctx)     ((ca_field_ptr) ((x)->field))
#define CA_FIELD_ULONG(x)    ((x)->field)

/* We always allocate QQ and QQ(i) */
#define CA_IS_QQ(x, ctx) (CA_FIELD(x, ctx) == (ctx)->field_qq)
#define CA_IS_QQ_I(x, ctx) (CA_FIELD(x, ctx) == (ctx)->field_qq_i)

/* Use the low two bits of the field pointer to encode special values. */
/* The field pointer with the mask removed is NULL for
   Unknown/Undefined/Uinf, and a normal field pointer for signed
   infinity (encoding the sign). */
#define CA_UNKNOWN        UWORD(1)
#define CA_UNDEFINED      UWORD(2)
#define CA_INF            UWORD(3)
#define CA_SPECIAL        (CA_UNKNOWN | CA_UNDEFINED | CA_INF)

#define CA_IS_SPECIAL(x)       (CA_FIELD_ULONG(x) & CA_SPECIAL)
#define CA_IS_UNKNOWN(x)       (CA_FIELD_ULONG(x) == CA_UNKNOWN)
#define CA_IS_UNDEFINED(x)     (CA_FIELD_ULONG(x) == CA_UNDEFINED)
#define CA_IS_INF(x)           ((CA_FIELD_ULONG(x) & CA_SPECIAL) == CA_INF)
#define CA_IS_UNSIGNED_INF(x)  (CA_FIELD_ULONG(x) == CA_INF)
#define CA_IS_SIGNED_INF(x)    (CA_IS_INF(x) && !CA_IS_UNSIGNED_INF(x))

#define CA_FIELD_UNSPECIAL(x, ctx) ((ca_field_ptr) (CA_FIELD_ULONG(x) & ~CA_SPECIAL))



/* Extension object **********************************************************/

typedef struct
{
    qqbar_struct x;        /* qqbar_t element */
    nf_struct * nf;        /* antic number field for fast arithmetic */
}
ca_ext_qqbar;

typedef struct
{
    ca_struct * args;       /* Function arguments x1, ..., xn. */
    slong nargs;            /* Number of function arguments n. */
    acb_struct enclosure;   /* Cached numerical enclosure of f(x1,...,xn) */
    slong prec;             /* Working precision of cached enclosure */
    qqbar_struct * qqbar;   /* Cached qqbar */
}
ca_ext_func_data;

typedef struct
{
    calcium_func_code head;   /* f = F_Pi, F_Exp, ... */
    ulong hash;
    slong depth;
    union {
        ca_ext_qqbar qqbar;
        ca_ext_func_data func_data;
    } data;
} ca_ext_struct;

typedef ca_ext_struct ca_ext_t[1];
typedef ca_ext_struct * ca_ext_ptr;
typedef const ca_ext_struct * ca_ext_srcptr;

#define CA_EXT_HEAD(x) ((x)->head)
#define CA_EXT_HASH(x) ((x)->hash)
#define CA_EXT_DEPTH(x) ((x)->depth)

#define CA_EXT_IS_QQBAR(x) ((x)->head == CA_QQBar)

#define CA_EXT_QQBAR(_x) (&((_x)->data.qqbar.x))
#define CA_EXT_QQBAR_NF(_x) ((_x)->data.qqbar.nf)

#define CA_EXT_FUNC_ARGS(x) ((x)->data.func_data.args)
#define CA_EXT_FUNC_NARGS(x) ((x)->data.func_data.nargs)
#define CA_EXT_FUNC_ENCLOSURE(x) (&((x)->data.func_data.enclosure))
#define CA_EXT_FUNC_PREC(x) ((x)->data.func_data.prec)

typedef struct
{
    ca_ext_struct ** items;
    slong length;
    slong alloc;

    slong hash_size;
    slong * hash_table;
}
ca_ext_cache_struct;

typedef ca_ext_cache_struct ca_ext_cache_t[1];

/* Field object **************************************************************/

typedef struct
{
    slong length;                /* Number of generators              */
    ca_ext_struct ** ext;        /* Generators                        */
    fmpz_mpoly_vec_struct ideal; /* Algebraic relations for reduction */
    ulong hash;
}
ca_field_struct;

typedef ca_field_struct ca_field_t[1];
typedef ca_field_struct * ca_field_ptr;
typedef const ca_field_struct * ca_field_srcptr;

#define CA_FIELD_LENGTH(K) ((K)->length)
#define CA_FIELD_EXT(K) ((K)->ext)
#define CA_FIELD_EXT_ELEM(K, i) ((K)->ext[i])
#define CA_FIELD_HASH(K) ((K)->hash)

#define CA_FIELD_IS_QQ(K) ((K)->length == 0)
#define CA_FIELD_IS_NF(K) ((K)->ideal.length == -1)
#define CA_FIELD_IS_GENERIC(K) (!CA_FIELD_IS_QQ(K) && !CA_FIELD_IS_NF(K))

#define CA_FIELD_NF(K) (((K)->ext[0]->data.qqbar.nf))
#define CA_FIELD_NF_QQBAR(K) (&((K)->ext[0]->data.qqbar.x))

#define CA_FIELD_IDEAL(K) (&((K)->ideal))
#define CA_FIELD_IDEAL_ELEM(K, i) fmpz_mpoly_vec_entry(CA_FIELD_IDEAL(K), i)
#define CA_FIELD_IDEAL_LENGTH(K) ((K)->ideal.length)
#define CA_FIELD_IDEAL_ALLOC(K) ((K)->ideal.alloc)
#define CA_FIELD_IDEAL_P(K) ((K)->ideal.p)

#define CA_MCTX_1(ctx) ((ctx)->mctx[0])
#define CA_FIELD_MCTX(K, ctx) ((ctx)->mctx[CA_FIELD_LENGTH(K) - 1])

typedef struct
{
    ca_field_struct ** items;
    slong length;
    slong alloc;

    slong hash_size;
    slong * hash_table;
}
ca_field_cache_struct;

typedef ca_field_cache_struct ca_field_cache_t[1];

/* Context object ************************************************************/

enum
{
    CA_OPT_VERBOSE,
    CA_OPT_PRINT_FLAGS,
    CA_OPT_MPOLY_ORD,
    CA_OPT_PREC_LIMIT,
    CA_OPT_QQBAR_DEG_LIMIT,
    CA_OPT_LOW_PREC,
    CA_OPT_SMOOTH_LIMIT,
    CA_OPT_LLL_PREC,
    CA_OPT_POW_LIMIT,
    CA_OPT_USE_GROEBNER,
    CA_OPT_GROEBNER_LENGTH_LIMIT,
    CA_OPT_GROEBNER_POLY_LENGTH_LIMIT,
    CA_OPT_GROEBNER_POLY_BITS_LIMIT,
    CA_OPT_VIETA_LIMIT,
    CA_OPT_TRIG_FORM,
    CA_OPT_NUM_OPTIONS
};

#define CA_TRIG_DIRECT       0
#define CA_TRIG_EXPONENTIAL  1
#define CA_TRIG_SINE_COSINE  2
#define CA_TRIG_TANGENT      3

typedef struct
{
    ca_ext_cache_struct ext_cache;              /* Cached extension objects */
    ca_field_cache_struct field_cache;          /* Cached extension fields  */
    ca_field_struct * field_qq;                 /* Quick access to QQ      */
    ca_field_struct * field_qq_i;               /* Quick access to QQ(i)   */
    fmpz_mpoly_ctx_struct ** mctx;              /* Cached contexts for multivariate polys */
    slong mctx_len;
    slong * options;
}
ca_ctx_struct;

typedef ca_ctx_struct ca_ctx_t[1];

#define CA_CTX_EXT_CACHE(ctx) (&((ctx)->ext_cache))
#define CA_CTX_FIELD_CACHE(ctx) (&((ctx)->field_cache))

#define CA_CTX_FIELD_WITH_INDEX(ctx, i) ((&((ctx)->field_cache))->items[i])

/* Context management */

void ca_ctx_init(ca_ctx_t ctx);
void ca_ctx_clear(ca_ctx_t ctx);
void ca_ctx_print(ca_ctx_t ctx);

CA_INLINE slong ca_ctx_get_option(ca_ctx_t ctx, slong i)
{
    return ctx->options[i];
}

CA_INLINE void ca_ctx_set_option(ca_ctx_t ctx, slong i, slong value)
{
    ctx->options[i] = value;
}

ca_field_ptr _ca_ctx_get_field_const(ca_ctx_t ctx, calcium_func_code func);
ca_field_ptr _ca_ctx_get_field_fx(ca_ctx_t ctx, calcium_func_code func, const ca_t x);
ca_field_ptr _ca_ctx_get_field_fxy(ca_ctx_t ctx, calcium_func_code func, const ca_t x, const ca_t y);

/* Numbers */

void ca_init(ca_t x, ca_ctx_t ctx);
void ca_clear(ca_t x, ca_ctx_t ctx);
void ca_swap(ca_t x, ca_t y, ca_ctx_t ctx);
void _ca_make_field_element(ca_t x, ca_field_srcptr field, ca_ctx_t ctx);

CA_INLINE void
_ca_make_fmpq(ca_t x, ca_ctx_t ctx)
{
    if (!CA_IS_QQ(x, ctx))
        _ca_make_field_element(x, ctx->field_qq, ctx);
}

CA_INLINE void
_ca_function_fx(ca_t res, calcium_func_code func, const ca_t x, ca_ctx_t ctx)
{
    ca_field_srcptr field = _ca_ctx_get_field_fx(ctx, func, x);
    _ca_make_field_element(res, field, ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_FIELD_MCTX(field, ctx));
}

CA_INLINE void
_ca_function_fxy(ca_t res, calcium_func_code func, const ca_t x, const ca_t y, ca_ctx_t ctx)
{
    ca_field_srcptr field = _ca_ctx_get_field_fxy(ctx, func, x, y);
    _ca_make_field_element(res, field, ctx);
    fmpz_mpoly_q_gen(CA_MPOLY_Q(res), 0, CA_FIELD_MCTX(field, ctx));
}

void ca_set(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_transfer(ca_t res, ca_ctx_t res_ctx, const ca_t src, ca_ctx_t src_ctx);

void ca_zero(ca_t x, ca_ctx_t ctx);
void ca_one(ca_t x, ca_ctx_t ctx);
void ca_neg_one(ca_t x, ca_ctx_t ctx);

void ca_set_si(ca_t x, slong v, ca_ctx_t ctx);
void ca_set_ui(ca_t x, ulong v, ca_ctx_t ctx);
void ca_set_fmpz(ca_t x, const fmpz_t v, ca_ctx_t ctx);
void ca_set_fmpq(ca_t x, const fmpq_t v, ca_ctx_t ctx);
void ca_set_d(ca_t res, double x, ca_ctx_t ctx);
void ca_set_d_d(ca_t res, double x, double y, ca_ctx_t ctx);

void ca_i(ca_t x, ca_ctx_t ctx);
void ca_neg_i(ca_t x, ca_ctx_t ctx);
void ca_pi(ca_t res, ca_ctx_t ctx);
void ca_pi_i(ca_t res, ca_ctx_t ctx);
void ca_euler(ca_t res, ca_ctx_t ctx);

void ca_unknown(ca_t x, ca_ctx_t ctx);

void ca_undefined(ca_t x, ca_ctx_t ctx);
void ca_uinf(ca_t x, ca_ctx_t ctx);
void ca_pos_inf(ca_t x, ca_ctx_t ctx);
void ca_neg_inf(ca_t x, ca_ctx_t ctx);
void ca_pos_i_inf(ca_t x, ca_ctx_t ctx);
void ca_neg_i_inf(ca_t x, ca_ctx_t ctx);

void ca_set_qqbar(ca_t res, const qqbar_t x, ca_ctx_t ctx);
int ca_can_evaluate_qqbar(const ca_t x, ca_ctx_t ctx);
int ca_get_qqbar(qqbar_t res, const ca_t x, ca_ctx_t ctx);
int ca_get_fmpq(fmpq_t res, const ca_t x, ca_ctx_t ctx);
int ca_get_fmpz(fmpz_t res, const ca_t x, ca_ctx_t ctx);

/* Symbolic expressions */

#define CA_FEXPR_SERIALIZATION 1

void ca_get_fexpr(fexpr_t res, const ca_t x, ulong flags, ca_ctx_t ctx);
int ca_set_fexpr(ca_t res, const fexpr_t expr, ca_ctx_t ctx);

/* Printing */

#define CA_PRINT_N         UWORD(1)
#define CA_PRINT_REPR      UWORD(2)
#define CA_PRINT_FIELD     UWORD(4)
#define CA_PRINT_DIGITS    UWORD(16)
#define CA_PRINT_DEFAULT   (CA_PRINT_N | CA_PRINT_REPR)
#define CA_PRINT_DEBUG     (CA_PRINT_N | CA_PRINT_REPR | CA_PRINT_FIELD)

#ifdef FLINT_HAVE_FILE
void ca_fprint(FILE * fp, const ca_t x, ca_ctx_t ctx);
#endif

void ca_print(const ca_t x, ca_ctx_t ctx);
void ca_printn(const ca_t x, slong n, ca_ctx_t ctx);
char * ca_get_str(const ca_t x, ca_ctx_t ctx);

/* Random generation */

void ca_randtest_same_nf(ca_t res, flint_rand_t state, const ca_t x, slong bits, slong den_bits, ca_ctx_t ctx);
void ca_randtest_rational(ca_t res, flint_rand_t state, slong bits, ca_ctx_t ctx);
void ca_randtest(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx);
void ca_randtest_special(ca_t res, flint_rand_t state, slong depth, slong bits, ca_ctx_t ctx);

/* Representation properties */

CA_INLINE int ca_is_special(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_SPECIAL(x);
}

CA_INLINE int ca_is_unknown(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_UNKNOWN(x);
}

CA_INLINE int ca_is_qq_elem(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_QQ(x, ctx);
}

CA_INLINE int ca_is_qq_elem_zero(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_QQ(x, ctx) && fmpq_is_zero(CA_FMPQ(x));
}

CA_INLINE int ca_is_qq_elem_one(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_QQ(x, ctx) && fmpq_is_one(CA_FMPQ(x));
}

CA_INLINE int ca_is_qq_elem_integer(const ca_t x, ca_ctx_t ctx)
{
    return CA_IS_QQ(x, ctx) && fmpz_is_one(CA_FMPQ_DENREF(x));
}

CA_INLINE int ca_is_nf_elem(const ca_t x, ca_ctx_t ctx)
{
    return !CA_IS_SPECIAL(x) && CA_FIELD_IS_NF(CA_FIELD(x, ctx));
}

CA_INLINE int ca_is_generic_elem(const ca_t x, ca_ctx_t ctx)
{
    return !CA_IS_SPECIAL(x) && CA_FIELD_IS_GENERIC(CA_FIELD(x, ctx));
}

int
ca_is_cyclotomic_nf_elem(slong * p, ulong * q, const ca_t x, ca_ctx_t ctx);

/* Value predicates and comparisons */

truth_t ca_is_zero_check_fast(const ca_t x, ca_ctx_t ctx);


truth_t ca_check_is_number(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_zero(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_one(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_neg_one(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_i(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_neg_i(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_algebraic(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_rational(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_integer(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_real(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_negative_real(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_imaginary(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_undefined(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_infinity(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_uinf(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_signed_inf(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_pos_inf(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_neg_inf(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_pos_i_inf(const ca_t x, ca_ctx_t ctx);
truth_t ca_check_is_neg_i_inf(const ca_t x, ca_ctx_t ctx);

truth_t ca_check_equal(const ca_t x, const ca_t y, ca_ctx_t ctx);
truth_t ca_check_lt(const ca_t x, const ca_t y, ca_ctx_t ctx);
truth_t ca_check_le(const ca_t x, const ca_t y, ca_ctx_t ctx);
truth_t ca_check_gt(const ca_t x, const ca_t y, ca_ctx_t ctx);
truth_t ca_check_ge(const ca_t x, const ca_t y, ca_ctx_t ctx);

int ca_equal_repr(const ca_t x, const ca_t y, ca_ctx_t ctx);
int ca_cmp_repr(const ca_t x, const ca_t y, ca_ctx_t ctx);
ulong ca_hash_repr(const ca_t x, ca_ctx_t ctx);

/* Field structure operations */

void ca_merge_fields(ca_t resx, ca_t resy, const ca_t x, const ca_t y, ca_ctx_t ctx);
void ca_condense_field(ca_t res, ca_ctx_t ctx);
ca_ext_ptr ca_is_gen_as_ext(const ca_t x, ca_ctx_t ctx);

/* Arithmetic */

/* todo: document */
void _ca_mpoly_q_reduce_ideal(fmpz_mpoly_q_t res, ca_field_srcptr field, ca_ctx_t ctx);
void _ca_mpoly_q_simplify_fraction_ideal(fmpz_mpoly_q_t res, ca_field_srcptr field, ca_ctx_t ctx);


void ca_neg(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_add_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx);
void ca_add_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx);
void ca_add_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx);
void ca_add_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx);
void ca_add(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx);

void ca_sub_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx);
void ca_sub_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx);
void ca_sub_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx);
void ca_sub_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx);
void ca_fmpq_sub(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx);
void ca_fmpz_sub(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx);
void ca_ui_sub(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx);
void ca_si_sub(ca_t res, slong x, const ca_t y, ca_ctx_t ctx);
void ca_sub(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx);

void ca_mul_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx);
void ca_mul_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx);
void ca_mul_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx);
void ca_mul_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx);
void ca_mul(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx);

void ca_inv(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_fmpq_div(ca_t res, const fmpq_t x, const ca_t y, ca_ctx_t ctx);
void ca_fmpz_div(ca_t res, const fmpz_t x, const ca_t y, ca_ctx_t ctx);
void ca_ui_div(ca_t res, ulong x, const ca_t y, ca_ctx_t ctx);
void ca_si_div(ca_t res, slong x, const ca_t y, ca_ctx_t ctx);

void ca_div_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx);
void ca_div_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx);
void ca_div_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx);
void ca_div_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx);
void ca_div(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx);

void ca_dot(ca_t res, const ca_t initial, int subtract,
    ca_srcptr x, slong xstep, ca_srcptr y, slong ystep, slong len, ca_ctx_t ctx);

void ca_fmpz_poly_evaluate(ca_t res, const fmpz_poly_t poly, const ca_t x, ca_ctx_t ctx);
void ca_fmpq_poly_evaluate(ca_t res, const fmpq_poly_t poly, const ca_t x, ca_ctx_t ctx);

void ca_fmpz_mpoly_evaluate_horner(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);
void ca_fmpz_mpoly_evaluate_iter(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);
void ca_fmpz_mpoly_evaluate(ca_t res, const fmpz_mpoly_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);

void ca_fmpz_mpoly_q_evaluate(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);

void ca_fmpz_mpoly_q_evaluate_no_division_by_zero(ca_t res, const fmpz_mpoly_q_t f, ca_srcptr x, const fmpz_mpoly_ctx_t mctx, ca_ctx_t ctx);
void ca_inv_no_division_by_zero(ca_t res, const ca_t x, ca_ctx_t ctx);


/* Powers and roots */

CA_INLINE void
ca_sqr(ca_t res, const ca_t x, ca_ctx_t ctx)
{
    ca_mul(res, x, x, ctx);
}

void ca_pow_fmpq(ca_t res, const ca_t x, const fmpq_t y, ca_ctx_t ctx);
void ca_pow_fmpz(ca_t res, const ca_t x, const fmpz_t y, ca_ctx_t ctx);
void ca_pow_ui(ca_t res, const ca_t x, ulong y, ca_ctx_t ctx);
void ca_pow_si(ca_t res, const ca_t x, slong y, ca_ctx_t ctx);
void ca_pow(ca_t res, const ca_t x, const ca_t y, ca_ctx_t ctx);

void ca_pow_si_arithmetic(ca_t res, const ca_t x, slong n, ca_ctx_t ctx);

void ca_sqrt_inert(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_sqrt_nofactor(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_sqrt_factor(ca_t res, const ca_t x, ulong flags, ca_ctx_t ctx);
void ca_sqrt(ca_t res, const ca_t x, ca_ctx_t ctx);

CA_INLINE void
ca_sqrt_ui(ca_t res, ulong n, ca_ctx_t ctx)
{
    ca_set_ui(res, n, ctx);
    ca_sqrt(res, res, ctx);
}

/* Complex parts */

void ca_conj_shallow(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_conj_deep(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_conj(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_abs(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_sgn(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_csgn(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_arg(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_re(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_im(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_floor(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_ceil(ca_t res, const ca_t x, ca_ctx_t ctx);

/* Elementary functions */

void ca_exp(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_log(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_sin_cos_exponential(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx);
void ca_sin_cos_direct_exp_hack(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx);
void ca_sin_cos_direct(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx);
void ca_sin_cos_tangent(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx);
void ca_sin_cos(ca_t res1, ca_t res2, const ca_t x, ca_ctx_t ctx);
void ca_sin(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_cos(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_tan_sine_cosine(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_tan_exponential(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_tan_direct(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_tan(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_cot(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_atan_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_atan_direct(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_atan(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_asin_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_asin_direct(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_asin(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_acos_logarithm(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_acos_direct(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_acos(ca_t res, const ca_t x, ca_ctx_t ctx);

/* Special functions */

void ca_erf(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_erfc(ca_t res, const ca_t x, ca_ctx_t ctx);
void ca_erfi(ca_t res, const ca_t x, ca_ctx_t ctx);

void ca_gamma(ca_t res, const ca_t x, ca_ctx_t ctx);

/* Numerical evaluation */

void ca_get_acb_raw(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx);
void ca_get_acb(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx);
void ca_get_acb_accurate_parts(acb_t res, const ca_t x, slong prec, ca_ctx_t ctx);

char * ca_get_decimal_str(const ca_t x, slong digits, ulong flags, ca_ctx_t ctx);

/* Factorisation */

#define CA_FACTOR_ZZ_NONE        0
#define CA_FACTOR_ZZ_SMOOTH      2
#define CA_FACTOR_ZZ_FULL        4
#define CA_FACTOR_POLY_NONE      0
#define CA_FACTOR_POLY_CONTENT   64
#define CA_FACTOR_POLY_SQF       128
#define CA_FACTOR_POLY_FULL      256

typedef struct
{
    ca_ptr base;
    ca_ptr exp;
    slong length;
    slong alloc;
}
ca_factor_struct;

typedef ca_factor_struct ca_factor_t[1];

void ca_factor_init(ca_factor_t fac, ca_ctx_t ctx);
void ca_factor_clear(ca_factor_t fac, ca_ctx_t ctx);
void ca_factor_one(ca_factor_t fac, ca_ctx_t ctx);
void ca_factor_print(const ca_factor_t fac, ca_ctx_t ctx);
void ca_factor_insert(ca_factor_t fac, const ca_t base, const ca_t exp, ca_ctx_t ctx);
void ca_factor_get_ca(ca_t res, const ca_factor_t fac, ca_ctx_t ctx);

void ca_factor(ca_factor_t res, const ca_t x, ulong flags, ca_ctx_t ctx);

/* Test helpers */

#define CA_TEST_PROPERTY(f, s, x, ctx, expected) \
     do { \
        truth_t t; \
        t = f(x, ctx); \
        if (t != expected) \
        { \
            flint_printf("FAIL\n"); \
            flint_printf("%s\n\n", s); \
            flint_printf("x = "); ca_print(x, ctx); flint_printf("\n\n"); \
            flint_printf("got = "); truth_print(t); flint_printf("\n\n"); \
            flint_printf("expected = "); truth_print(expected); flint_printf("\n\n"); \
            flint_abort(); \
        } \
     } while (0) \

#ifdef __cplusplus
}
#endif

#endif
