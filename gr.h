#ifndef GR_H
#define GR_H

#ifdef GR_INLINES_C
#define GR_INLINE FLINT_DLL
#else
#define GR_INLINE static __inline__
#endif

#include <string.h>
#include <assert.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/fmpz_poly.h"
#include "flint/fmpq_poly.h"
#include "flint/fmpz_mpoly.h"

#ifdef __cplusplus
 extern "C" {
#endif

#ifndef CALCIUM_H

typedef enum
{
    T_TRUE,
    T_FALSE,
    T_UNKNOWN
} truth_t;

#endif

GR_INLINE truth_t truth_and(truth_t x, truth_t y)
{
    if (x == T_FALSE || y == T_FALSE)
        return T_FALSE;
    if (x == T_TRUE && y == T_TRUE)
        return T_TRUE;
    return T_UNKNOWN;
}

GR_INLINE truth_t truth_or(truth_t x, truth_t y)
{
    if (x == T_TRUE || y == T_TRUE)
        return T_TRUE;
    if (x == T_FALSE && y == T_FALSE)
        return T_FALSE;
    return T_UNKNOWN;
}

GR_INLINE truth_t truth_not(truth_t x)
{
    if (x == T_TRUE)
        return T_FALSE;
    if (x == T_FALSE)
        return T_TRUE;
    return T_UNKNOWN;
}

GR_INLINE void truth_println(truth_t x)
{
    if (x == T_TRUE) flint_printf("T_TRUE\n");
    if (x == T_FALSE) flint_printf("T_FALSE\n");
    if (x == T_UNKNOWN) flint_printf("T_UNKNOWN\n");
}

typedef int (*gr_funcptr)(void);

/* Copied from Calcium: stream interface allows simple file or string IO. */

typedef struct
{
    FILE * fp;
    char * s;
    slong len;
    slong alloc;
}
gr_stream_struct;

typedef gr_stream_struct gr_stream_t[1];

void gr_stream_init_file(gr_stream_t out, FILE * fp);
void gr_stream_init_str(gr_stream_t out);
void gr_stream_write(gr_stream_t out, const char * s);
void gr_stream_write_si(gr_stream_t out, slong x);
void gr_stream_write_ui(gr_stream_t out, ulong x);
void gr_stream_write_free(gr_stream_t out, char * s);
void gr_stream_write_fmpz(gr_stream_t out, const fmpz_t x);

#define GR_SUCCESS 0
#define GR_DOMAIN 1
#define GR_UNABLE 2
#define GR_TEST_FAIL 4

#define WARN_UNUSED_RESULT __attribute__((warn_unused_result))

#define GR_MUST_SUCCEED(expr) do { if ((expr) != GR_SUCCESS) { flint_printf("GR_MUST_SUCCEED failure: %s", __FILE__); flint_abort(); } } while (0)

typedef void * gr_ptr;
typedef const void * gr_srcptr;
typedef void * gr_ctx_ptr;

#define GR_ENTRY(vec, i, size) ((void *) (((char *) (vec)) + ((i) * (size))))

GR_INLINE int gr_not_implemented(void) { return GR_UNABLE; }

typedef enum
{
    GR_METHOD_CTX_WRITE,
    GR_METHOD_CTX_CLEAR,

    /* general ring properties */
    GR_METHOD_CTX_IS_RING,
    GR_METHOD_CTX_IS_COMMUTATIVE_RING,
    GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,
    GR_METHOD_CTX_IS_FIELD,
    GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
    GR_METHOD_CTX_IS_FINITE,
    GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
    GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,

    /* ring properties related to orderings and norms */
    GR_METHOD_CTX_IS_ORDERED_RING,

    /* group properties */
    GR_METHOD_CTX_IS_MULTIPLICATIVE_GROUP,

    /* context properties represented to the representation */
    GR_METHOD_CTX_IS_EXACT,            /* we have no inexact elements */
    GR_METHOD_CTX_IS_CANONICAL,        /* we have no non-canonical representations */

    GR_METHOD_INIT,
    GR_METHOD_CLEAR,
    GR_METHOD_SWAP,
    GR_METHOD_RANDTEST,
    GR_METHOD_WRITE,

    GR_METHOD_ZERO,
    GR_METHOD_ONE,
    GR_METHOD_NEG_ONE,

    GR_METHOD_IS_ZERO,
    GR_METHOD_IS_ONE,
    GR_METHOD_IS_NEG_ONE,

    GR_METHOD_EQUAL,

    GR_METHOD_SET,
    GR_METHOD_SET_UI,
    GR_METHOD_SET_SI,
    GR_METHOD_SET_FMPZ,
    GR_METHOD_SET_FMPQ,
    GR_METHOD_SET_OTHER,
    GR_METHOD_SET_STR,

    GR_METHOD_GET_SI,
    GR_METHOD_GET_UI,
    GR_METHOD_GET_FMPZ,
    GR_METHOD_GET_FMPQ,

    GR_METHOD_IS_INTEGER,
    GR_METHOD_IS_RATIONAL,

    GR_METHOD_NEG,

    GR_METHOD_ADD,
    GR_METHOD_ADD_UI,
    GR_METHOD_ADD_SI,
    GR_METHOD_ADD_FMPZ,
    GR_METHOD_ADD_FMPQ,

    GR_METHOD_SUB,
    GR_METHOD_SUB_UI,
    GR_METHOD_SUB_SI,
    GR_METHOD_SUB_FMPZ,
    GR_METHOD_SUB_FMPQ,

    GR_METHOD_MUL,
    GR_METHOD_MUL_UI,
    GR_METHOD_MUL_SI,
    GR_METHOD_MUL_FMPZ,
    GR_METHOD_MUL_FMPQ,

    GR_METHOD_ADDMUL,
    GR_METHOD_ADDMUL_UI,
    GR_METHOD_ADDMUL_SI,
    GR_METHOD_ADDMUL_FMPZ,
    GR_METHOD_ADDMUL_FMPQ,

    GR_METHOD_SUBMUL,
    GR_METHOD_SUBMUL_UI,
    GR_METHOD_SUBMUL_SI,
    GR_METHOD_SUBMUL_FMPZ,
    GR_METHOD_SUBMUL_FMPQ,

    GR_METHOD_MUL_TWO,
    GR_METHOD_SQR,

    GR_METHOD_IS_INVERTIBLE,
    GR_METHOD_INV,

    GR_METHOD_DIV,
    GR_METHOD_DIV_UI,
    GR_METHOD_DIV_SI,
    GR_METHOD_DIV_FMPZ,
    GR_METHOD_DIV_FMPQ,

    GR_METHOD_DIVEXACT,
    GR_METHOD_DIVEXACT_UI,
    GR_METHOD_DIVEXACT_SI,
    GR_METHOD_DIVEXACT_FMPZ,
    GR_METHOD_DIVEXACT_FMPQ,

    GR_METHOD_POW_UI,
    GR_METHOD_POW_SI,
    GR_METHOD_POW_FMPZ,
    GR_METHOD_POW_FMPQ,
    GR_METHOD_POW,

    GR_METHOD_IS_SQUARE,
    GR_METHOD_SQRT,
    GR_METHOD_RSQRT,

    /* todo: test the following */
    GR_METHOD_ABS,
    GR_METHOD_CONJ,
    GR_METHOD_RE,
    GR_METHOD_IM,
    GR_METHOD_SGN,  /* todo: implement */
    GR_METHOD_CSGN, /* todo: implement */

    GR_METHOD_CMP,
    GR_METHOD_CMPABS,

    GR_METHOD_EXP,
    GR_METHOD_LOG,
    GR_METHOD_SIN,
    GR_METHOD_COS,
    GR_METHOD_ATAN,

    /* Finite field methods */
    GR_METHOD_CTX_FQ_PRIME,
    GR_METHOD_CTX_FQ_DEGREE,
    GR_METHOD_CTX_FQ_ORDER,
    GR_METHOD_FQ_GEN,
    GR_METHOD_FQ_FROBENIUS,
    GR_METHOD_FQ_MULTIPLICATIVE_ORDER,
    GR_METHOD_FQ_NORM,
    GR_METHOD_FQ_TRACE,
    GR_METHOD_FQ_IS_PRIMITIVE,
    GR_METHOD_FQ_PTH_ROOT,

    /* Vector methods */
    GR_METHOD_VEC_INIT,
    GR_METHOD_VEC_CLEAR,
    GR_METHOD_VEC_SWAP,
    GR_METHOD_VEC_SET,
    GR_METHOD_VEC_ZERO,
    GR_METHOD_VEC_EQUAL,
    GR_METHOD_VEC_IS_ZERO,
    GR_METHOD_VEC_NEG,
    GR_METHOD_VEC_ADD,
    GR_METHOD_VEC_SUB,
    GR_METHOD_VEC_SCALAR_MUL,
    GR_METHOD_VEC_SCALAR_ADDMUL,
    GR_METHOD_VEC_SCALAR_SUBMUL,
    GR_METHOD_VEC_SCALAR_ADDMUL_SI,
    GR_METHOD_VEC_SCALAR_SUBMUL_SI,

    GR_METHOD_VEC_DOT,
    GR_METHOD_VEC_DOT_REV,

    GR_METHOD_VEC_DOT_UI,
    GR_METHOD_VEC_DOT_SI,
    GR_METHOD_VEC_DOT_FMPZ,

    GR_METHOD_VEC_SET_POWERS,

    /* Polynomial methods */
    GR_METHOD_POLY_MULLOW,

    /* Matrix methods */
    GR_METHOD_MAT_MUL,

    GR_METHOD_TAB_SIZE
}
gr_method;

typedef gr_funcptr gr_static_method_table[GR_METHOD_TAB_SIZE];

typedef struct
{
    gr_method index;
    gr_funcptr function;
}
gr_method_tab_input;

void gr_method_tab_init(gr_funcptr * methods, gr_method_tab_input * tab);

/* Identify specific rings/fields. */

typedef enum
{
    GR_CTX_FMPZ,
    GR_CTX_FMPQ,

    GR_CTX_FMPZ_MOD,
    GR_CTX_NMOD8,

    GR_CTX_FQ,
    GR_CTX_FQ_NMOD,
    GR_CTX_FQ_ZECH,

    GR_CTX_REAL_ALGEBRAIC_QQBAR,
    GR_CTX_COMPLEX_ALGEBRAIC_QQBAR,

    GR_CTX_REAL_ALGEBRAIC_CA,
    GR_CTX_COMPLEX_ALGEBRAIC_CA,

    GR_CTX_RR_CA,
    GR_CTX_CC_CA,

    GR_CTX_RR_ARB,
    GR_CTX_CC_ACB,

    GR_CTX_GR_POLY,
    GR_CTX_GR_MPOLY,
    GR_CTX_GR_MAT,

    GR_CTX_PSL2Z,

    GR_CTX_WHICH_STRUCTURE_TAB_SIZE
}
gr_which_structure;

typedef struct
{
    ulong flags;
    ulong which_ring;
    ssize_t sizeof_elem;
    void * elem_ctx;
    gr_funcptr * methods;
    ulong size_limit;
}
gr_ctx_struct;

typedef gr_ctx_struct gr_ctx_t[1];

GR_INLINE slong gr_ctx_sizeof_elem(gr_ctx_t ctx)
{
    return ctx->sizeof_elem;
}

/* Typedefs for method function pointers. */

typedef void ((*gr_method_init_clear_op)(gr_ptr, gr_ctx_ptr));
typedef void ((*gr_method_swap_op)(gr_ptr, gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_ctx)(gr_ctx_ptr));
typedef truth_t ((*gr_method_ctx_predicate)(gr_ctx_ptr));
typedef int ((*gr_method_ctx_stream)(gr_stream_t, gr_ctx_ptr));
typedef int ((*gr_method_stream_in)(gr_stream_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_randtest)(gr_ptr, flint_rand_t state, gr_ctx_ptr));
typedef int ((*gr_method_constant_op)(gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_constant_op_get_si)(slong *, gr_ctx_ptr));
typedef int ((*gr_method_constant_op_get_fmpz)(fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_si)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_ui)(gr_ptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpz)(gr_ptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpq)(gr_ptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_other)(gr_ptr, gr_srcptr, gr_ctx_ptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_str)(gr_ptr, const char *, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_fmpz)(fmpz_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_get_int)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef truth_t ((*gr_method_unary_predicate)(gr_srcptr, gr_ctx_ptr));
typedef truth_t ((*gr_method_binary_predicate)(gr_srcptr, gr_srcptr, gr_ctx_ptr));

typedef void ((*gr_method_vec_init_clear_op)(gr_ptr, slong, gr_ctx_ptr));
typedef void ((*gr_method_vec_swap_op)(gr_ptr, gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_constant_op)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_op)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_si)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr));
typedef truth_t ((*gr_method_vec_predicate)(gr_srcptr, slong, gr_ctx_ptr));
typedef truth_t ((*gr_method_vec_vec_predicate)(gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_op)(gr_ptr, gr_srcptr, int, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_si_op)(gr_ptr, gr_srcptr, int, gr_srcptr, const slong *, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_ui_op)(gr_ptr, gr_srcptr, int, gr_srcptr, const ulong *, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_fmpz_op)(gr_ptr, gr_srcptr, int, gr_srcptr, const fmpz *, slong, gr_ctx_ptr));

typedef int ((*gr_method_poly_binary_trunc_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, slong, slong, gr_ctx_ptr));


/* Macros to retrieve methods (with correct call signature) from context object. */

#define GR_CTX_OP(ctx, NAME) (((gr_method_ctx *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_STREAM(ctx, NAME) (((gr_method_ctx_stream *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CTX_PREDICATE(ctx, NAME) (((gr_method_ctx_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_STREAM_IN(ctx, NAME) (((gr_method_stream_in *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_RANDTEST(ctx, NAME) (((gr_method_randtest *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_INIT_CLEAR_OP(ctx, NAME) (((gr_method_init_clear_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_SWAP_OP(ctx, NAME) (((gr_method_swap_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP(ctx, NAME) (((gr_method_constant_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP_GET_SI(ctx, NAME) (((gr_method_constant_op_get_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP_GET_FMPZ(ctx, NAME) (((gr_method_constant_op_get_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP(ctx, NAME) (((gr_method_unary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_SI(ctx, NAME) (((gr_method_unary_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_UI(ctx, NAME) (((gr_method_unary_op_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPZ(ctx, NAME) (((gr_method_unary_op_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPQ(ctx, NAME) (((gr_method_unary_op_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_OTHER(ctx, NAME) (((gr_method_unary_op_other *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_STR(ctx, NAME) (((gr_method_unary_op_str *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_FMPZ(ctx, NAME) (((gr_method_unary_op_get_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP(ctx, NAME) (((gr_method_binary_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_SI(ctx, NAME) (((gr_method_binary_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_UI(ctx, NAME) (((gr_method_binary_op_ui *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPZ(ctx, NAME) (((gr_method_binary_op_fmpz *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPQ(ctx, NAME) (((gr_method_binary_op_fmpq *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_GET_INT(ctx, NAME) (((gr_method_binary_op_get_int *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_PREDICATE(ctx, NAME) (((gr_method_unary_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_PREDICATE(ctx, NAME) (((gr_method_binary_predicate *) ctx->methods)[GR_METHOD_ ## NAME])

#define GR_VEC_INIT_CLEAR_OP(ctx, NAME) (((gr_method_vec_init_clear_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SWAP_OP(ctx, NAME) (((gr_method_vec_swap_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_CONSTANT_OP(ctx, NAME) (((gr_method_vec_constant_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_OP(ctx, NAME) (((gr_method_vec_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_OP(ctx, NAME) (((gr_method_vec_vec_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP(ctx, NAME) (((gr_method_vec_scalar_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_SI(ctx, NAME) (((gr_method_vec_scalar_op_si *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_vec_predicate *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_DOT_OP(ctx, NAME) (((gr_method_vec_dot_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_DOT_SI_OP(ctx, NAME) (((gr_method_vec_dot_si_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_DOT_UI_OP(ctx, NAME) (((gr_method_vec_dot_ui_op *) ctx->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_DOT_FMPZ_OP(ctx, NAME) (((gr_method_vec_dot_fmpz_op *) ctx->methods)[GR_METHOD_ ## NAME])

#define GR_POLY_BINARY_TRUNC_OP(ctx, NAME) (((gr_method_poly_binary_trunc_op *) ctx->methods)[GR_METHOD_ ## NAME])


/* Wrappers to call methods. */

GR_INLINE int gr_ctx_clear(gr_ctx_t ctx) { return GR_CTX_OP(ctx, CTX_CLEAR)(ctx); }
GR_INLINE int gr_ctx_write(gr_stream_t out, gr_ctx_t ctx) { return GR_CTX_STREAM(ctx, CTX_WRITE)(out, ctx); }

GR_INLINE truth_t gr_ctx_is_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_RING)(ctx); }
GR_INLINE truth_t gr_ctx_is_commutative_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_COMMUTATIVE_RING)(ctx); }
GR_INLINE truth_t gr_ctx_is_integral_domain(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_INTEGRAL_DOMAIN)(ctx); }
GR_INLINE truth_t gr_ctx_is_field(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FIELD)(ctx); }

GR_INLINE truth_t gr_ctx_is_unique_factorization_domain(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_UNIQUE_FACTORIZATION_DOMAIN)(ctx); }
GR_INLINE truth_t gr_ctx_is_finite(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FINITE)(ctx); }
GR_INLINE truth_t gr_ctx_is_finite_characteristic(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_FINITE_CHARACTERISTIC)(ctx); }
GR_INLINE truth_t gr_ctx_is_algebraically_closed(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_ALGEBRAICALLY_CLOSED)(ctx); }
GR_INLINE truth_t gr_ctx_is_ordered_ring(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_ORDERED_RING)(ctx); }

GR_INLINE truth_t gr_ctx_is_multiplicative_group(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_MULTIPLICATIVE_GROUP)(ctx); }

GR_INLINE truth_t gr_ctx_is_exact(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_EXACT)(ctx); }
GR_INLINE truth_t gr_ctx_is_canonical(gr_ctx_t ctx) { return GR_CTX_PREDICATE(ctx, CTX_IS_CANONICAL)(ctx); }

GR_INLINE void gr_init(gr_ptr res, gr_ctx_t ctx) { GR_INIT_CLEAR_OP(ctx, INIT)(res, ctx); }
GR_INLINE void gr_clear(gr_ptr res, gr_ctx_t ctx) { GR_INIT_CLEAR_OP(ctx, CLEAR)(res, ctx); }
GR_INLINE void gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) { GR_SWAP_OP(ctx, SWAP)(x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_randtest(gr_ptr x, flint_rand_t state, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST)(x, state, ctx); }
GR_INLINE /* todo: warn? */ int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) { return GR_STREAM_IN(ctx, WRITE)(out, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_zero(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ZERO)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ONE)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_neg_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, NEG_ONE)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SET)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_SI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_ui(gr_ptr res, ulong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_UI)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fmpz(gr_ptr res, const fmpz_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPZ(ctx, SET_FMPZ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_fmpq(gr_ptr res, const fmpq_t x, gr_ctx_t ctx) { return GR_UNARY_OP_FMPQ(ctx, SET_FMPQ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_other(gr_ptr res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx) { return GR_UNARY_OP_OTHER(ctx, SET_OTHER)(res, x, x_ctx, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_set_str(gr_ptr res, const char * x, gr_ctx_t ctx) { return GR_UNARY_OP_STR(ctx, SET_STR)(res, x, ctx); }

GR_INLINE truth_t gr_is_zero(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ZERO)(x, ctx); }
GR_INLINE truth_t gr_is_one(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ONE)(x, ctx); }
GR_INLINE truth_t gr_is_neg_one(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_NEG_ONE)(x, ctx); }

GR_INLINE truth_t gr_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_PREDICATE(ctx, EQUAL)(x, y, ctx); }
GR_INLINE truth_t gr_not_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return truth_not(GR_BINARY_PREDICATE(ctx, EQUAL)(x, y, ctx)); }

GR_INLINE WARN_UNUSED_RESULT int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NEG)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADD)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, ADD_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADD_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, ADD_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_add_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, ADD_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUB)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, SUB_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, SUB_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, SUB_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sub_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, SUB_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, MUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, MUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, MUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, MUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_mul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, MUL_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_addmul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADDMUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, ADDMUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADDMUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, ADDMUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_addmul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, ADDMUL_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_submul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUBMUL)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, SUBMUL_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, SUBMUL_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, SUBMUL_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_submul_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, SUBMUL_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_mul_two(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, MUL_TWO)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sqr(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SQR)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIV)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, DIV_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, DIV_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, DIV_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_div_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, DIV_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, INV)(res, x, ctx); }
GR_INLINE truth_t gr_is_invertible(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_INVERTIBLE)(x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, POW)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, POW_UI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, POW_SI)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_fmpz(gr_ptr res, gr_srcptr x, const fmpz_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPZ(ctx, POW_FMPZ)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_pow_fmpq(gr_ptr res, gr_srcptr x, const fmpq_t y, gr_ctx_t ctx) { return GR_BINARY_OP_FMPQ(ctx, POW_FMPQ)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SQRT)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RSQRT)(res, x, ctx); }
GR_INLINE truth_t gr_is_square(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_SQUARE)(x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_abs(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ABS)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_conj(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CONJ)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_re(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, RE)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_im(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, IM)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SGN)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_csgn(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, CSGN)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_cmp(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP_GET_INT(ctx, CMP)(res, x, y, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_cmpabs(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP_GET_INT(ctx, CMPABS)(res, x, y, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_exp(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, EXP)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_log(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, LOG)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_sin(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SIN)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_cos(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, COS)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_atan(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, ATAN)(res, x, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_prime(fmpz_t res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_FMPZ(ctx, CTX_FQ_PRIME)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_degree(slong * res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_SI(ctx, CTX_FQ_DEGREE)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_ctx_fq_order(fmpz_t res, gr_ctx_t ctx) { return GR_CONSTANT_OP_GET_FMPZ(ctx, CTX_FQ_ORDER)(res, ctx); }

GR_INLINE WARN_UNUSED_RESULT int gr_fq_gen(fmpz_t res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, FQ_GEN)(res, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, FQ_FROBENIUS)(res, x, e, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_MULTIPLICATIVE_ORDER)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_NORM)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP_GET_FMPZ(ctx, FQ_TRACE)(res, x, ctx); }
GR_INLINE WARN_UNUSED_RESULT truth_t gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, FQ_IS_PRIMITIVE)(x, ctx); }
GR_INLINE WARN_UNUSED_RESULT int gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, FQ_PTH_ROOT)(res, x, ctx); }

GR_INLINE void _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx) { GR_VEC_INIT_CLEAR_OP(ctx, VEC_INIT)(vec, len, ctx); }
GR_INLINE void _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx) { GR_VEC_INIT_CLEAR_OP(ctx, VEC_CLEAR)(vec, len, ctx); }
GR_INLINE void _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx) { GR_VEC_SWAP_OP(ctx, VEC_SWAP)(vec1, vec2, len, ctx); }

GR_INLINE WARN_UNUSED_RESULT int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_ZERO)(vec, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_SET)(res, src, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_NEG)(res, src, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_ADD)(res, src1, src2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_SUB)(res, src1, src2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_scalar_mul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_MUL)(vec1, vec2, len, c, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_ADDMUL)(vec1, vec2, len, c, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_SUBMUL)(vec1, vec2, len, c, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_ADDMUL_SI)(vec1, vec2, len, c, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_SUBMUL_SI)(vec1, vec2, len, c, ctx); }
GR_INLINE truth_t _gr_vec_equal(gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_PREDICATE(ctx, VEC_EQUAL)(vec1, vec2, len, ctx); }
GR_INLINE truth_t _gr_vec_is_zero(gr_srcptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_PREDICATE(ctx, VEC_IS_ZERO)(vec, len, ctx); }

GR_INLINE WARN_UNUSED_RESULT int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT_REV)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_dot_si(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const slong * vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_SI_OP(ctx, VEC_DOT_SI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_dot_ui(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const ulong * vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_UI_OP(ctx, VEC_DOT_UI)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_INLINE WARN_UNUSED_RESULT int _gr_vec_dot_fmpz(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, const fmpz * vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_FMPZ_OP(ctx, VEC_DOT_FMPZ)(res, initial, subtract, vec1, vec2, len, ctx); }

GR_INLINE WARN_UNUSED_RESULT int _gr_vec_set_powers(gr_ptr res, gr_srcptr x, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_SET_POWERS)(res, x, len, ctx); }


GR_INLINE WARN_UNUSED_RESULT int
_gr_poly_mullow(gr_ptr res, gr_srcptr poly1, slong len1, gr_srcptr poly2, slong len2, slong len, gr_ctx_t ctx)
{
    return GR_POLY_BINARY_TRUNC_OP(ctx, POLY_MULLOW)(res, poly1, len1, poly2, len2, len, ctx);
}

/* some useful generic functions, currently not overloadable */
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate_horner(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate_rectangular(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int _gr_fmpz_poly_evaluate(gr_ptr res, const fmpz * f, slong len, gr_srcptr x, gr_ctx_t ctx);
WARN_UNUSED_RESULT int gr_fmpz_poly_evaluate(gr_ptr res, const fmpz_poly_t f, gr_srcptr x, gr_ctx_t ctx);

WARN_UNUSED_RESULT int gr_fmpz_mpoly_evaluate(gr_ptr res, const fmpz_mpoly_t f, gr_srcptr x, const fmpz_mpoly_ctx_t mctx, gr_ctx_t ctx);

/* todo: could allow overloading this as well */
/* todo: worth warning about unused result? */
GR_INLINE int
_gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;

    status = GR_SUCCESS;
    for (i = 0; i < len; i++)
    {
        if (n_randint(state, 2))
            status |= gr_zero(GR_ENTRY(res, i, sz), ctx);
        else
            status |= gr_randtest(GR_ENTRY(res, i, sz), state, ctx);
    }

    return status;
}

/* todo: warn unused? */
int gr_ctx_print(gr_ctx_t ctx);
int gr_ctx_println(gr_ctx_t ctx);
int gr_print(gr_srcptr x, gr_ctx_t ctx);
int gr_println(gr_srcptr x, gr_ctx_t ctx);
int gr_ctx_get_str(char ** s, gr_ctx_t ctx);
int gr_get_str(char ** s, gr_srcptr x, gr_ctx_t ctx);

/* Macros for allocating temporary variables on the stack. */
/* todo: use vector init/clear functions when provided */

#define GR_TMP_VEC_ALLOC_MAX_STACK 1024
#define GR_TMP_ALLOC(size) ((size <= GR_TMP_VEC_ALLOC_MAX_STACK) ? alloca(size) : flint_malloc(size))
#define GR_TMP_FREE(ptr, size) do { if (size > GR_TMP_VEC_ALLOC_MAX_STACK) flint_free(ptr); } while (0)
#define GR_TMP_ALLOC_SMALL(size) alloca(size)

#define GR_TMP_INIT_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_init_clear_op vec_init = GR_VEC_INIT_CLEAR_OP(ctx, VEC_INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        (vec) = (gr_ptr) GR_TMP_ALLOC((len) * _gr_elem_size); \
        vec_init((vec), (len), (ctx)); \
    } while (0)

#define GR_TMP_CLEAR_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_init_clear_op vec_clear = GR_VEC_INIT_CLEAR_OP(ctx, VEC_CLEAR); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        vec_clear((vec), (len), (ctx)); \
        GR_TMP_FREE(vec, (len) * _gr_elem_size); \
    } while (0)

#define GR_TMP_INIT(x1, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(1 * _gr_elem_size); \
        init(x1, (ctx)); \
    } while (0)

#define GR_TMP_INIT2(x1, x2, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(2 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
    } while (0)

#define GR_TMP_INIT3(x1, x2, x3, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(3 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
    } while (0)

#define GR_TMP_INIT4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(4 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        x4 = (gr_ptr) ((char *) x3 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
        init(x4, (ctx)); \
    } while (0)

#define GR_TMP_INIT5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_init_clear_op init = GR_INIT_CLEAR_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC_SMALL(5 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        x4 = (gr_ptr) ((char *) x3 + _gr_elem_size); \
        x5 = (gr_ptr) ((char *) x4 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
        init(x4, (ctx)); \
        init(x5, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR(x1, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR2(x1, x2, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR3(x1, x2, x3, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_init_clear_op clear = GR_INIT_CLEAR_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
        clear(x5, (ctx)); \
    } while (0)

GR_INLINE gr_ptr gr_heap_init(gr_ctx_t ctx)
{
    gr_ptr ptr;
    ptr = (gr_ptr) flint_malloc(ctx->sizeof_elem);
    gr_init(ptr, ctx);
    return ptr;
}

GR_INLINE void gr_heap_clear(gr_ptr x, gr_ctx_t ctx)
{
    gr_clear(x, ctx);
    flint_free(x);
}

/* Some generic implementations */

truth_t gr_generic_ctx_predicate(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_true(gr_ctx_t ctx);
truth_t gr_generic_ctx_predicate_false(gr_ctx_t ctx);

/* Some base rings */

void gr_ctx_init_random(gr_ctx_t ctx, flint_rand_t state);

void gr_ctx_init_fmpz(gr_ctx_t ctx);
void gr_ctx_init_fmpq(gr_ctx_t ctx);

void gr_ctx_init_fmpz_mod(gr_ctx_t ctx, const fmpz_t n);
void gr_ctx_fmpz_mod_set_primality(gr_ctx_t ctx, truth_t is_prime);

void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n);

void gr_ctx_init_real_qqbar(gr_ctx_t ctx);
void gr_ctx_init_complex_qqbar(gr_ctx_t ctx);

void gr_ctx_init_real_arb(gr_ctx_t ctx, slong prec);
void gr_ctx_init_complex_acb(gr_ctx_t ctx, slong prec);

void gr_ctx_arb_set_prec(gr_ctx_t ctx, slong prec);
slong gr_ctx_arb_get_prec(gr_ctx_t ctx);


void gr_ctx_init_real_ca(gr_ctx_t ctx);
void gr_ctx_init_complex_ca(gr_ctx_t ctx);
void gr_ctx_init_real_algebraic_ca(gr_ctx_t ctx);
void gr_ctx_init_complex_algebraic_ca(gr_ctx_t ctx);

void gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);
void gr_ctx_init_fq_nmod(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);
void gr_ctx_init_fq_zech(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var);

/* Groups */

void gr_ctx_init_psl2z(gr_ctx_t ctx);

/* Generic polynomial ring */

typedef struct
{
    gr_ctx_struct * base_ring;
    slong degree_limit;
    char * var;
}
polynomial_ctx_t;

#define POLYNOMIAL_CTX(ring_ctx) ((polynomial_ctx_t *)((ring_ctx)->elem_ctx))
#define POLYNOMIAL_ELEM_CTX(ring_ctx) (POLYNOMIAL_CTX(ring_ctx)->base_ring)

void gr_ctx_init_polynomial(gr_ctx_t ctx, gr_ctx_t base_ring);

/* Generic matrix ring */

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
}
matrix_ctx_t;

#define MATRIX_CTX(ring_ctx) ((matrix_ctx_t *)((ring_ctx)->elem_ctx))

void gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n);

/* Multivariate */

void gr_ctx_init_mpoly(gr_ctx_t ctx, gr_ctx_t base_ring, slong nvars, const ordering_t ord);

/* Coercions */

int gr_ctx_cmp_coercion(gr_ctx_t ctx1, gr_ctx_t ctx2);

/* Testing */

#define GR_TEST_VERBOSE 8
#define GR_TEST_ALWAYS_ABLE 16

/* todo: just have gr_test_structure() */
void gr_test_ring(gr_ctx_t R, slong iters, int test_flags);
void gr_test_multiplicative_group(gr_ctx_t R, slong iters, int test_flags);

#ifdef __cplusplus
}
#endif

#endif
