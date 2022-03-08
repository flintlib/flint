#ifndef GR_H
#define GR_H

#ifdef GR_INLINES_C
#define GR_INLINE FLINT_DLL
#else
#define GR_INLINE static __inline__
#endif

#include <string.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/nmod_vec.h"
#include "flint/ulong_extras.h"
#include "flint/profiler.h"

#ifdef __cplusplus
 extern "C" {
#endif

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

GR_INLINE void gr_stream_init_file(gr_stream_t out, FILE * fp)
{
    out->fp = fp;
    out->s = NULL;
}

GR_INLINE void gr_stream_init_str(gr_stream_t out)
{
    out->fp = NULL;
    out->s = flint_malloc(16);
    out->s[0] = '\0';
    out->len = 0;
    out->alloc = 16;
}

GR_INLINE void gr_stream_write(gr_stream_t out, const char * s)
{
    if (out->fp != NULL)
    {
        fprintf(out->fp, "%s", s);
    }
    else
    {
        slong len, alloc;

        len = strlen(s);
        alloc = out->len + len + 1;

        if (alloc > out->alloc)
        {
            alloc = FLINT_MAX(alloc, out->alloc * 2);
            out->s = realloc(out->s, alloc);
            out->alloc = alloc;
        }

        memcpy(out->s + out->len, s, len + 1);
        out->len += len;
    }
}

GR_INLINE void
gr_stream_write_si(gr_stream_t out, slong x)
{
    if (out->fp != NULL)
    {
        flint_fprintf(out->fp, "%wd", x);
    }
    else
    {
        char tmp[22];
        if (sizeof(slong) == sizeof(long))
            sprintf(tmp, "%ld", x);
        else
            flint_sprintf(tmp, "%wd", x);
        gr_stream_write(out, tmp);
    }
}

GR_INLINE void
gr_stream_write_free(gr_stream_t out, char * s)
{
    gr_stream_write(out, s);
    flint_free(s);
}

/*
All functions implement error handling and return a status code.

GR_SUCCESS  - The operation finished as expected.

GR_DOMAIN   - Domain error: the input is not valid for this operations
              (e.g. division by zero, square root of a non-square) or the
              result does not fit in the range of the output type
              (e.g. get_ui() with input >= 2^64).

GR_UNABLE   - The operation could not be performed for implementation reasons
              (e.g. too large input, missing algorithm). If this flag is
              set, there is also potentially a domain error (but this is
              unknown).

GR_WRONG    - Test failure (used only in test code).

For uniformity, all functions return a status code, even functions
that should never fail (in which case we might throw in an assert).
Flags can be OR'ed to avoid complex control flow.
*/

#define GR_SUCCESS 0
#define GR_DOMAIN 1
#define GR_UNABLE 2
#define GR_WRONG 4

typedef void * gr_ptr;
typedef const void * gr_srcptr;
typedef void * gr_ctx_ptr;

GR_INLINE int gr_not_implemented(void) { return GR_UNABLE; }

typedef enum
{
    GR_METHOD_INIT,
    GR_METHOD_CLEAR,
    GR_METHOD_SWAP,
    GR_METHOD_RANDTEST,
    GR_METHOD_WRITE,

    GR_METHOD_ZERO,
    GR_METHOD_ONE,
    GR_METHOD_IS_ZERO,
    GR_METHOD_IS_ONE,

    GR_METHOD_EQUAL,

    GR_METHOD_SET,
    GR_METHOD_SET_SI,
    GR_METHOD_SET_UI,
    GR_METHOD_SET_FMPZ,
    GR_METHOD_SET_FMPQ,

    GR_METHOD_GET_SI,

    GR_METHOD_NEG,

    GR_METHOD_ADD,
    GR_METHOD_ADD_SI,
    GR_METHOD_ADD_UI,
    GR_METHOD_ADD_FMPZ,
    GR_METHOD_ADD_FMPQ,

    GR_METHOD_SUB,
    GR_METHOD_SUB_SI,
    GR_METHOD_SUB_UI,
    GR_METHOD_SUB_FMPZ,
    GR_METHOD_SUB_FMPQ,

    GR_METHOD_MUL,
    GR_METHOD_MUL_SI,
    GR_METHOD_MUL_UI,
    GR_METHOD_MUL_FMPZ,
    GR_METHOD_MUL_FMPQ,

    GR_METHOD_IS_INVERTIBLE,
    GR_METHOD_INV,

    GR_METHOD_DIV,
    GR_METHOD_DIV_SI,
    GR_METHOD_DIV_UI,
    GR_METHOD_DIV_FMPZ,
    GR_METHOD_DIV_FMPQ,

    GR_METHOD_DIVEXACT,
    GR_METHOD_DIVEXACT_SI,
    GR_METHOD_DIVEXACT_UI,
    GR_METHOD_DIVEXACT_FMPZ,
    GR_METHOD_DIVEXACT_FMPQ,

    GR_METHOD_POW_UI,
    GR_METHOD_POW_SI,
    GR_METHOD_POW_FMPZ,
    GR_METHOD_POW_FMPQ,
    GR_METHOD_POW,

    GR_METHOD_IS_SQUARE,
    GR_METHOD_SQRT,

    GR_METHOD_EXP,
    GR_METHOD_LOG,

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
    GR_METHOD_VEC_SCALAR_ADDMUL,
    GR_METHOD_VEC_SCALAR_SUBMUL,
    GR_METHOD_VEC_SCALAR_ADDMUL_SI,
    GR_METHOD_VEC_SCALAR_SUBMUL_SI,
    GR_METHOD_VEC_DOT,
    GR_METHOD_VEC_DOT_REV,

    GR_METHOD_TAB_SIZE,
}
gr_method;




typedef struct
{
    gr_funcptr * methods;
}
gr_method_tab_t;

typedef struct
{
    gr_method index;
    gr_funcptr function;
}
gr_method_tab_input;

void gr_method_tab_init(gr_method_tab_t * methods, gr_method_tab_input * tab);

/* Flags. These are not cumulative. (?) */

#define GR_COMMUTATIVE_RING  UWORD(1)
#define GR_INTEGRAL_DOMAIN   UWORD(2)
#define GR_FIELD             UWORD(4)

typedef struct
{
    ulong flags;
    ssize_t sizeof_elem;
    void * elem_ctx;
    char * debug_string;
    gr_method_tab_t * methods2;
}
gr_ctx_struct;

typedef gr_ctx_struct gr_ctx_t[1];


/* Typedefs for method function pointers. */

typedef int ((*gr_method_stream_in)(gr_stream_t, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_randtest)(gr_ptr, flint_rand_t state, const void * options, gr_ctx_ptr));
typedef int ((*gr_method_constant_op)(gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_swap_op)(gr_ptr, gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op)(gr_ptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_si)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_ui)(gr_ptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpz)(gr_ptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_fmpq)(gr_ptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_op_get_si)(slong * gr_ptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op)(gr_ptr, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_si)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_ui)(gr_ptr, gr_srcptr, ulong, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpz)(gr_ptr, gr_srcptr, const fmpz_t, gr_ctx_ptr));
typedef int ((*gr_method_binary_op_fmpq)(gr_ptr, gr_srcptr, const fmpq_t, gr_ctx_ptr));
typedef int ((*gr_method_unary_predicate)(int *, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_binary_predicate)(int *, gr_srcptr, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_vec_constant_op)(gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_swap_op)(gr_ptr, gr_ptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_op)(gr_ptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_vec_op)(gr_ptr, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op)(gr_ptr, gr_srcptr, slong, gr_srcptr, gr_ctx_ptr));
typedef int ((*gr_method_vec_scalar_op_si)(gr_ptr, gr_srcptr, slong, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_predicate)(int *, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_vec_predicate)(int *, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));
typedef int ((*gr_method_vec_dot_op)(gr_ptr, gr_srcptr, int, gr_srcptr, gr_srcptr, slong, gr_ctx_ptr));

/* Macros to retrieve methods (with correct call signature) from context object. */

#define GR_STREAM_IN(ctx, NAME) (((gr_method_stream_in *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_RANDTEST(ctx, NAME) (((gr_method_randtest *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_SWAP_OP(ctx, NAME) (((gr_method_swap_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_CONSTANT_OP(ctx, NAME) (((gr_method_constant_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP(ctx, NAME) (((gr_method_unary_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_SI(ctx, NAME) (((gr_method_unary_op_si *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_UI(ctx, NAME) (((gr_method_unary_op_ui *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPZ(ctx, NAME) (((gr_method_unary_op_fmpz *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_FMPQ(ctx, NAME) (((gr_method_unary_op_fmpq *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_OP_GET_SI(ctx, NAME) (((gr_method_unary_op_get_si *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP(ctx, NAME) (((gr_method_binary_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_SI(ctx, NAME) (((gr_method_binary_op_si *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_UI(ctx, NAME) (((gr_method_binary_op_ui *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPZ(ctx, NAME) (((gr_method_binary_op_fmpz *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_OP_FMPQ(ctx, NAME) (((gr_method_binary_op_fmpq *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_UNARY_PREDICATE(ctx, NAME) (((gr_method_unary_predicate *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_BINARY_PREDICATE(ctx, NAME) (((gr_method_binary_predicate *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_CONSTANT_OP(ctx, NAME) (((gr_method_vec_constant_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SWAP_OP(ctx, NAME) (((gr_method_vec_swap_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_OP(ctx, NAME) (((gr_method_vec_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_OP(ctx, NAME) (((gr_method_vec_vec_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP(ctx, NAME) (((gr_method_vec_scalar_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_SCALAR_OP_SI(ctx, NAME) (((gr_method_vec_scalar_op_si *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_predicate *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_VEC_PREDICATE(ctx, NAME) (((gr_method_vec_vec_predicate *) ctx->methods2->methods)[GR_METHOD_ ## NAME])
#define GR_VEC_DOT_OP(ctx, NAME) (((gr_method_vec_dot_op *) ctx->methods2->methods)[GR_METHOD_ ## NAME])

/* Wrappers to call methods. */

GR_INLINE int gr_init(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, INIT)(res, ctx); }
GR_INLINE int gr_clear(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, CLEAR)(res, ctx); }
GR_INLINE int gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) { return GR_SWAP_OP(ctx, SWAP)(x, y, ctx); }
GR_INLINE int gr_randtest(gr_ptr x, flint_rand_t state, const void * options, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST)(x, state, options, ctx); }
GR_INLINE int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) { return GR_STREAM_IN(ctx, WRITE)(out, x, ctx); }
GR_INLINE int gr_zero(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ZERO)(res, ctx); }
GR_INLINE int gr_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ONE)(res, ctx); }
GR_INLINE int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_SI)(res, x, ctx); }
GR_INLINE int gr_is_zero(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ZERO)(res, x, ctx); }
GR_INLINE int gr_is_one(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ONE)(res, x, ctx); }
GR_INLINE int gr_equal(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_PREDICATE(ctx, EQUAL)(res, x, y, ctx); }
GR_INLINE int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SET)(res, x, ctx); }
GR_INLINE int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NEG)(res, x, ctx); }
GR_INLINE int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADD)(res, x, y, ctx); }
GR_INLINE int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADD_SI)(res, x, y, ctx); }
GR_INLINE int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUB)(res, x, y, ctx); }
GR_INLINE int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, MUL)(res, x, y, ctx); }
GR_INLINE int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, MUL_SI)(res, x, y, ctx); }
GR_INLINE int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIV)(res, x, y, ctx); }
GR_INLINE int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, INV)(res, x, ctx); }
GR_INLINE int gr_is_invertible(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_INVERTIBLE)(res, x, ctx); }
GR_INLINE int gr_pow_ui(gr_ptr res, gr_srcptr x, ulong y, gr_ctx_t ctx) { return GR_BINARY_OP_UI(ctx, POW_UI)(res, x, y, ctx); }
GR_INLINE int gr_pow_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, POW_SI)(res, x, y, ctx); }


/* Macros for allocating temporary variables on the stack. */
/* todo: use vector init functions when provided */

#define GR_TMP_START TMP_INIT; TMP_START;
#define GR_TMP_ALLOC TMP_ALLOC
#define GR_TMP_END TMP_END

#define GR_TMP_INIT_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_constant_op vec_init = GR_VEC_CONSTANT_OP(ctx, VEC_INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        (vec) = (gr_ptr) GR_TMP_ALLOC((len) * _gr_elem_size); \
        vec_init((vec), (len), (ctx)); \
    } while (0)

#define GR_TMP_CLEAR_VEC(vec, len, ctx) \
    do { \
        gr_method_vec_constant_op vec_clear = GR_VEC_CONSTANT_OP(ctx, VEC_CLEAR); \
        vec_clear((vec), (len), (ctx)); \
    } while (0)

#define GR_TMP_INIT1(x1, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(1 * _gr_elem_size); \
        init(x1, (ctx)); \
    } while (0)

#define GR_TMP_INIT2(x1, x2, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(2 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
    } while (0)

#define GR_TMP_INIT3(x1, x2, x3, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(3 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        init(x1, (ctx)); \
        init(x2, (ctx)); \
        init(x3, (ctx)); \
    } while (0)

#define GR_TMP_INIT4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(4 * _gr_elem_size); \
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
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(5 * _gr_elem_size); \
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

#define GR_TMP_CLEAR1(x1, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR2(x1, x2, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR3(x1, x2, x3, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
    } while (0)

#define GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)); \
        clear(x2, (ctx)); \
        clear(x3, (ctx)); \
        clear(x4, (ctx)); \
        clear(x5, (ctx)); \
    } while (0)

#define GR_ENTRY(vec, i, size) ((void *) (((char *) (vec)) + ((i) * (size))))

GR_INLINE int _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_INIT)(vec, len, ctx); }
GR_INLINE int _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_CLEAR)(vec, len, ctx); }
GR_INLINE int _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_SWAP_OP(ctx, VEC_SWAP)(vec1, vec2, len, ctx); }
GR_INLINE int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_ZERO)(vec, len, ctx); }
GR_INLINE int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_SET)(res, src, len, ctx); }
GR_INLINE int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_NEG)(res, src, len, ctx); }
GR_INLINE int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_ADD)(res, src1, src2, len, ctx); }
GR_INLINE int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_SUB)(res, src1, src2, len, ctx); }
GR_INLINE int _gr_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_ADDMUL)(vec1, vec2, len, c, ctx); }
GR_INLINE int _gr_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_SUBMUL)(vec1, vec2, len, c, ctx); }
GR_INLINE int _gr_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_ADDMUL_SI)(vec1, vec2, len, c, ctx); }
GR_INLINE int _gr_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_SUBMUL_SI)(vec1, vec2, len, c, ctx); }
GR_INLINE int _gr_vec_equal(int * res, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_PREDICATE(ctx, VEC_EQUAL)(res, vec1, vec2, len, ctx); }
GR_INLINE int _gr_vec_is_zero(int * res, gr_srcptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_PREDICATE(ctx, VEC_IS_ZERO)(res, vec, len, ctx); }
GR_INLINE int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }
GR_INLINE int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT_REV)(res, initial, subtract, vec1, vec2, len, ctx); }

/* todo: could allow overloading this as well */
GR_INLINE int
_gr_vec_randtest(gr_ptr res, flint_rand_t state, slong len, void * options, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;

    status = GR_SUCCESS;
    for (i = 0; i < len; i++)
    {
        status |= gr_randtest(GR_ENTRY(res, i, sz), state, options, ctx);
    }

    return status;
}

GR_INLINE int
gr_print(gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_write(out, x, ctx);
    return GR_SUCCESS;
}

GR_INLINE int
gr_println(gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_write(out, x, ctx);
    gr_stream_write(out, "\n");
    return GR_SUCCESS;
}

/* Todo: generic vectors */

typedef struct
{
    gr_ptr entries;
    slong length;
    slong alloc;
}
gr_vec_struct;

typedef gr_vec_struct gr_vec_t[1];


/* Todo: generic matrices */

typedef struct
{
    gr_ptr entries;
    slong r;
    slong c;
    gr_ptr * rows;
}
gr_mat_struct;

typedef gr_mat_struct gr_mat_t[1];


#define GR_MAT_ENTRY(mat,i,j,ctx) GR_ENTRY((mat)->rows[i], j, sz)
#define gr_mat_nrows(mat, ctx) ((mat)->r)
#define gr_mat_ncols(mat, ctx) ((mat)->c)

GR_INLINE int
gr_mat_init(gr_mat_t mat, slong rows, slong cols, gr_ctx_t ctx)
{
    slong i, sz;

    sz = ctx->sizeof_elem;

    if (rows != 0)
        mat->rows = flint_malloc(rows * sizeof(gr_ptr));
    else
        mat->rows = NULL;

    if (rows != 0 && cols != 0)
    {
        mat->entries = (gr_ptr) flint_malloc(flint_mul_sizes(rows, cols) * sz);

        _gr_vec_init(mat->entries, rows * cols, ctx);

        for (i = 0; i < rows; i++)
            mat->rows[i] = GR_ENTRY(mat->entries, i * cols, sz);
    }
    else
    {
        mat->entries = NULL;
        for (i = 0; i < rows; i++)
            mat->rows[i] = NULL;
    }

    mat->r = rows;
    mat->c = cols;

    return GR_SUCCESS;
}

GR_INLINE int
gr_mat_clear(gr_mat_t mat, gr_ctx_t ctx)
{
    if (mat->entries != NULL)
    {
        _gr_vec_clear(mat->entries, mat->r * mat->c, ctx);

        flint_free(mat->entries);
        flint_free(mat->rows);
    }
    else if (mat->r != 0)
    {
        flint_free(mat->rows);
    }

    return GR_SUCCESS;
}

GR_INLINE int
gr_mat_randtest(gr_mat_t mat, flint_rand_t state, void * options, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    status = GR_SUCCESS;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_randtest(mat->rows[i], state, c, options, ctx);
    }

    return status;
}

GR_INLINE int
gr_mat_equal(int * res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status, equal, this_equal;
    slong i, r, c;

    status = GR_SUCCESS;

    r = gr_mat_nrows(mat1, ctx);
    c = gr_mat_ncols(mat1, ctx);

    if (r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx))
    {
        res[0] = 0;
        return GR_SUCCESS;
    }

    if (r == 0 || c == 0)
    {
        res[0] = 1;
        return GR_SUCCESS;
    }

    equal = 1;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_equal(&this_equal, mat1->rows[i], mat2->rows[i], c, ctx);
        equal = equal && this_equal;
    }

    res[0] = equal;
    return status;
}

GR_INLINE int
gr_mat_zero(gr_mat_t res, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    status = GR_SUCCESS;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_zero(res->rows[i], c, ctx);
    }

    return status;
}

GR_INLINE int
gr_mat_set_si(gr_mat_t res, slong v, gr_ctx_t ctx)
{
    int status;
    slong i, r, c, sz;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);
    sz = ctx->sizeof_elem;

    status = gr_mat_zero(res, ctx);

    for (i = 0; i < FLINT_MIN(r, c); i++)
        status |= gr_set_si(GR_MAT_ENTRY(res, i, i, sz), v, ctx);

    return status;
}

GR_INLINE int
gr_mat_one(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_set_si(res, 1, ctx);
}

GR_INLINE int
gr_mat_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (r != gr_mat_nrows(mat, ctx) ||
        c != gr_mat_ncols(mat, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;

    if (res != mat)
    {
        for (i = 0; i < r; i++)
        {
            status |= _gr_vec_set(res->rows[i], mat->rows[i], c, ctx);
        }
    }

    return status;
}

GR_INLINE int
gr_mat_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (r != gr_mat_nrows(mat, ctx) ||
        c != gr_mat_ncols(mat, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;

    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_neg(res->rows[i], mat->rows[i], c, ctx);
    }

    return status;
}

GR_INLINE int
gr_mat_swap_entrywise(gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(mat1, ctx);
    c = gr_mat_ncols(mat1, ctx);

    if (r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;

    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_swap(mat1->rows[i], mat2->rows[i], c, ctx);
    }

    return status;
}

GR_INLINE int
gr_mat_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (r != gr_mat_nrows(mat1, ctx) ||
        c != gr_mat_ncols(mat1, ctx) ||
        r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_add(res->rows[i], mat1->rows[i], mat2->rows[i], c, ctx);
    }

    return status;
}

GR_INLINE int
gr_mat_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    int status;
    slong i, r, c;

    r = gr_mat_nrows(res, ctx);
    c = gr_mat_ncols(res, ctx);

    if (r != gr_mat_nrows(mat1, ctx) ||
        c != gr_mat_ncols(mat1, ctx) ||
        r != gr_mat_nrows(mat2, ctx) ||
        c != gr_mat_ncols(mat2, ctx))
    {
        return GR_DOMAIN;
    }

    status = GR_SUCCESS;
    for (i = 0; i < r; i++)
    {
        status |= _gr_vec_sub(res->rows[i], mat1->rows[i], mat2->rows[i], c, ctx);
    }

    return status;
}

// todo use stream
// todo status
GR_INLINE int
gr_mat_print(const gr_mat_t mat, gr_ctx_t ctx)
{
    int status;
    slong r, c;
    slong i, j;
    slong sz;

    sz = ctx->sizeof_elem;
    r = gr_mat_nrows(mat, ctx);
    c = gr_mat_ncols(mat, ctx);

    status = GR_SUCCESS;
    flint_printf("[");

    for (i = 0; i < r; i++)
    {
        flint_printf("[");

        for (j = 0; j < c; j++)
        {
            status |= gr_print(GR_MAT_ENTRY(mat, i, j, sz), ctx);
            if (j < c - 1)
                flint_printf(", ");
        }

        if (i < r - 1)
            flint_printf("],\n");
        else
            flint_printf("]");
    }
    flint_printf("]\n");
    return status;
}

GR_INLINE int
gr_mat_mul_classical(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong ar, ac, br, bc, i, j, sz;
    int status;

    ar = gr_mat_nrows(A, ctx);
    ac = gr_mat_ncols(A, ctx);
    br = gr_mat_nrows(B, ctx);
    bc = gr_mat_ncols(B, ctx);

    if (ac != br || ar != gr_mat_nrows(C, ctx) || bc != gr_mat_ncols(C, ctx))
        return GR_DOMAIN;

    if (br == 0)
    {
        return gr_mat_zero(C, ctx);
    }

    status = GR_SUCCESS;

    if (A == C || B == C)
    {
        gr_mat_t T;
        status |= gr_mat_init(T, ar, bc, ctx);
        status |= gr_mat_mul_classical(T, A, B, ctx);
        status |= gr_mat_swap_entrywise(T, C, ctx);
        status |= gr_mat_clear(T, ctx);
        return status;
    }

    sz = ctx->sizeof_elem;

    if (br == 1)
    {
        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                status |= gr_mul(GR_MAT_ENTRY(C, i, j, sz),
                                 GR_MAT_ENTRY(A, i, 0, sz),
                                 GR_MAT_ENTRY(B, 0, j, sz), ctx);
            }
        }
    }
    else
    {
        gr_ptr tmp;
        TMP_INIT;

        TMP_START;
        tmp = TMP_ALLOC(sz * br * bc);

        for (i = 0; i < br; i++)
            for (j = 0; j < bc; j++)
                memcpy(GR_ENTRY(tmp, j * br + i, sz), GR_MAT_ENTRY(B, i, j, sz), sz);

        for (i = 0; i < ar; i++)
        {
            for (j = 0; j < bc; j++)
            {
                status |= _gr_vec_dot(GR_MAT_ENTRY(C, i, j, sz), NULL, 0,
                    GR_MAT_ENTRY(A, i, 0, sz), GR_ENTRY(tmp, j * br, sz), br, ctx);
            }
        }

        TMP_END;
    }

    return status;
}

/* Todo: generic polynomials */

typedef struct
{
    gr_ptr coeffs;
    slong length;
    slong alloc;
}
gr_poly_struct;

typedef gr_poly_struct gr_poly_t[1];





/*
A ring to test.
*/

#define NMOD8_CTX(ring_ctx) (((nmod_t *)((ring_ctx)->elem_ctx))[0])

typedef unsigned char nmod8_struct;
typedef nmod8_struct nmod8_t[1];

GR_INLINE int
nmod8_init(nmod8_t x, gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_clear(nmod8_t x, gr_ctx_t ctx)
{
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_swap(nmod8_t x, nmod8_t y, gr_ctx_t ctx)
{
    nmod8_t t;
    *t = *x;
    *x = *y;
    *y = *t;
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_randtest(nmod8_t res, flint_rand_t state, const void * options, gr_ctx_t ctx)
{
    res[0] = n_randtest(state) % NMOD8_CTX(ctx).n;
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_write(gr_stream_t out, const nmod8_t x, gr_ctx_t ctx)
{
    gr_stream_write_si(out, x[0]);
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_zero(nmod8_t x, gr_ctx_t ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_one(nmod8_t x, gr_ctx_t ctx)
{
    x[0] = (NMOD8_CTX(ctx).n != 0);
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_set_si(nmod8_t res, slong v, gr_ctx_t ctx)
{
    ulong t;
    nmod_t mod = NMOD8_CTX(ctx);
    t = FLINT_ABS(v);
    NMOD_RED(t, t, mod);
    if (v < 0)
        t = nmod_neg(t, mod);
    res[0] = t;
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_is_zero(int * res, const nmod8_t x, gr_ctx_t ctx)
{
    res[0] = (x[0] == 0);
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_is_one(int * res, const nmod8_t x, gr_ctx_t ctx)
{
    res[0] = (x[0] == (NMOD8_CTX(ctx).n != 1));
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_equal(int * res, const nmod8_t x, const nmod8_t y, gr_ctx_t ctx)
{
    res[0] = (x[0] == y[0]);
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_set(nmod8_t res, const nmod8_t x, gr_ctx_t ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_neg(nmod8_t res, const nmod8_t x, gr_ctx_t ctx)
{
    res[0] = nmod_neg(x[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_add(nmod8_t res, const nmod8_t x, const nmod8_t y, gr_ctx_t ctx)
{
    res[0] = nmod_add(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_add_si(nmod8_t res, const nmod8_t x, slong y, gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_add(res, x, t, ctx);
}

GR_INLINE int
nmod8_sub(nmod8_t res, const nmod8_t x, const nmod8_t y, gr_ctx_t ctx)
{
    res[0] = nmod_sub(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_mul(nmod8_t res, const nmod8_t x, const nmod8_t y, gr_ctx_t ctx)
{
    res[0] = nmod_mul(x[0], y[0], NMOD8_CTX(ctx));
    return GR_SUCCESS;
}

GR_INLINE int
nmod8_mul_si(nmod8_t res, const nmod8_t x, slong y, gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_mul(res, x, t, ctx);
}

GR_INLINE int
nmod8_inv(nmod8_t res, const nmod8_t x, gr_ctx_t ctx)
{
    ulong r, g;

    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    if (g == 1)
    {
        res[0] = r;
        return GR_SUCCESS;
    }
    else
    {
        res[0] = 0;
        return GR_DOMAIN;
    }
}

GR_INLINE int
nmod8_div(nmod8_t res, const nmod8_t x, const nmod8_t y, gr_ctx_t ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);
    nmod8_mul(res, x, t, ctx);
    return status;
}

GR_INLINE int
nmod8_div_si(nmod8_t res, const nmod8_t x, slong y, gr_ctx_t ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

GR_INLINE int
nmod8_is_invertible(int * res, const nmod8_t x, gr_ctx_t ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], NMOD8_CTX(ctx).n);
    res[0] = (g == 1);
    return GR_SUCCESS;
}

void gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n);
void gr_ctx_clear_nmod8(gr_ctx_t ctx);



/* Matrix ring to test */

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
}
matrix_ctx_t;

#define MATRIX_CTX(ring_ctx) ((matrix_ctx_t *)((ring_ctx)->elem_ctx))


GR_INLINE int
matrix_init(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_init(res, MATRIX_CTX(ctx)->n, MATRIX_CTX(ctx)->n, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_clear(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_clear(res, MATRIX_CTX(ctx)->base_ring);
}

/* TODO UNIFY */
GR_INLINE int
matrix_write(gr_stream_t out, gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_print(res, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_randtest(gr_mat_t res, flint_rand_t state, void * options, gr_ctx_t ctx)
{
    return gr_mat_randtest(res, state, options, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_equal(int * equal, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_equal(equal, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_set(res, mat, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_zero(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_zero(res, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_one(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_one(res, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_neg(res, mat, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_add(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_sub(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

GR_INLINE int
matrix_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_mul_classical(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

/* todo */
#define matrix_swap gr_not_implemented
#define matrix_is_zero gr_not_implemented
#define matrix_is_one gr_not_implemented
#define matrix_set_si gr_not_implemented

void gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n);
void gr_ctx_clear_matrix(gr_ctx_t ctx);

#ifdef __cplusplus
}
#endif

#endif
