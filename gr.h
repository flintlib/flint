#include <string.h>
#include "flint/flint.h"
#include "flint/fmpz.h"
#include "flint/fmpq.h"
#include "flint/nmod_vec.h"
#include "flint/ulong_extras.h"
#include "flint/profiler.h"

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

void gr_stream_init_file(gr_stream_t out, FILE * fp)
{
    out->fp = fp;
    out->s = NULL;
}

void gr_stream_init_str(gr_stream_t out)
{
    out->fp = NULL;
    out->s = flint_malloc(16);
    out->s[0] = '\0';
    out->len = 0;
    out->alloc = 16;
}

void gr_stream_write(gr_stream_t out, const char * s)
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

void
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

void
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

int gr_not_implemented(void) { return GR_UNABLE; }

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

const gr_method_tab_input gr_generic_methods[];

void
gr_method_tab_init(gr_method_tab_t * methods, gr_method_tab_input * tab)
{
    slong i;

    methods->methods = flint_malloc(sizeof(void *) * GR_METHOD_TAB_SIZE);

    for (i = 0; i < GR_METHOD_TAB_SIZE; i++)
        methods->methods[i] = (gr_funcptr) gr_not_implemented;

    /* Assign generic methods as fallbacks */
    for (i = 0; ; i++)
    {
        if (gr_generic_methods[i].function == NULL)
            break;

        if (gr_generic_methods[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods->methods[gr_generic_methods[i].index] = gr_generic_methods[i].function;
    }

    for (i = 0; ; i++)
    {
        if (tab[i].function == NULL)
            break;

        if (tab[i].index >= GR_METHOD_TAB_SIZE)
            abort();

        methods->methods[tab[i].index] = tab[i].function;
    }
}


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

int gr_init(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, INIT)(res, ctx->elem_ctx); }
int gr_clear(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, CLEAR)(res, ctx->elem_ctx); }
int gr_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx) { return GR_SWAP_OP(ctx, SWAP)(x, y, ctx->elem_ctx); }
int gr_randtest(gr_ptr x, flint_rand_t state, const void * options, gr_ctx_t ctx) { return GR_RANDTEST(ctx, RANDTEST)(x, state, options, ctx->elem_ctx); }
int gr_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx) { return GR_STREAM_IN(ctx, WRITE)(out, x, ctx->elem_ctx); }
int gr_zero(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ZERO)(res, ctx->elem_ctx); }
int gr_one(gr_ptr res, gr_ctx_t ctx) { return GR_CONSTANT_OP(ctx, ONE)(res, ctx->elem_ctx); }
int gr_set_si(gr_ptr res, slong x, gr_ctx_t ctx) { return GR_UNARY_OP_SI(ctx, SET_SI)(res, x, ctx->elem_ctx); }
int gr_is_zero(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ZERO)(res, x, ctx->elem_ctx); }
int gr_is_one(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_ONE)(res, x, ctx->elem_ctx); }
int gr_equal(int * res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_PREDICATE(ctx, EQUAL)(res, x, y, ctx->elem_ctx); }
int gr_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, SET)(res, x, ctx->elem_ctx); }
int gr_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, NEG)(res, x, ctx->elem_ctx); }
int gr_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, ADD)(res, x, y, ctx->elem_ctx); }
int gr_add_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, ADD_SI)(res, x, y, ctx->elem_ctx); }
int gr_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, SUB)(res, x, y, ctx->elem_ctx); }
int gr_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, MUL)(res, x, y, ctx->elem_ctx); }
int gr_mul_si(gr_ptr res, gr_srcptr x, slong y, gr_ctx_t ctx) { return GR_BINARY_OP_SI(ctx, MUL_SI)(res, x, y, ctx->elem_ctx); }
int gr_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx) { return GR_BINARY_OP(ctx, DIV)(res, x, y, ctx->elem_ctx); }
int gr_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_OP(ctx, INV)(res, x, ctx->elem_ctx); }
int gr_is_invertible(int * res, gr_srcptr x, gr_ctx_t ctx) { return GR_UNARY_PREDICATE(ctx, IS_INVERTIBLE)(res, x, ctx->elem_ctx); }

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
        init(x1, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_INIT2(x1, x2, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(2 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        init(x1, (ctx)->elem_ctx); \
        init(x2, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_INIT3(x1, x2, x3, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(3 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        init(x1, (ctx)->elem_ctx); \
        init(x2, (ctx)->elem_ctx); \
        init(x3, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_INIT4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_constant_op init = GR_CONSTANT_OP(ctx, INIT); \
        ssize_t _gr_elem_size = (ctx)->sizeof_elem; \
        x1 = (gr_ptr) GR_TMP_ALLOC(4 * _gr_elem_size); \
        x2 = (gr_ptr) ((char *) x1 + _gr_elem_size); \
        x3 = (gr_ptr) ((char *) x2 + _gr_elem_size); \
        x4 = (gr_ptr) ((char *) x3 + _gr_elem_size); \
        init(x1, (ctx)->elem_ctx); \
        init(x2, (ctx)->elem_ctx); \
        init(x3, (ctx)->elem_ctx); \
        init(x4, (ctx)->elem_ctx); \
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
        init(x1, (ctx)->elem_ctx); \
        init(x2, (ctx)->elem_ctx); \
        init(x3, (ctx)->elem_ctx); \
        init(x4, (ctx)->elem_ctx); \
        init(x5, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_CLEAR1(x1, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_CLEAR2(x1, x2, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)->elem_ctx); \
        clear(x2, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_CLEAR3(x1, x2, x3, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)->elem_ctx); \
        clear(x2, (ctx)->elem_ctx); \
        clear(x3, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_CLEAR4(x1, x2, x3, x4, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)->elem_ctx); \
        clear(x2, (ctx)->elem_ctx); \
        clear(x3, (ctx)->elem_ctx); \
        clear(x4, (ctx)->elem_ctx); \
    } while (0)

#define GR_TMP_CLEAR5(x1, x2, x3, x4, x5, ctx) \
    do { \
        gr_method_constant_op clear = GR_CONSTANT_OP(ctx, CLEAR); \
        clear(x1, (ctx)->elem_ctx); \
        clear(x2, (ctx)->elem_ctx); \
        clear(x3, (ctx)->elem_ctx); \
        clear(x4, (ctx)->elem_ctx); \
        clear(x5, (ctx)->elem_ctx); \
    } while (0)


/* Generic vector functions */

#define GR_ENTRY(vec, i, size) ((void *) (((char *) (vec)) + ((i) * (size))))

int
gr_generic_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_init(GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
gr_generic_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_clear(GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
gr_generic_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_swap(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);

    return status;
}

int
gr_generic_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_zero(GR_ENTRY(vec, i, sz), ctx);

    return status;
}

int
gr_generic_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_set(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

int
gr_generic_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_neg(GR_ENTRY(res, i, sz), GR_ENTRY(src, i, sz), ctx);

    return status;
}

int
gr_generic_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_add(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz), ctx);

    return status;
}

int
gr_generic_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx)
{
    int status;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    for (i = 0; i < len; i++)
        status |= gr_sub(GR_ENTRY(res, i, sz), GR_ENTRY(src1, i, sz), GR_ENTRY(src2, i, sz), ctx);

    return status;
}

int
gr_generic_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= gr_add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_mul(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= gr_sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= gr_add(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    for (i = 0; i < len; i++)
    {
        status |= gr_mul_si(t, GR_ENTRY(vec2, i, sz), c, ctx);
        status |= gr_sub(GR_ENTRY(vec1, i, sz), GR_ENTRY(vec1, i, sz), t, ctx);
    }

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_equal(int * res, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    int status, equal, this_equal;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    equal = 1;
    for (i = 0; i < len; i++)
    {
        status |= gr_equal(&this_equal, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
        equal = equal && this_equal;
    }

    res[0] = equal;
    return status;
}

int
gr_generic_vec_is_zero(int * res, gr_srcptr vec, slong len, gr_ctx_t ctx)
{
    int status, equal, this_equal;
    slong i, sz;

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    equal = 1;
    for (i = 0; i < len; i++)
    {
        status |= gr_is_zero(&this_equal, GR_ENTRY(vec, i, sz), ctx);
        equal = equal && this_equal;
    }

    res[0] = equal;
    return status;
}

int
gr_generic_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    if (initial == NULL)
    {
        status |= gr_mul(res, vec1, vec2, ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= gr_mul(t, vec1, vec2, ctx);
        status |= gr_add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= gr_mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, i, sz), ctx);
        status |= gr_add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

int
gr_generic_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx)
{
    GR_TMP_START;
    int status;
    slong i, sz;
    gr_ptr t;

    if (len <= 0)
    {
        if (initial == NULL)
            return gr_zero(res, ctx);
        else
            return gr_set(res, initial, ctx);
    }

    sz = ctx->sizeof_elem;
    status = GR_SUCCESS;

    GR_TMP_INIT1(t, ctx);

    if (initial == NULL)
    {
        status |= gr_mul(res, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
    }
    else
    {
        if (subtract)
            status |= gr_neg(res, initial, ctx);
        else
            status |= gr_set(res, initial, ctx);

        status |= gr_mul(t, vec1, GR_ENTRY(vec2, len - 1, sz), ctx);
        status |= gr_add(res, res, t, ctx);
    }

    for (i = 1; i < len; i++)
    {
        status |= gr_mul(t, GR_ENTRY(vec1, i, sz), GR_ENTRY(vec2, len - 1 - i, sz), ctx);
        status |= gr_add(res, res, t, ctx);
    }

    if (subtract)
        status |= gr_neg(res, res, ctx);

    GR_TMP_CLEAR1(t, ctx);
    GR_TMP_END;
    return status;
}

/* Generic method implementations */
const gr_method_tab_input gr_generic_methods[] =
{
    {GR_METHOD_VEC_INIT,                (gr_funcptr) gr_generic_vec_init},
    {GR_METHOD_VEC_CLEAR,               (gr_funcptr) gr_generic_vec_clear},
    {GR_METHOD_VEC_SWAP,                (gr_funcptr) gr_generic_vec_swap},
    {GR_METHOD_VEC_ZERO,                (gr_funcptr) gr_generic_vec_zero},
    {GR_METHOD_VEC_SET,                 (gr_funcptr) gr_generic_vec_set},
    {GR_METHOD_VEC_NEG,                 (gr_funcptr) gr_generic_vec_neg},
    {GR_METHOD_VEC_ADD,                 (gr_funcptr) gr_generic_vec_add},
    {GR_METHOD_VEC_SUB,                 (gr_funcptr) gr_generic_vec_sub},
    {GR_METHOD_VEC_SCALAR_ADDMUL,       (gr_funcptr) gr_generic_vec_scalar_addmul},
    {GR_METHOD_VEC_SCALAR_SUBMUL,       (gr_funcptr) gr_generic_vec_scalar_submul},
    {GR_METHOD_VEC_SCALAR_ADDMUL_SI,    (gr_funcptr) gr_generic_vec_scalar_addmul_si},
    {GR_METHOD_VEC_SCALAR_SUBMUL_SI,    (gr_funcptr) gr_generic_vec_scalar_submul_si},
    {GR_METHOD_VEC_EQUAL,               (gr_funcptr) gr_generic_vec_equal},
    {GR_METHOD_VEC_IS_ZERO,             (gr_funcptr) gr_generic_vec_is_zero},
    {GR_METHOD_VEC_DOT,                 (gr_funcptr) gr_generic_vec_dot},
    {GR_METHOD_VEC_DOT_REV,             (gr_funcptr) gr_generic_vec_dot_rev},
    {0,                                 (gr_funcptr) NULL},
};

int _gr_vec_init(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_INIT)(vec, len, ctx); }
int _gr_vec_clear(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_CLEAR)(vec, len, ctx); }
int _gr_vec_swap(gr_ptr vec1, gr_ptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_SWAP_OP(ctx, VEC_SWAP)(vec1, vec2, len, ctx); }
int _gr_vec_zero(gr_ptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_CONSTANT_OP(ctx, VEC_ZERO)(vec, len, ctx); }
int _gr_vec_set(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_SET)(res, src, len, ctx); }
int _gr_vec_neg(gr_ptr res, gr_srcptr src, slong len, gr_ctx_t ctx) { return GR_VEC_OP(ctx, VEC_NEG)(res, src, len, ctx); }
int _gr_vec_add(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_ADD)(res, src1, src2, len, ctx); }
int _gr_vec_sub(gr_ptr res, gr_srcptr src1, gr_srcptr src2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_OP(ctx, VEC_SUB)(res, src1, src2, len, ctx); }
int _gr_vec_scalar_addmul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_ADDMUL)(vec1, vec2, len, c, ctx); }
int _gr_vec_scalar_submul(gr_ptr vec1, gr_srcptr vec2, slong len, gr_srcptr c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP(ctx, VEC_SCALAR_SUBMUL)(vec1, vec2, len, c, ctx); }
int _gr_vec_scalar_addmul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_ADDMUL_SI)(vec1, vec2, len, c, ctx); }
int _gr_vec_scalar_submul_si(gr_ptr vec1, gr_srcptr vec2, slong len, slong c, gr_ctx_t ctx) { return GR_VEC_SCALAR_OP_SI(ctx, VEC_SCALAR_SUBMUL_SI)(vec1, vec2, len, c, ctx); }
int _gr_vec_equal(int * res, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_VEC_PREDICATE(ctx, VEC_EQUAL)(res, vec1, vec2, len, ctx); }
int _gr_vec_is_zero(int * res, gr_srcptr vec, slong len, gr_ctx_t ctx) { return GR_VEC_PREDICATE(ctx, VEC_IS_ZERO)(res, vec, len, ctx); }
int _gr_vec_dot(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT)(res, initial, subtract, vec1, vec2, len, ctx); }
int _gr_vec_dot_rev(gr_ptr res, gr_srcptr initial, int subtract, gr_srcptr vec1, gr_srcptr vec2, slong len, gr_ctx_t ctx) { return GR_VEC_DOT_OP(ctx, VEC_DOT_REV)(res, initial, subtract, vec1, vec2, len, ctx); }

/* todo: could allow overloading this as well */
int
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


/* Some generic functions. */

/* Assumes exp >= 2; res and tmp not not aliased with x. */
static int
_gr_generic_pow_ui_binexp(gr_ptr res, gr_ptr tmp, gr_srcptr x, ulong exp, gr_ctx_t ctx)
{
    gr_ptr R, S, T;
    int status;
    int zeros;
    ulong bit;

    status = GR_SUCCESS;

    /* Determine parity due to swaps */
    zeros = 0;
    bit = exp;
    while (bit > 1)
    {
        zeros += !(bit & 1);
        bit >>= 1;
    }

    if (zeros % 2)
    {
        R = res;
        S = tmp;
    }
    else
    {
        R = tmp;
        S = res;
    }

    bit = UWORD(1) << (FLINT_BIT_COUNT(exp) - 2);

    status |= gr_mul(R, x, x, ctx);

    if (bit & exp)
    {
        status |= gr_mul(S, R, x, ctx);
        T = R;
        R = S;
        S = T;
    }

    while (bit >>= 1)
    {
        status |= gr_mul(S, R, R, ctx);

        if (bit & exp)
        {
            status |= gr_mul(R, S, x, ctx);
        }
        else
        {
            T = R;
            R = S;
            S = T;
        }
    }

    return status;
}

int
gr_generic_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    int status;

    if (e == 0)
    {
        return gr_one(res, ctx);
    }
    else if (e == 1)
    {
        return gr_set(res, x, ctx);
    }
    else if (e == 2)
    {
        return gr_mul(res, x, x, ctx);
    }
    else if (e == 4)
    {
        status = gr_mul(res, x, x, ctx);
        status |= gr_mul(res, res, res, ctx);
        return status;
    }
    else
    {
        gr_ptr t, u;
        GR_TMP_START;

        if (res == x)
        {
            GR_TMP_INIT2(t, u, ctx);
            gr_set(u, x, ctx);
            status = _gr_generic_pow_ui_binexp(res, t, u, e, ctx);
            GR_TMP_CLEAR2(t, u, ctx);
        }
        else
        {
            GR_TMP_INIT1(t, ctx);
            status = _gr_generic_pow_ui_binexp(res, t, x, e, ctx);
            GR_TMP_CLEAR1(t, ctx);
        }

        GR_TMP_END;
        return status;
    }
}

int
gr_generic_pow_si(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    if (e >= 0)
    {
        return gr_generic_pow_ui(res, x, e, ctx);
    }
    else
    {
        int status;

        /* todo: some heuristic for when we want to invert before/after powering */
        status = gr_inv(res, x, ctx);
        if (status == GR_SUCCESS)
            status = gr_generic_pow_ui(res, x, -e, ctx);

        return status;
    }
}

/* tbd */
int
gr_pow_ui(gr_ptr res, gr_srcptr x, ulong e, gr_ctx_t ctx)
{
    return gr_generic_pow_ui(res, x, e, ctx);
}

int
gr_print(gr_srcptr x, gr_ctx_t ctx)
{
    gr_stream_t out;
    gr_stream_init_file(out, stdout);
    gr_write(out, x, ctx);
    return GR_SUCCESS;
}

int
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

int
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

int
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

int
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

int
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

int
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

int
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

int
gr_mat_one(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_set_si(res, 1, ctx);
}

int
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

int
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

int
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

int
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

int
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
int
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

int
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

typedef unsigned char nmod8_struct;
typedef nmod8_struct nmod8_t[1];

int
nmod8_init(nmod8_t x, nmod_t * ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

int
nmod8_clear(nmod8_t x, nmod_t * ctx)
{
    return GR_SUCCESS;
}

int
nmod8_swap(nmod8_t x, nmod8_t y, nmod_t * ctx)
{
    nmod8_t t;
    *t = *x;
    *x = *y;
    *y = *t;
    return GR_SUCCESS;
}

int
nmod8_randtest(nmod8_t res, flint_rand_t state, const void * options, nmod_t * ctx)
{
    res[0] = n_randtest(state) % ctx->n;
    return GR_SUCCESS;
}

int
nmod8_write(gr_stream_t out, const nmod8_t x, nmod_t * ctx)
{
    gr_stream_write_si(out, x[0]);
    return GR_SUCCESS;
}

int
nmod8_zero(nmod8_t x, nmod_t * ctx)
{
    x[0] = 0;
    return GR_SUCCESS;
}

int
nmod8_one(nmod8_t x, nmod_t * ctx)
{
    x[0] = (ctx->n != 0);
    return GR_SUCCESS;
}

int
nmod8_set_si(nmod8_t res, slong v, nmod_t * ctx)
{
    ulong t;
    t = FLINT_ABS(v);
    NMOD_RED(t, t, *ctx);
    if (v < 0)
        t = nmod_neg(t, *ctx);
    res[0] = t;
    return GR_SUCCESS;
}

int
nmod8_is_zero(int * res, const nmod8_t x, nmod_t * ctx)
{
    res[0] = (x[0] == 0);
    return GR_SUCCESS;
}

int
nmod8_is_one(int * res, const nmod8_t x, nmod_t * ctx)
{
    res[0] = (x[0] == (ctx->n != 1));
    return GR_SUCCESS;
}

int
nmod8_equal(int * res, const nmod8_t x, const nmod8_t y, nmod_t * ctx)
{
    res[0] = (x[0] == y[0]);
    return GR_SUCCESS;
}

int
nmod8_set(nmod8_t res, const nmod8_t x, nmod_t * ctx)
{
    res[0] = x[0];
    return GR_SUCCESS;
}

int
nmod8_neg(nmod8_t res, const nmod8_t x, nmod_t * ctx)
{
    res[0] = nmod_neg(x[0], *ctx);
    return GR_SUCCESS;
}

int
nmod8_add(nmod8_t res, const nmod8_t x, const nmod8_t y, nmod_t * ctx)
{
    res[0] = nmod_add(x[0], y[0], *ctx);
    return GR_SUCCESS;
}

int
nmod8_add_si(nmod8_t res, const nmod8_t x, slong y, nmod_t * ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_add(res, x, t, ctx);
}

int
nmod8_sub(nmod8_t res, const nmod8_t x, const nmod8_t y, nmod_t * ctx)
{
    res[0] = nmod_sub(x[0], y[0], *ctx);
    return GR_SUCCESS;
}

int
nmod8_mul(nmod8_t res, const nmod8_t x, const nmod8_t y, nmod_t * ctx)
{
    res[0] = nmod_mul(x[0], y[0], *ctx);
    return GR_SUCCESS;
}

int
nmod8_mul_si(nmod8_t res, const nmod8_t x, slong y, nmod_t * ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_mul(res, x, t, ctx);
}

int
nmod8_inv(nmod8_t res, const nmod8_t x, nmod_t * ctx)
{
    ulong r, g;

    g = n_gcdinv(&r, x[0], ctx->n);
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

int
nmod8_div(nmod8_t res, const nmod8_t x, const nmod8_t y, nmod_t * ctx)
{
    nmod8_t t;
    int status;

    status = nmod8_inv(t, y, ctx);
    nmod8_mul(res, x, t, ctx);
    return status;
}

int
nmod8_div_si(nmod8_t res, const nmod8_t x, slong y, nmod_t * ctx)
{
    nmod8_t t;
    nmod8_set_si(t, y, ctx);
    return nmod8_div(res, x, t, ctx);
}

int
nmod8_is_invertible(int * res, const nmod8_t x, nmod_t * ctx)
{
    ulong r, g;
    g = n_gcdinv(&r, x[0], ctx->n);
    res[0] = (g == 1);
    return GR_SUCCESS;
}



/* todo: thread safe */
int _nmod8_methods2_initialized = 0;
gr_method_tab_t _nmod8_methods2;

gr_method_tab_input nmod8_methods2[] =
{
    {GR_METHOD_INIT,            (gr_funcptr) nmod8_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) nmod8_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) nmod8_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) nmod8_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) nmod8_write},
    {GR_METHOD_ZERO,            (gr_funcptr) nmod8_zero},
    {GR_METHOD_ONE,             (gr_funcptr) nmod8_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) nmod8_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) nmod8_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) nmod8_equal},
    {GR_METHOD_SET,             (gr_funcptr) nmod8_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) nmod8_set_si},
    {GR_METHOD_NEG,             (gr_funcptr) nmod8_neg},
    {GR_METHOD_ADD,             (gr_funcptr) nmod8_add},
    {GR_METHOD_ADD_SI,          (gr_funcptr) nmod8_add_si},
    {GR_METHOD_SUB,             (gr_funcptr) nmod8_sub},
    {GR_METHOD_MUL,             (gr_funcptr) nmod8_mul},
    {GR_METHOD_MUL_SI,          (gr_funcptr) nmod8_mul_si},
    {GR_METHOD_DIV,             (gr_funcptr) nmod8_div},
    {GR_METHOD_DIV_SI,          (gr_funcptr) nmod8_div_si},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) nmod8_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) nmod8_inv},
    {0,                         (gr_funcptr) NULL},
};


void
gr_ctx_init_nmod8(gr_ctx_t ctx, unsigned char n)
{
    ctx->flags = GR_COMMUTATIVE_RING;
    ctx->sizeof_elem = sizeof(nmod8_struct);
    ctx->elem_ctx = flint_malloc(sizeof(nmod_t));  /* This could be something more interesting */
    nmod_init(ctx->elem_ctx, n);

    if (!_nmod8_methods2_initialized)
    {
        gr_method_tab_init(&_nmod8_methods2, nmod8_methods2);
        _nmod8_methods2_initialized = 1;
    }

    ctx->methods2 = &_nmod8_methods2;

    ctx->debug_string = "nmod8 ring";
}

void
gr_ctx_clear_nmod8(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
}

/* Matrix ring to test */

typedef struct
{
    gr_ctx_struct * base_ring;
    slong n;
}
matrix_ctx_t;

int
matrix_init(gr_mat_t res, matrix_ctx_t * ctx)
{
    return gr_mat_init(res, ctx->n, ctx->n, ctx->base_ring);
}

int
matrix_clear(gr_mat_t res, matrix_ctx_t * ctx)
{
    return gr_mat_clear(res, ctx->base_ring);
}

/* TODO UNIFY */
int
matrix_write(gr_stream_t out, gr_mat_t res, matrix_ctx_t * ctx)
{
    return gr_mat_print(res, ctx->base_ring);
}

int
matrix_randtest(gr_mat_t res, flint_rand_t state, void * options, matrix_ctx_t * ctx)
{
    return gr_mat_randtest(res, state, options, ctx->base_ring);
}

int
matrix_equal(int * equal, const gr_mat_t mat1, const gr_mat_t mat2, matrix_ctx_t * ctx)
{
    return gr_mat_equal(equal, mat1, mat2, ctx->base_ring);
}

int
matrix_set(gr_mat_t res, const gr_mat_t mat, matrix_ctx_t * ctx)
{
    return gr_mat_set(res, mat, ctx->base_ring);
}

int
matrix_zero(gr_mat_t res, matrix_ctx_t * ctx)
{
    return gr_mat_zero(res, ctx->base_ring);
}

int
matrix_one(gr_mat_t res, matrix_ctx_t * ctx)
{
    return gr_mat_one(res, ctx->base_ring);
}

int
matrix_neg(gr_mat_t res, const gr_mat_t mat, matrix_ctx_t * ctx)
{
    return gr_mat_neg(res, mat, ctx->base_ring);
}

int
matrix_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, matrix_ctx_t * ctx)
{
    return gr_mat_add(res, mat1, mat2, ctx->base_ring);
}

int
matrix_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, matrix_ctx_t * ctx)
{
    return gr_mat_sub(res, mat1, mat2, ctx->base_ring);
}

int
matrix_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, matrix_ctx_t * ctx)
{
    return gr_mat_mul_classical(res, mat1, mat2, ctx->base_ring);
}

/* todo */
#define matrix_swap gr_not_implemented
#define matrix_is_zero gr_not_implemented
#define matrix_is_one gr_not_implemented
#define matrix_set_si gr_not_implemented

/* todo: thread safe */
int _matrix_methods2_initialized = 0;
gr_method_tab_t _matrix_methods2;

gr_method_tab_input matrix_methods2[] =
{
    {GR_METHOD_INIT,        (gr_funcptr) matrix_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) matrix_clear},
    {GR_METHOD_RANDTEST,    (gr_funcptr) matrix_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) matrix_write},
    {GR_METHOD_ZERO,        (gr_funcptr) matrix_zero},
    {GR_METHOD_ONE,         (gr_funcptr) matrix_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) matrix_equal},
    {GR_METHOD_SET,         (gr_funcptr) matrix_set},
    {GR_METHOD_NEG,         (gr_funcptr) matrix_neg},
    {GR_METHOD_ADD,         (gr_funcptr) matrix_add},
    {GR_METHOD_SUB,         (gr_funcptr) matrix_sub},
    {GR_METHOD_MUL,         (gr_funcptr) matrix_mul},
    {0,                     (gr_funcptr) NULL},
};

/* rename: gr_ctx_init_gr_mat */
void
gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    ctx->flags = 0;
    ctx->sizeof_elem = sizeof(gr_mat_struct);
    ctx->elem_ctx = flint_malloc(sizeof(matrix_ctx_t));
    ((matrix_ctx_t *) ctx->elem_ctx)->base_ring = (gr_ctx_struct *) base_ring;
    ((matrix_ctx_t *) ctx->elem_ctx)->n = n;

    if (!_matrix_methods2_initialized)
    {
        gr_method_tab_init(&_matrix_methods2, matrix_methods2);
        _matrix_methods2_initialized = 1;
    }

    ctx->methods2 = &_matrix_methods2;

    ctx->debug_string = "matrix ring";
}

void
gr_ctx_clear_matrix(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
}

