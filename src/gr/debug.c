#include <stdio.h>
#include "gr.h"
#include "profiler.h"

typedef struct
{
    gr_ptr elem;
    gr_ctx_struct * elem_ctx;
    int cleared;
}
gr_debug_wrap_struct;

typedef struct
{
    gr_ctx_struct * cctx;
    int flags;
    float unable_probability;
    flint_rand_struct * rand_state;
}
_gr_debug_ctx_t;

#define GR_DEBUG_CTX(ctx) (((_gr_debug_ctx_t *)(ctx)))
#define GR_DEBUG_ELEM_CTX(ctx) (GR_DEBUG_CTX(ctx)->cctx)
#define GR_DEBUG_CTX_FLAGS(ctx) (GR_DEBUG_CTX(ctx)->flags)
#define GR_DEBUG_CTX_RAND_STATE(ctx) (GR_DEBUG_CTX(ctx)->rand_state)
#define GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) (GR_DEBUG_CTX(ctx)->unable_probability)

#define GR_DEBUG_ELEM(xx) ((GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP) ? (((gr_debug_wrap_struct *) (xx))->elem) : (xx))

#define CHECK_TRUTH(t) \
    if ((GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_CHECK_ALWAYS_ABLE) && ((t) == T_UNKNOWN)) \
        flint_throw(FLINT_ERROR, "unexpected T_UNKNOWN");

#define CHECK_STATUS \
    if ((GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_CHECK_ALWAYS_ABLE) && (status & GR_UNABLE)) \
        flint_throw(FLINT_ERROR, "unexpected GR_UNABLE");

#define CHECK_WRAPPER(x) \
    if ((GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP) && (((const gr_debug_wrap_struct *) x)->elem_ctx != GR_DEBUG_ELEM_CTX(ctx))) \
    { \
        flint_printf("\n%{gr_ctx} at %p\n", ctx, ctx); \
        flint_printf("%{gr_ctx} at %p\n", GR_DEBUG_ELEM_CTX(ctx), GR_DEBUG_ELEM_CTX(ctx)); \
        flint_printf("%{gr_ctx} at %p\n", ((const gr_debug_wrap_struct *) x)->elem_ctx, ((const gr_debug_wrap_struct *) x)->elem_ctx); \
        flint_throw(FLINT_ERROR, "ctx of wrapped element does not match ctx of debug context"); \
    }

#define DEBUG_PRINT_MSG(msg) \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE) \
        flint_printf("%s\n", msg);

#define DEBUG_PRINT_ELEM(name, x) \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE) \
    { \
        flint_printf("%s = %{gr}\n", name, GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
    }

#define DEBUG_PRINT_STATUS \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE) \
    { \
        if (status == (GR_DOMAIN | GR_UNABLE)) flint_printf("status = GR_DOMAIN | GR_UNABLE\n"); \
        else if (status == GR_DOMAIN) flint_printf("status = GR_DOMAIN\n"); \
        else if (status == GR_UNABLE) flint_printf("status = GR_UNABLE\n"); \
        else if (status == GR_SUCCESS) flint_printf("status = GR_SUCCESS\n"); \
        else flint_printf("status = ???"); \
    }

#define DEBUG_PRINT_TRUTH(name, t) \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE) \
    { \
        flint_printf("%s = %{truth}\n", name, t); \
    }

#define RANDOMLY_UNABLE \
    if (GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) != 0.0) \
    { \
        if (n_randlimb(GR_DEBUG_CTX_RAND_STATE(ctx)) < GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) * (double) UWORD_MAX) \
        { \
            DEBUG_PRINT_MSG("randomly returning GR_UNABLE"); \
            return GR_UNABLE; \
        } \
    }

#define RANDOMLY_UNKNOWN \
    if (GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) != 0.0) \
    { \
        if (n_randlimb(GR_DEBUG_CTX_RAND_STATE(ctx)) < GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) * (double) UWORD_MAX) \
        { \
            DEBUG_PRINT_MSG("randomly returning T_UNKNOWN"); \
            return T_UNKNOWN; \
        } \
    }




static void
_gr_debug_ctx_clear(gr_ctx_t ctx)
{
    flint_rand_clear(GR_DEBUG_CTX_RAND_STATE(ctx));
    flint_free(GR_DEBUG_CTX_RAND_STATE(ctx));
}

static int
_gr_debug_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    status |= gr_stream_write(out, "Debugging context (flags = ");
    status |= gr_stream_write_si(out, GR_DEBUG_CTX_FLAGS(ctx));
    status |= gr_stream_write(out, ", unable_probability = ");
    char tmp[1024];
    sprintf(tmp, "%.3f", GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx));
    status |= gr_stream_write(out, tmp);
    status |= gr_stream_write(out, ") for ");
    status |= gr_ctx_write(out, GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

static truth_t
_gr_debug_ctx_is_threadsafe(gr_ctx_t ctx)
{
    /* currently threadsafe because we access the random state */
    if (GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) != 0.0)
        return T_FALSE;

    return gr_ctx_is_threadsafe(GR_DEBUG_ELEM_CTX(ctx));
}

static void
_gr_debug_init(gr_ptr x, gr_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_DEBUG_ELEM_CTX(ctx);

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE)
        flint_printf("gr_init: %p, %{gr_ctx}\n", x, cctx);

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP)
    {
        gr_debug_wrap_struct * xx = x;

        xx->elem = flint_malloc(cctx->sizeof_elem);
        gr_init(xx->elem, cctx);
        xx->elem_ctx = cctx;
        xx->cleared = 0;
    }
    else
    {
        gr_init(x, GR_DEBUG_ELEM_CTX(ctx));
    }
}

static void
_gr_debug_clear(gr_ptr x, gr_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_DEBUG_ELEM_CTX(ctx);

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_VERBOSE)
        flint_printf("gr_clear: %p, %{gr_ctx}\n", x, cctx);

    CHECK_WRAPPER(x)

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP)
    {
        gr_debug_wrap_struct * xx = x;

        if (xx->cleared)
            flint_throw(FLINT_ERROR, "gr_clear: element cleared twice");

        gr_clear(xx->elem, cctx);
        flint_free(xx->elem);
        xx->cleared = 1;
    }
    else
    {
        gr_clear(x, cctx);
    }
}

static void
_gr_debug_swap(gr_ptr x, gr_ptr y, gr_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_DEBUG_ELEM_CTX(ctx);

    CHECK_WRAPPER(x)
    CHECK_WRAPPER(y)

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP)
    {
        gr_debug_wrap_struct * xx = x;
        gr_debug_wrap_struct * yy = y;
        FLINT_SWAP(gr_debug_wrap_struct, *xx, *yy);
    }
    else
    {
        gr_swap(x, y, cctx);
    }
}

static void
_gr_debug_set_shallow(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    gr_ctx_struct * cctx = GR_DEBUG_ELEM_CTX(ctx);

    CHECK_WRAPPER(x)

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP)
    {
        gr_debug_wrap_struct * rres = res;
        const gr_debug_wrap_struct * xx = x;
        *rres = *xx;
    }
    else
    {
        gr_set_shallow(res, x, cctx);
    }

    CHECK_WRAPPER(res)
}

static int
_gr_debug_randtest(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;
    CHECK_WRAPPER(res)
    status = gr_randtest(GR_DEBUG_ELEM(res), state, GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

static int
_gr_debug_randtest_small(gr_ptr res, flint_rand_t state, gr_ctx_t ctx)
{
    int status;
    CHECK_WRAPPER(res)
    status = gr_randtest_small(GR_DEBUG_ELEM(res), state, GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

static int
_gr_debug_write(gr_stream_t out, gr_srcptr x, gr_ctx_t ctx)
{
    int status;
    CHECK_WRAPPER(x)
    status = gr_write(out, GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

#define DEBUG_UNARY_PREDICATE_WITHOUT_UNABLE(name) \
    truth_t res; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(x)"); \
    DEBUG_PRINT_ELEM("x", x); \
    CHECK_WRAPPER(x) \
    res = name(GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
    DEBUG_PRINT_TRUTH("res", res) \
    CHECK_TRUTH(res) \
    DEBUG_PRINT_MSG(""); \
    return res;

#define DEBUG_UNARY_PREDICATE(name) \
    truth_t res; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(x)"); \
    DEBUG_PRINT_ELEM("x", x); \
    CHECK_WRAPPER(x) \
    RANDOMLY_UNKNOWN \
    res = name(GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
    DEBUG_PRINT_TRUTH("res", res) \
    CHECK_TRUTH(res) \
    DEBUG_PRINT_MSG(""); \
    return res;

#define DEBUG_BINARY_PREDICATE(name) \
    truth_t res; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(x, y)"); \
    DEBUG_PRINT_ELEM("x", x); \
    DEBUG_PRINT_ELEM("y", y); \
    CHECK_WRAPPER(x) \
    CHECK_WRAPPER(y) \
    RANDOMLY_UNKNOWN \
    res = name(GR_DEBUG_ELEM(x), GR_DEBUG_ELEM(y), GR_DEBUG_ELEM_CTX(ctx)); \
    DEBUG_PRINT_TRUTH("res", res) \
    CHECK_TRUTH(res) \
    DEBUG_PRINT_MSG(""); \
    return res;

#define DEBUG_UNARY_OP_WITHOUT_UNABLE(name) \
    int status; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(res, x)"); \
    if (res == x) DEBUG_PRINT_MSG("res is aliased with x"); \
    DEBUG_PRINT_ELEM("x   (before)", x); \
    if (res != x) DEBUG_PRINT_ELEM("res (before)", res); \
    CHECK_WRAPPER(res) \
    CHECK_WRAPPER(x) \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_TIMING) \
    { \
        TIMEIT_ONCE_START; \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
        TIMEIT_ONCE_STOP; \
    } \
    else \
    { \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
    } \
    DEBUG_PRINT_ELEM("x   (after) ", x); \
    DEBUG_PRINT_ELEM("res (after) ", res); \
    DEBUG_PRINT_STATUS \
    DEBUG_PRINT_MSG(""); \
    return status;

#define DEBUG_UNARY_OP(name) \
    int status; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(res, x)"); \
    DEBUG_PRINT_ELEM("x   (before)", x); \
    if (res == x) DEBUG_PRINT_MSG("res is aliased with x"); \
    if (res != x) DEBUG_PRINT_ELEM("res (before)", res); \
    CHECK_WRAPPER(res) \
    CHECK_WRAPPER(x) \
    RANDOMLY_UNABLE \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_TIMING) \
    { \
        TIMEIT_ONCE_START; \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
        TIMEIT_ONCE_STOP; \
    } \
    else \
    { \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM_CTX(ctx)); \
    } \
    DEBUG_PRINT_ELEM("res (after) ", res); \
    DEBUG_PRINT_STATUS \
    DEBUG_PRINT_MSG(""); \
    return status;

#define DEBUG_BINARY_OP(name) \
    int status; \
    DEBUG_PRINT_MSG(""); \
    DEBUG_PRINT_MSG(#name "(res, x, y)"); \
    if (x == y)   DEBUG_PRINT_MSG("y is aliased with x"); \
    if (x != y) DEBUG_PRINT_ELEM("y   (before)", y); \
    if (res == x) DEBUG_PRINT_MSG("res is aliased with x"); \
    if (res == y) DEBUG_PRINT_MSG("res is aliased with y"); \
    if (res != x && res != y) DEBUG_PRINT_ELEM("res (before)", res); \
    DEBUG_PRINT_ELEM("x   (before)", x); \
    CHECK_WRAPPER(res) \
    CHECK_WRAPPER(x) \
    CHECK_WRAPPER(y) \
    RANDOMLY_UNABLE \
    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_TIMING) \
    { \
        TIMEIT_ONCE_START; \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM(y), GR_DEBUG_ELEM_CTX(ctx)); \
        TIMEIT_ONCE_STOP; \
    } \
    else \
    { \
        status = name(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM(x), GR_DEBUG_ELEM(y), GR_DEBUG_ELEM_CTX(ctx)); \
    } \
    DEBUG_PRINT_ELEM("res (after) ", res); \
    DEBUG_PRINT_STATUS \
    DEBUG_PRINT_MSG(""); \
    return status;

static int
_gr_debug_zero(gr_ptr res, gr_ctx_t ctx)
{
    int status;
    CHECK_WRAPPER(res)
    /* RANDOMLY_UNABLE */
    status = gr_zero(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

static int
_gr_debug_one(gr_ptr res, gr_ctx_t ctx)
{
    int status;
    CHECK_WRAPPER(res)
    /* RANDOMLY_UNABLE */
    status = gr_one(GR_DEBUG_ELEM(res), GR_DEBUG_ELEM_CTX(ctx));
    return status;
}

static truth_t
_gr_debug_is_zero(gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_PREDICATE_WITHOUT_UNABLE(gr_is_zero)
}

static truth_t
_gr_debug_is_one(gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_PREDICATE(gr_is_one)
}

static truth_t
_gr_debug_equal(gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_PREDICATE(gr_equal)
}

static int
_gr_debug_set(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_OP_WITHOUT_UNABLE(gr_set)
}

static int
_gr_debug_neg(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_OP_WITHOUT_UNABLE(gr_neg)
}

static int
_gr_debug_add(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_OP(gr_add)
}

static int
_gr_debug_sub(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_OP(gr_sub)
}

static int
_gr_debug_mul(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_OP(gr_mul)
}

static int
_gr_debug_div(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_OP(gr_div)
}

static int
_gr_debug_pow(gr_ptr res, gr_srcptr x, gr_srcptr y, gr_ctx_t ctx)
{
    DEBUG_BINARY_OP(gr_pow)
}


static truth_t
_gr_debug_is_invertible(gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_PREDICATE(gr_is_invertible)
}

static int
_gr_debug_inv(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_OP(gr_inv)
}

static int
_gr_debug_sqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_OP(gr_sqrt)
}

static int
_gr_debug_rsqrt(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    DEBUG_UNARY_OP(gr_rsqrt)
}

int _gr_debug_methods_initialized = 0;

gr_static_method_table _gr_debug_methods;

gr_method_tab_input _gr_debug_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_debug_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_debug_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_EXACT,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_THREADSAFE,  (gr_funcptr) _gr_debug_ctx_is_threadsafe},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_debug_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_debug_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_debug_swap},
    {GR_METHOD_SET_SHALLOW,     (gr_funcptr) _gr_debug_set_shallow},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_debug_randtest},
    {GR_METHOD_RANDTEST_SMALL,  (gr_funcptr) _gr_debug_randtest_small},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_debug_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_debug_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_debug_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_debug_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_debug_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_debug_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_debug_set},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_debug_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_debug_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_debug_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_debug_mul},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_debug_div},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_debug_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_debug_inv},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_debug_sqrt},
    {GR_METHOD_RSQRT,           (gr_funcptr) _gr_debug_rsqrt},
    {GR_METHOD_POW,             (gr_funcptr) _gr_debug_pow},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_debug(gr_ctx_t ctx, gr_ctx_t elem_ctx, int flags, double unable_probability)
{
    ctx->which_ring = GR_CTX_DEBUG;

    GR_DEBUG_ELEM_CTX(ctx) = elem_ctx;
    GR_DEBUG_CTX_FLAGS(ctx) = flags;
    GR_DEBUG_CTX_RAND_STATE(ctx) = flint_malloc(sizeof(flint_rand_struct));
    flint_rand_init(GR_DEBUG_CTX_RAND_STATE(ctx));

    if (GR_DEBUG_CTX_FLAGS(ctx) & GR_DEBUG_WRAP)
        ctx->sizeof_elem = sizeof(gr_debug_wrap_struct);
    else
        ctx->sizeof_elem = GR_DEBUG_ELEM_CTX(ctx)->sizeof_elem;

    GR_DEBUG_CTX_UNABLE_PROBABILITY(ctx) = unable_probability;

    ctx->methods = _gr_debug_methods;
    if (!_gr_debug_methods_initialized)
    {
        gr_method_tab_init(_gr_debug_methods, _gr_debug_methods_input);
        _gr_debug_methods_initialized = 1;
    }
}

