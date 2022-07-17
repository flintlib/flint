#include "gr.h"
#include "flint/fq.h"
#include "flint/fq_poly.h"
#include "flint/fq_mat.h"

#define FQ_CTX(ring_ctx) ((fq_ctx_struct *)((ring_ctx)->elem_ctx))

void
_gr_fq_ctx_clear(gr_ctx_t ctx)
{
    fq_ctx_clear(FQ_CTX(ctx));
    flint_free(ctx->elem_ctx);
}

int
_gr_fq_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Finite field (fq)");
    return GR_SUCCESS;
}

void
_gr_fq_init(fq_t x, const gr_ctx_t ctx)
{
    fq_init(x, FQ_CTX(ctx));
}

void
_gr_fq_clear(fq_t x, const gr_ctx_t ctx)
{
    fq_clear(x, FQ_CTX(ctx));
}

void
_gr_fq_swap(fq_t x, fq_t y, const gr_ctx_t ctx)
{
    fq_t t;
    *t = *x;
    *x = *y;
    *y = *t;
}

int
_gr_fq_randtest(fq_t res, flint_rand_t state, const gr_ctx_t ctx)
{
    fq_randtest(res, state, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_write(gr_stream_t out, const fq_t x, const gr_ctx_t ctx)
{
    gr_stream_write_free(out, fq_get_str_pretty(x, FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_fq_zero(fq_t x, const gr_ctx_t ctx)
{
    fq_zero(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_one(fq_t x, const gr_ctx_t ctx)
{
    fq_one(x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_si(fq_t res, slong v, const gr_ctx_t ctx)
{
    fq_set_si(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_ui(fq_t res, ulong v, const gr_ctx_t ctx)
{
    fq_set_ui(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_set_fmpz(fq_t res, const fmpz_t v, const gr_ctx_t ctx)
{
    fq_set_fmpz(res, v, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_is_zero(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_zero(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_is_one(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_one(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_equal(const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    return fq_equal(x, y, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_set(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    fq_set(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_neg(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    fq_neg(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_add(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_add(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_sub(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_sub(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    fq_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_si(fq_t res, const fq_t x, slong y, const gr_ctx_t ctx)
{
    fq_mul_si(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_ui(fq_t res, const fq_t x, ulong y, const gr_ctx_t ctx)
{
    fq_mul_ui(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_mul_fmpz(fq_t res, const fq_t x, const fmpz_t y, const gr_ctx_t ctx)
{
    fq_mul_fmpz(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_inv(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    if (fq_is_zero(x, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_inv(res, x, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}

int
_gr_fq_div(fq_t res, const fq_t x, const fq_t y, const gr_ctx_t ctx)
{
    if (fq_is_zero(y, FQ_CTX(ctx)))
    {
        return GR_DOMAIN;
    }
    else
    {
        fq_t t;
        fq_init(t, FQ_CTX(ctx));
        fq_inv(t, y, FQ_CTX(ctx));
        fq_mul(res, x, t, FQ_CTX(ctx));
        fq_clear(t, FQ_CTX(ctx));
        return GR_SUCCESS;
    }
}


truth_t
_gr_fq_is_invertible(const fq_t x, const gr_ctx_t ctx)
{
    return (!fq_is_zero(x, FQ_CTX(ctx))) ? T_TRUE : T_FALSE;
}

truth_t
_gr_fq_is_square(const fq_t x, const gr_ctx_t ctx)
{
    return fq_is_square(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_sqrt(fq_t res, const fq_t x, const gr_ctx_t ctx)
{
    if (fq_sqrt(res, x, FQ_CTX(ctx)))
    {
        return GR_SUCCESS;
    }
    else
    {
        return GR_DOMAIN;
    }
}

int
_gr_ctx_fq_prime(fmpz_t p, gr_ctx_t ctx)
{
    fmpz_set(p, fq_ctx_prime(FQ_CTX(ctx)));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_degree(slong * deg, gr_ctx_t ctx)
{
    *deg = fq_ctx_degree(FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_ctx_fq_order(fmpz_t q, gr_ctx_t ctx)
{
    fq_ctx_order(q, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_gen(gr_ptr res, gr_ctx_t ctx)
{
    fq_gen(res, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_frobenius(gr_ptr res, gr_srcptr x, slong e, gr_ctx_t ctx)
{
    fq_frobenius(res, x, e, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_multiplicative_order(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    int ret;
    ret = fq_multiplicative_order(res, x, FQ_CTX(ctx));

    if (ret == 1)
        return GR_SUCCESS;

    /* todo: better solution? */
    return GR_DOMAIN;
}

int
_gr_fq_norm(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_norm(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_trace(fmpz_t res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_trace(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

truth_t
_gr_fq_is_primitive(gr_srcptr x, gr_ctx_t ctx)
{
    return fq_is_primitive(x, FQ_CTX(ctx)) ? T_TRUE : T_FALSE;
}

int
_gr_fq_pth_root(gr_ptr res, gr_srcptr x, gr_ctx_t ctx)
{
    fq_pth_root(res, x, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int
_gr_fq_poly_mullow(fq_struct * res,
    const fq_struct * poly1, slong len1,
    const fq_struct * poly2, slong len2, slong n, gr_ctx_t ctx)
{
    if (len1 >= len2)
        _fq_poly_mullow(res, poly1, len1, poly2, len2, n, FQ_CTX(ctx));
    else
        _fq_poly_mullow(res, poly2, len2, poly1, len1, n, FQ_CTX(ctx));

    return GR_SUCCESS;
}

int
_gr_fq_mat_mul(fq_mat_t res, const fq_mat_t x, const fq_mat_t y, gr_ctx_t ctx)
{
    fq_mat_mul(res, x, y, FQ_CTX(ctx));
    return GR_SUCCESS;
}

int _fq_methods_initialized = 0;

gr_static_method_table _fq_methods;

gr_method_tab_input _fq_methods_input[] =
{
    {GR_METHOD_CTX_CLEAR,       (gr_funcptr) _gr_fq_ctx_clear},
    {GR_METHOD_CTX_WRITE,       (gr_funcptr) _gr_fq_ctx_write},
    {GR_METHOD_CTX_IS_RING,     (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_UNIQUE_FACTORIZATION_DOMAIN,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_FINITE_CHARACTERISTIC,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_ALGEBRAICALLY_CLOSED,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_ORDERED_RING,
                                (gr_funcptr) gr_generic_ctx_predicate_false},
    {GR_METHOD_CTX_IS_EXACT,    (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_CTX_IS_CANONICAL,
                                (gr_funcptr) gr_generic_ctx_predicate_true},
    {GR_METHOD_INIT,            (gr_funcptr) _gr_fq_init},
    {GR_METHOD_CLEAR,           (gr_funcptr) _gr_fq_clear},
    {GR_METHOD_SWAP,            (gr_funcptr) _gr_fq_swap},
    {GR_METHOD_RANDTEST,        (gr_funcptr) _gr_fq_randtest},
    {GR_METHOD_WRITE,           (gr_funcptr) _gr_fq_write},
    {GR_METHOD_ZERO,            (gr_funcptr) _gr_fq_zero},
    {GR_METHOD_ONE,             (gr_funcptr) _gr_fq_one},
    {GR_METHOD_IS_ZERO,         (gr_funcptr) _gr_fq_is_zero},
    {GR_METHOD_IS_ONE,          (gr_funcptr) _gr_fq_is_one},
    {GR_METHOD_EQUAL,           (gr_funcptr) _gr_fq_equal},
    {GR_METHOD_SET,             (gr_funcptr) _gr_fq_set},
    {GR_METHOD_SET_SI,          (gr_funcptr) _gr_fq_set_si},
    {GR_METHOD_SET_UI,          (gr_funcptr) _gr_fq_set_ui},
    {GR_METHOD_SET_FMPZ,        (gr_funcptr) _gr_fq_set_fmpz},
    {GR_METHOD_NEG,             (gr_funcptr) _gr_fq_neg},
    {GR_METHOD_ADD,             (gr_funcptr) _gr_fq_add},
    {GR_METHOD_SUB,             (gr_funcptr) _gr_fq_sub},
    {GR_METHOD_MUL,             (gr_funcptr) _gr_fq_mul},
    {GR_METHOD_MUL_UI,          (gr_funcptr) _gr_fq_mul_ui},
    {GR_METHOD_MUL_SI,          (gr_funcptr) _gr_fq_mul_si},
    {GR_METHOD_MUL_FMPZ,        (gr_funcptr) _gr_fq_mul_fmpz},
    {GR_METHOD_IS_INVERTIBLE,   (gr_funcptr) _gr_fq_is_invertible},
    {GR_METHOD_INV,             (gr_funcptr) _gr_fq_inv},
    {GR_METHOD_DIV,             (gr_funcptr) _gr_fq_div},
    {GR_METHOD_IS_SQUARE,       (gr_funcptr) _gr_fq_is_square},
    {GR_METHOD_SQRT,            (gr_funcptr) _gr_fq_sqrt},

    {GR_METHOD_CTX_FQ_PRIME,            (gr_funcptr) _gr_ctx_fq_prime},
    {GR_METHOD_CTX_FQ_DEGREE,           (gr_funcptr) _gr_ctx_fq_degree},
    {GR_METHOD_CTX_FQ_ORDER,            (gr_funcptr) _gr_ctx_fq_order},
    {GR_METHOD_FQ_GEN,                  (gr_funcptr) _gr_fq_gen},
    {GR_METHOD_FQ_FROBENIUS,            (gr_funcptr) _gr_fq_frobenius},
    {GR_METHOD_FQ_MULTIPLICATIVE_ORDER, (gr_funcptr) _gr_fq_multiplicative_order},
    {GR_METHOD_FQ_NORM,                 (gr_funcptr) _gr_fq_norm},
    {GR_METHOD_FQ_TRACE,                (gr_funcptr) _gr_fq_trace},
    {GR_METHOD_FQ_IS_PRIMITIVE,         (gr_funcptr) _gr_fq_is_primitive},
    {GR_METHOD_FQ_PTH_ROOT,             (gr_funcptr) _gr_fq_pth_root},

    {GR_METHOD_POLY_MULLOW,     (gr_funcptr) _gr_fq_poly_mullow},
    {GR_METHOD_MAT_MUL,         (gr_funcptr) _gr_fq_mat_mul},
    {0,                         (gr_funcptr) NULL},
};

void
gr_ctx_init_fq(gr_ctx_t ctx, const fmpz_t p, slong d, const char * var)
{
    ctx->flags = 0;
    ctx->which_ring = GR_CTX_FQ;
    ctx->sizeof_elem = sizeof(fq_struct);
    ctx->elem_ctx = flint_malloc(sizeof(fq_ctx_struct));
    ctx->size_limit = WORD_MAX;

    fq_ctx_init(FQ_CTX(ctx), p, d, var == NULL ? "a" : var);
    ctx->methods = _fq_methods;

    if (!_fq_methods_initialized)
    {
        gr_method_tab_init(_fq_methods, _fq_methods_input);
        _fq_methods_initialized = 1;
    }
}
