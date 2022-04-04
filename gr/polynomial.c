/* Polynomials over generic rings */

#include "gr.h"
#include "gr_poly.h"

void
polynomial_init(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_init(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int polynomial_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_stream_write(out, "Ring of polynomials over ");
    gr_ctx_write(out, POLYNOMIAL_ELEM_CTX(ctx));
    return GR_SUCCESS;
}

void
polynomial_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
}

truth_t
polynomial_ctx_is_commutative_ring(gr_ctx_t ctx)
{
    return gr_ctx_is_commutative_ring(POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_ctx_is_integral_domain(gr_ctx_t ctx)
{
    return gr_ctx_is_integral_domain(POLYNOMIAL_ELEM_CTX(ctx));
}

void
polynomial_clear(gr_poly_t res, gr_ctx_t ctx)
{
    gr_poly_clear(res, POLYNOMIAL_ELEM_CTX(ctx));
}

void
polynomial_swap(gr_poly_t poly1, gr_poly_t poly2, gr_ctx_t ctx)
{
    gr_poly_swap(poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_write(gr_stream_t out, gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_write(out, poly, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_randtest(gr_poly_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_poly_randtest(res, state, n_randint(state, 5), POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_equal(const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_equal(poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set(gr_poly_t res, const gr_poly_t mat, gr_ctx_t ctx)
{
    return gr_poly_set(res, mat, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_si(gr_poly_t res, slong v, gr_ctx_t ctx)
{
    return gr_poly_set_si(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_ui(gr_poly_t res, ulong v, gr_ctx_t ctx)
{
    return gr_poly_set_ui(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_fmpz(gr_poly_t res, const fmpz_t v, gr_ctx_t ctx)
{
    return gr_poly_set_fmpz(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_set_fmpq(gr_poly_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return gr_poly_set_fmpq(res, v, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_zero(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_zero(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_one(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_one(res, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_neg_one(gr_poly_t res, gr_ctx_t ctx)
{
    return gr_poly_neg_one(res, POLYNOMIAL_ELEM_CTX(ctx));
}

/*
truth_t
polynomial_is_zero(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_zero(poly, POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_is_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_one(poly, POLYNOMIAL_ELEM_CTX(ctx));
}

truth_t
polynomial_is_neg_one(const gr_poly_t poly, gr_ctx_t ctx)
{
    return gr_poly_is_neg_one(poly, POLYNOMIAL_ELEM_CTX(ctx));
}
*/

int
polynomial_neg(gr_poly_t res, const gr_poly_t mat, gr_ctx_t ctx)
{
    return gr_poly_neg(res, mat, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_add(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_add(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_sub(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    return gr_poly_sub(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int
polynomial_mul(gr_poly_t res, const gr_poly_t poly1, const gr_poly_t poly2, gr_ctx_t ctx)
{
    if (POLYNOMIAL_CTX(ctx)->degree_limit != WORD_MAX)
    {
        if (poly1->length != 0 && poly2->length != 0 &&
            poly1->length + poly2->length > POLYNOMIAL_CTX(ctx)->degree_limit)
            return GR_UNABLE;
    }

    return gr_poly_mul(res, poly1, poly2, POLYNOMIAL_ELEM_CTX(ctx));
}

int _gr_poly_methods_initialized = 0;

gr_static_method_table _gr_poly_methods;

gr_method_tab_input _gr_poly_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) polynomial_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) polynomial_ctx_clear},

    {GR_METHOD_CTX_IS_COMMUTATIVE_RING, (gr_funcptr) polynomial_ctx_is_commutative_ring},
    {GR_METHOD_CTX_IS_INTEGRAL_DOMAIN,  (gr_funcptr) polynomial_ctx_is_integral_domain},
    {GR_METHOD_CTX_IS_FIELD,            (gr_funcptr) gr_generic_ctx_predicate_false},

    {GR_METHOD_INIT,        (gr_funcptr) polynomial_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) polynomial_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) polynomial_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) polynomial_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) polynomial_write},
    {GR_METHOD_ZERO,        (gr_funcptr) polynomial_zero},
    {GR_METHOD_ONE,         (gr_funcptr) polynomial_one},
    {GR_METHOD_NEG_ONE,     (gr_funcptr) polynomial_neg_one},
/*
    {GR_METHOD_IS_ZERO,     (gr_funcptr) polynomial_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) polynomial_is_one},
    {GR_METHOD_IS_NEG_ONE,  (gr_funcptr) polynomial_is_neg_one},
*/
    {GR_METHOD_EQUAL,       (gr_funcptr) polynomial_equal},
    {GR_METHOD_SET,         (gr_funcptr) polynomial_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) polynomial_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) polynomial_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) polynomial_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) polynomial_set_fmpq},
    {GR_METHOD_NEG,         (gr_funcptr) polynomial_neg},
    {GR_METHOD_ADD,         (gr_funcptr) polynomial_add},
    {GR_METHOD_SUB,         (gr_funcptr) polynomial_sub},
    {GR_METHOD_MUL,         (gr_funcptr) polynomial_mul},
    {0,                     (gr_funcptr) NULL},
};

void
gr_ctx_init_polynomial(gr_ctx_t ctx, gr_ctx_t base_ring)
{
    ctx->flags = 0;

    if (base_ring->flags & GR_COMMUTATIVE_RING)
        ctx->flags |= GR_COMMUTATIVE_RING;

    ctx->which_ring = GR_WHICH_RING_CUSTOM;

    ctx->sizeof_elem = sizeof(gr_poly_struct);
    ctx->elem_ctx = flint_malloc(sizeof(polynomial_ctx_t));
    ctx->size_limit = WORD_MAX;

    ((polynomial_ctx_t *) ctx->elem_ctx)->base_ring = (gr_ctx_struct *) base_ring;
    ((polynomial_ctx_t *) ctx->elem_ctx)->degree_limit = WORD_MAX;

    ctx->methods = _gr_poly_methods;

    if (!_gr_poly_methods_initialized)
    {
        gr_method_tab_init(_gr_poly_methods, _gr_poly_methods_input);
        _gr_poly_methods_initialized = 1;
    }
}
