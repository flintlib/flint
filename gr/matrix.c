/* Matrices over generic rings */

#include "gr.h"
#include "gr_mat.h"

void
matrix_init(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_init(res, MATRIX_CTX(ctx)->n, MATRIX_CTX(ctx)->n, MATRIX_CTX(ctx)->base_ring);
}

int matrix_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    slong n = MATRIX_CTX(ctx)->n;
    gr_ctx_ptr elem_ctx = MATRIX_CTX(ctx)->base_ring;

    gr_stream_write(out, "Ring of ");
    gr_stream_write_si(out, n);
    gr_stream_write(out, " x ");
    gr_stream_write_si(out, n);
    gr_stream_write(out, " matrices over ");
    gr_ctx_write(out, elem_ctx);
    return GR_SUCCESS;
}

void
matrix_ctx_clear(gr_ctx_t ctx)
{
    flint_free(ctx->elem_ctx);
}

void
matrix_clear(gr_mat_t res, gr_ctx_t ctx)
{
    gr_mat_clear(res, MATRIX_CTX(ctx)->base_ring);
}

void
matrix_swap(gr_mat_t mat1, gr_mat_t mat2, gr_ctx_t ctx)
{
    gr_mat_swap(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_write(gr_stream_t out, gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_write(out, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_randtest(gr_mat_t res, flint_rand_t state, gr_ctx_t ctx)
{
    return gr_mat_randtest(res, state, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_equal(const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_equal(mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_set(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_si(gr_mat_t res, slong v, gr_ctx_t ctx)
{
    return gr_mat_set_si(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_ui(gr_mat_t res, ulong v, gr_ctx_t ctx)
{
    return gr_mat_set_ui(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_fmpz(gr_mat_t res, const fmpz_t v, gr_ctx_t ctx)
{
    return gr_mat_set_fmpz(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_set_fmpq(gr_mat_t res, const fmpq_t v, gr_ctx_t ctx)
{
    return gr_mat_set_fmpq(res, v, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_zero(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_zero(res, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_one(gr_mat_t res, gr_ctx_t ctx)
{
    return gr_mat_one(res, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_zero(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_zero(mat, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_one(mat, MATRIX_CTX(ctx)->base_ring);
}

truth_t
matrix_is_neg_one(const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_is_neg_one(mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_neg(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_neg(res, mat, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_add(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_add(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_sub(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_sub(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_mul(gr_mat_t res, const gr_mat_t mat1, const gr_mat_t mat2, gr_ctx_t ctx)
{
    return gr_mat_mul_classical(res, mat1, mat2, MATRIX_CTX(ctx)->base_ring);
}

int
matrix_inv(gr_mat_t res, const gr_mat_t mat, gr_ctx_t ctx)
{
    return gr_mat_inv(res, mat, MATRIX_CTX(ctx)->base_ring);
}


int _gr_mat_methods_initialized = 0;

gr_static_method_table _gr_mat_methods;

gr_method_tab_input _gr_mat_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) matrix_ctx_write},
    {GR_METHOD_CTX_CLEAR,   (gr_funcptr) matrix_ctx_clear},
    {GR_METHOD_INIT,        (gr_funcptr) matrix_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) matrix_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) matrix_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) matrix_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) matrix_write},
    {GR_METHOD_ZERO,        (gr_funcptr) matrix_zero},
    {GR_METHOD_ONE,         (gr_funcptr) matrix_one},
    {GR_METHOD_IS_ZERO,     (gr_funcptr) matrix_is_zero},
    {GR_METHOD_IS_ONE,      (gr_funcptr) matrix_is_one},
    {GR_METHOD_IS_NEG_ONE,  (gr_funcptr) matrix_is_neg_one},
    {GR_METHOD_EQUAL,       (gr_funcptr) matrix_equal},
    {GR_METHOD_SET,         (gr_funcptr) matrix_set},
    {GR_METHOD_SET_UI,      (gr_funcptr) matrix_set_ui},
    {GR_METHOD_SET_SI,      (gr_funcptr) matrix_set_si},
    {GR_METHOD_SET_FMPZ,    (gr_funcptr) matrix_set_fmpz},
    {GR_METHOD_SET_FMPQ,    (gr_funcptr) matrix_set_fmpq},
    {GR_METHOD_NEG,         (gr_funcptr) matrix_neg},
    {GR_METHOD_ADD,         (gr_funcptr) matrix_add},
    {GR_METHOD_SUB,         (gr_funcptr) matrix_sub},
    {GR_METHOD_MUL,         (gr_funcptr) matrix_mul},
    {GR_METHOD_INV,         (gr_funcptr) matrix_inv},
    {0,                     (gr_funcptr) NULL},
};

/* rename: gr_ctx_init_gr_mat */
void
gr_ctx_init_matrix(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    ctx->flags = 0;

    if (base_ring->flags & GR_FINITE_RING)
        ctx->flags |= GR_FINITE_RING;

    ctx->which_ring = GR_WHICH_RING_CUSTOM;

    ctx->sizeof_elem = sizeof(gr_mat_struct);
    ctx->elem_ctx = flint_malloc(sizeof(matrix_ctx_t));
    ctx->size_limit = WORD_MAX;

    ((matrix_ctx_t *) ctx->elem_ctx)->base_ring = (gr_ctx_struct *) base_ring;
    ((matrix_ctx_t *) ctx->elem_ctx)->n = n;

    ctx->methods = _gr_mat_methods;

    if (!_gr_mat_methods_initialized)
    {
        gr_method_tab_init(_gr_mat_methods, _gr_mat_methods_input);
        _gr_mat_methods_initialized = 1;
    }
}
