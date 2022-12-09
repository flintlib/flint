/* Vectors over generic rings */

#include "gr.h"
#include "gr_vec.h"
#include "gr_poly.h"

#define ENTRY_CTX(ctx) (VECTOR_CTX(ctx)->base_ring)

int
_gr_vec_check_resize(gr_vec_t res, slong n, gr_ctx_t ctx)
{
    if (VECTOR_CTX(ctx)->all_sizes)
    {
        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
        return GR_SUCCESS;
    }
    else
    {
        if (n != VECTOR_CTX(ctx)->n)
            return GR_DOMAIN;

        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
        return GR_SUCCESS;
    }
}

void
vector_gr_vec_init(gr_vec_t res, gr_ctx_t ctx)
{
    gr_vec_init(res, VECTOR_CTX(ctx)->n, ENTRY_CTX(ctx));
}

int vector_gr_vec_ctx_write(gr_stream_t out, gr_ctx_t ctx)
{
    gr_ctx_ptr elem_ctx = ENTRY_CTX(ctx);

    if (VECTOR_CTX(ctx)->all_sizes)
    {
        gr_stream_write(out, "Vectors (any length) over ");
    }
    else
    {
        gr_stream_write(out, "Space of length ");
        gr_stream_write_si(out, VECTOR_CTX(ctx)->n);
        gr_stream_write(out, " vectors over ");
    }

    gr_ctx_write(out, elem_ctx);
    return GR_SUCCESS;
}

/* todo: public */
truth_t gr_ctx_vector_gr_vec_is_fixed_size(gr_ctx_t ctx)
{
    return (VECTOR_CTX(ctx)->all_sizes) ? T_FALSE : T_TRUE;
}

void
vector_gr_vec_clear(gr_vec_t res, gr_ctx_t ctx)
{
    gr_vec_clear(res, ENTRY_CTX(ctx));
}

void
vector_gr_vec_swap(gr_vec_t vec1, gr_vec_t vec2, gr_ctx_t ctx)
{
    gr_poly_swap((gr_poly_struct *) vec1, (gr_poly_struct *) vec2, ENTRY_CTX(ctx));
}

int
vector_gr_vec_write(gr_stream_t out, gr_vec_t vec, gr_ctx_t ctx)
{
    slong i;

    gr_stream_write(out, "[");
    for (i = 0; i < vec->length; i++)
    {
        gr_write(out, gr_vec_entry_ptr(vec, i, ENTRY_CTX(ctx)), ENTRY_CTX(ctx));
        if (i < vec->length - 1)
            gr_stream_write(out, ", ");
    }
    gr_stream_write(out, "]");
    return GR_SUCCESS;
}

int
vector_gr_vec_randtest(gr_vec_t res, flint_rand_t state, gr_ctx_t ctx)
{
    slong i, n;
    int status;

    if (VECTOR_CTX(ctx)->all_sizes)
    {
        n = n_randint(state, 7);
        gr_vec_set_length(res, n, ENTRY_CTX(ctx));
    }
    else
        n = res->length;


    status = GR_SUCCESS;
    for (i = 0; i < n; i++)
        status |= gr_randtest(gr_vec_entry_ptr(res, i, ENTRY_CTX(ctx)), state, ENTRY_CTX(ctx));

    return status;
}

truth_t
vector_gr_vec_equal(const gr_vec_t vec1, const gr_vec_t vec2, gr_ctx_t ctx)
{
    return gr_poly_equal((gr_poly_struct *) vec1, (gr_poly_struct *) vec2, ENTRY_CTX(ctx));
}

int
vector_gr_vec_set(gr_vec_t res, const gr_vec_t vec, gr_ctx_t ctx)
{
    return gr_vec_set(res, vec, ENTRY_CTX(ctx));
}

/* todo: convert from matrices, ...? */
int
vector_gr_vec_set_other(gr_vec_t res, gr_srcptr x, gr_ctx_t x_ctx, gr_ctx_t ctx)
{
    if (x_ctx == ctx)
    {
        return vector_gr_vec_set(res, x, ctx);
    }
    else if (x_ctx->which_ring == GR_CTX_GR_VEC)
    {
        const gr_vec_struct * xvec = x;
        slong i, j;
        int status = GR_SUCCESS;
        slong sz, xsz;

        if (res->length != xvec->length)
        {
            if (VECTOR_CTX(ctx)->all_sizes)
                gr_vec_set_length(res, xvec->length, ENTRY_CTX(ctx));
            else
                return GR_DOMAIN;
        }

        sz = ENTRY_CTX(ctx)->sizeof_elem;
        xsz = ENTRY_CTX(x_ctx)->sizeof_elem;

        for (i = 0; i < res->length; i++)
        {
            for (j = 0; j < xvec->length; j++)
            {
                status = gr_set_other(GR_ENTRY(res->entries, i, sz),
                            GR_ENTRY(xvec->entries, i, xsz),
                            ENTRY_CTX(x_ctx),
                            ENTRY_CTX(ctx));

                if (status != GR_SUCCESS)
                    return status;
            }
        }

        return GR_SUCCESS;
    }

    return GR_UNABLE;
}

int
vector_gr_vec_neg(gr_vec_t res, const gr_vec_t src, gr_ctx_t ctx)
{
    return gr_poly_neg((gr_poly_struct *) res, (gr_poly_struct *) src, ENTRY_CTX(ctx));
}


int _gr_vec_methods_initialized = 0;

gr_static_method_table _gr_vec_methods;

gr_method_tab_input _gr_vec_methods_input[] =
{
    {GR_METHOD_CTX_WRITE,   (gr_funcptr) vector_gr_vec_ctx_write},
    {GR_METHOD_INIT,        (gr_funcptr) vector_gr_vec_init},
    {GR_METHOD_CLEAR,       (gr_funcptr) vector_gr_vec_clear},
    {GR_METHOD_SWAP,        (gr_funcptr) vector_gr_vec_swap},
    {GR_METHOD_RANDTEST,    (gr_funcptr) vector_gr_vec_randtest},
    {GR_METHOD_WRITE,       (gr_funcptr) vector_gr_vec_write},
    {GR_METHOD_EQUAL,       (gr_funcptr) vector_gr_vec_equal},
    {GR_METHOD_SET,         (gr_funcptr) vector_gr_vec_set},
    {GR_METHOD_SET_OTHER,   (gr_funcptr) vector_gr_vec_set_other},
    {GR_METHOD_NEG,         (gr_funcptr) vector_gr_vec_neg},
    {0,                     (gr_funcptr) NULL},
};

void
_gr_ctx_init_vector(gr_ctx_t ctx, gr_ctx_t base_ring, int all_sizes, slong n)
{
    ctx->which_ring = GR_CTX_GR_VEC;
    ctx->sizeof_elem = sizeof(gr_vec_struct);
    ctx->size_limit = WORD_MAX;

    if (n < 0) flint_abort();

    VECTOR_CTX(ctx)->base_ring = (gr_ctx_struct *) base_ring;
    VECTOR_CTX(ctx)->all_sizes = all_sizes;
    VECTOR_CTX(ctx)->n = n;

    ctx->methods = _gr_vec_methods;

    if (!_gr_vec_methods_initialized)
    {
        gr_method_tab_init(_gr_vec_methods, _gr_vec_methods_input);
        _gr_vec_methods_initialized = 1;
    }
}

void gr_ctx_init_vector_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring)
{
    _gr_ctx_init_vector(ctx, base_ring, 1, 0);
}

void gr_ctx_init_vector_space_gr_vec(gr_ctx_t ctx, gr_ctx_t base_ring, slong n)
{
    _gr_ctx_init_vector(ctx, base_ring, 0, n);
}
