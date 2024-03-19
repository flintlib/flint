#include <stdlib.h>
#include "gr_sparse_vec.h"

static int
slong_cmp(const void * a, const void * b)
{
    slong ax = *((slong *) a);
    slong bx = *((slong *) b);
    return ax - bx;
}


int
gr_sparse_vec_randtest(gr_sparse_vec_t vec, slong nnz, int replacement, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, sz;
    sz = ctx->sizeof_elem;

    if (nnz < 0 || nnz > vec->length)
        return GR_DOMAIN;

    // Make space
    gr_sparse_vec_fit_nnz(vec, nnz, ctx);

    if (replacement)
    {
        // Randomly sample nnz columns with replacement, and then sort and prune duplicates
        for (i = 0; i < nnz; ++i)
            vec->inds[i] = n_randint(state, vec->length);
        qsort(vec->inds, nnz, sizeof(slong), slong_cmp);

        j = 0;
        for (i = 0; i < nnz; ++i)
            if (i == 0 || vec->inds[i] != vec->inds[i-1])
                vec->inds[j++] = vec->inds[i];
        vec->nnz = j;        
    }
    else
    {
        // Randomly sample nnz columns without replacement, then sort
        for (i = 0; i < vec->length; ++i)
        {
            j = i < nnz ? i : n_randint(state, i+1);
            if (j < nnz) vec->inds[j] = i;
        }
        if (nnz < vec->length)
            qsort(vec->inds, nnz, sizeof(slong), slong_cmp);
        vec->nnz = nnz;
    }

    for (i = 0; i < vec->nnz; ++i)
        status |= gr_randtest_not_zero(GR_ENTRY(vec->nzs, i, sz), state, ctx);

    return status;
}

int
gr_sparse_vec_randtest_prob(gr_sparse_vec_t vec, double prob, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, sz;
    sz = ctx->sizeof_elem;

    if (prob < 0 || prob > 1)
        return GR_DOMAIN;

    // Handle corner cases
    if (prob == 0)
    {
        gr_sparse_vec_zero(vec, ctx);
        return GR_SUCCESS;
    }
    if (prob == 1)
    {
        status |= gr_sparse_vec_randtest(vec, vec->length, 0, state, ctx);
        return status;
    }

    // Allocate space for expected number of nonzeroes, and expand as needed
    gr_sparse_vec_fit_nnz(vec, prob * vec->length, ctx);

    // TODO: for low probability, should be able to do this faster
    vec->nnz = 0;
    for (i = 0; i < vec->length; ++i)
    {
        if (n_randint(state, 0) < 2 * prob * WORD_MAX)
        {
            if (vec->nnz == vec->alloc)
                gr_sparse_vec_fit_nnz(vec, vec->alloc * 2, ctx);
            status |= gr_randtest_not_zero(GR_ENTRY(vec->nzs, vec->nnz, sz), state, ctx);
            vec->inds[vec->nnz++] = i;
        }
    }

    return status;
}
