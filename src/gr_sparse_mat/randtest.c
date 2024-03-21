#include <stdlib.h>
#include "gr_sparse_mat.h"

static int
slong_cmp(const void * a, const void * b)
{
    slong ax = *((slong *) a);
    slong bx = *((slong *) b);
    return ax - bx;
}


int gr_coo_mat_randtest(gr_coo_mat_t mat, slong nnz, int replacement, truth_t is_canonical, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, sz, len;
    sz = ctx->sizeof_elem;
    len = flint_mul_sizes(mat->r, mat->c);

    if (nnz < 0 || nnz > len)
        return GR_DOMAIN;

    // Make space
    gr_coo_mat_fit_nnz(mat, nnz, ctx);

    if (replacement)
    {
        // Randomly sample nnz locations with replacement
        for (i = 0; i < nnz; ++i)
            mat->cols[i] = n_randint(state, len);

        // If canonical, sort and compress
        if (is_canonical == T_TRUE)
        {
            qsort(mat->cols, nnz, sizeof(slong), slong_cmp);

            j = 0;
            for (i = 0; i < nnz; ++i)
                if (i == 0 || mat->cols[i] != mat->cols[i-1])
                    mat->cols[j++] = mat->cols[i];
            mat->nnz = j;        
        }
        else
            mat->nnz = nnz;
    }
    else
    {
        // Randomly sample nnz columns without replacement
        for (i = 0; i < len; ++i)
        {
            j = i < nnz ? i : n_randint(state, i+1);
            if (j < nnz) mat->cols[j] = i;
        }
        if (is_canonical == T_TRUE && nnz < len)
            qsort(mat->cols, nnz, sizeof(slong), slong_cmp);
        mat->nnz = nnz;
    }

    for (i = 0; i < nnz; ++i)
    {
        mat->rows[i] = mat->cols[i] / mat->c;
        mat->cols[i] %= mat->c;
    }
    mat->is_canonical = is_canonical;

    for (i = 0; i < mat->nnz; ++i)
    {
        if (is_canonical == T_TRUE)
            status |= gr_randtest_not_zero(GR_ENTRY(mat->nzs, i, sz), state, ctx);
        else
            status |= gr_randtest(GR_ENTRY(mat->nzs, i, sz), state, ctx);
    }
    return status;
}

int
gr_coo_mat_randtest_prob(gr_coo_mat_t mat, double prob, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, sz, len;
    sz = ctx->sizeof_elem;
    len = flint_mul_sizes(mat->r, mat->c);

    if (prob < 0 || prob > 1)
        return GR_DOMAIN;

    // Handle corner cases
    if (prob == 0)
    {
        gr_coo_mat_zero(mat, ctx);
        return GR_SUCCESS;
    }
    if (prob == 1)
    {
        status |= gr_coo_mat_randtest(mat, len, 0, T_TRUE, state, ctx);
        return status;
    }

    // Allocate space for expected number of nonzeroes, and expand as needed
    gr_coo_mat_fit_nnz(mat, prob * len, ctx);

    // TODO: for low probability, should be able to do this faster
    mat->nnz = 0;
    for (i = 0; i < mat->r; ++i)
    {
        for (j = 0; j < mat->c; ++j)
        {
            if (n_randint(state, 0) < 2 * prob * WORD_MAX)
            {
                if (mat->nnz == mat->alloc)
                    gr_coo_mat_fit_nnz(mat, mat->alloc * 2, ctx);
                status |= gr_randtest_not_zero(GR_ENTRY(mat->nzs, mat->nnz, sz), state, ctx);
                mat->rows[mat->nnz] = i;
                mat->cols[mat->nnz++] = j;
            }
        }
    }
    mat->is_canonical = T_TRUE;
    return status;
}
