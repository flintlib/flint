#include "gr_sparse_vec.h"

static int
int_sorter(const void * a, const void * b)
{
    return ((slong*)a) - ((slong*)b);
}


int
_gr_vec_randtest(gr_sparse_vec_t vec, double density, slong len, flint_rand_t state, gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    slong i, j, sz, nnz;
    sz = ctx->sizeof_elem;
    nnz = (slong)(density * (double)len);
    vec->length = len;
    /* We choose nnz columns and deupe.  For densities above about 1/sqrt(nnz), */
    /* this is not efficient and won't get the density right, but this isn't important for this function */
    if (nnz > vec->alloc)
    {
        vec->cols = flint_realloc(vec->cols, nnz*sizeof(slong));
    }
    for (i = 0; i < nnz; i++)
    {
        vec->cols[i] = n_randint(state, len);
    }
    qsort(vec->cols, nnz, sizeof(slong), int_sorter);
    /* Contract the col list to remove duplicates */
    j = 0;
    for (i = 0; i < nnz; i++)
    {
        if (i == 0 || vec->cols[j-1] != vec->cols[i])
        {
            vec->cols[j] = vec->cols[i];
            j++;
        }
    }
    nnz = j; /* The actual number of distinct nnz */
    /* Make alloc exactly the new nnz to get cols and entries back in sync */
    vec->cols = flint_realloc(vec->cols, nnz*sizeof(slong));
    if (vec->alloc > nnz)
    {
        _gr_vec_clear(GR_ENTRY(vec->entries, nnz, sz), vec->alloc - nnz, ctx);
        vec->entries = flint_realloc(vec->entries, nnz * sz);
    }
    else if (vec->alloc < nnz)
    {
        vec->entries = flint_realloc(vec->entries, nnz * sz);
        _gr_vec_init(GR_ENTRY(vec->entries, vec->alloc, sz), nnz - vec->alloc, ctx);
    }
    vec->alloc = nnz;
    vec->nnz = nnz;
    for (i = 0; i < nnz; i++)
    {
        status |= gr_randtest_not_zero(GR_ENTRY(vec->entries, i, sz), state, ctx);
    }
    return status;
}
