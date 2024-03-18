#include "test_helpers.h"
#include "gr_sparse_vec.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int
test_init(gr_ctx_t ctx)
{
    gr_sparse_vec_t vec;
    int status = GR_SUCCESS;
    gr_sparse_vec_init(vec, 5, ctx);
    gr_sparse_vec_clear(vec, ctx);
    return status;
}


int
test_init_from_entries_sorted(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    gr_sparse_vec_t vec;
    gr_vec_t dvec;
    gr_ptr temp;
    GR_TMP_INIT(temp, ctx);
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong N = 5;
    ulong inds[5] = {0, 2, 4, 6, 8};
    gr_vec_init(dvec, N, ctx);
    status |= _gr_vec_randtest(GR_VEC_ENTRY(dvec, 0, sz), state, N, ctx);
    gr_sparse_vec_init(vec, 2*N, ctx);
    status |= gr_sparse_vec_set_from_entries_sorted_deduped(vec, inds, GR_VEC_ENTRY(dvec, 0, sz), N, ctx);
    for (i = 0; i < 2*N; i++)
    {
        status |= gr_sparse_vec_find_entry(temp, vec, i, ctx);
        if ((i & 1) == 1 && (T_TRUE != gr_is_zero(temp, ctx))) /* Odd inds should all be zero */
            return GR_TEST_FAIL;
        if ((i & 1) == 0 && (T_TRUE != gr_equal(temp, GR_VEC_ENTRY(dvec, i>>1, sz), ctx)))
            return GR_TEST_FAIL;
    }
    GR_TMP_CLEAR(temp, ctx);
    gr_vec_clear(dvec, ctx);
    gr_sparse_vec_clear(vec, ctx);
    return status;
}


int
test_init_from_entries_unsorted_internal(ulong *inds, gr_srcptr entries, slong len, slong n_entries, gr_ctx_t ctx)
{
    /* Yes, this is quadratic time just to be careful */
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong i,j;
    gr_sparse_vec_t vec;
    gr_ptr temp, temp2;
    GR_TMP_INIT2(temp, temp2, ctx);
    gr_sparse_vec_init(vec, len, ctx);
    status |= gr_sparse_vec_set_from_entries(vec, inds, entries, n_entries, ctx);
    if (status != GR_SUCCESS)
        return GR_TEST_FAIL;
    /* Scan through every column in the sparse vector */
    for (i = 0; i < len; i++)
    {
        /* Scan through the entries to figure out what the correct value should be */
        status |= gr_zero(temp, ctx);
        for (j = 0; j < n_entries; j++)
        {
            if (inds[j] == i)
                status |= gr_add(temp, temp, GR_ENTRY(entries, j, sz), ctx);
        }
        /* Check it */
        status |= gr_sparse_vec_find_entry(temp2, vec, i, ctx);
        if (T_TRUE != gr_equal(temp, temp2, ctx))
            return GR_TEST_FAIL;
    }
    GR_TMP_CLEAR(temp, ctx);
    GR_TMP_CLEAR(temp2, ctx);
    gr_sparse_vec_clear(vec, ctx);
    return status;
}



int
test_init_from_entries_unsorted(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong N = 5;
    ulong inds[5] = {8, 4, 3, 8, 1};
    gr_vec_t dvec;
    gr_vec_init(dvec, N, ctx);
    /* First test against random entries */
    status |= _gr_vec_randtest(GR_VEC_ENTRY(dvec, 0, sz), state, N, ctx);
    status |= test_init_from_entries_unsorted_internal(inds, GR_VEC_ENTRY(dvec, 0, sz), 2*N, N, ctx);

    /* Next test against some adversarial entries */
    slong entries_si[5] = {5, 0, 2, -5, 1};
    for (i = 0; i < N; i++)
        status |= gr_set_si(GR_VEC_ENTRY(dvec, i, sz), entries_si[i], ctx);
    status |= test_init_from_entries_unsorted_internal(inds, GR_VEC_ENTRY(dvec, 0, sz), 2*N, N, ctx);

    gr_vec_clear(dvec, ctx);
    return status;
}




TEST_FUNCTION_START(gr_sparse_vec_creation, state)
{   
    gr_ctx_t ctx;
    gr_ctx_init_random(ctx, state);
    CHECK_TEST(test_init(ctx), "Init");
    CHECK_TEST(test_init_from_entries_sorted(state, ctx), "Init from entries sorted");
    CHECK_TEST(test_init_from_entries_unsorted(state, ctx), "Init from entries unsorted");
    gr_ctx_clear(ctx);
    TEST_FUNCTION_END(state);
}
