#include "test_helpers.h"
#include "gr_sparse_vec.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int
test_init(gr_ctx_t ctx)
{
    gr_sparse_vec_t vec;
    int status = GR_SUCCESS;
    status |= gr_sparse_vec_init(vec, 5, ctx);
    status |= gr_sparse_vec_clear(vec, ctx);
    return status;
}

int
test_init_from_entries(flint_rand_t state, gr_ctx_t ctx)
{
    slong i;
    gr_sparse_vec_t vec;
    gr_vec_t dvec;
    gr_ptr temp;
    GR_TMP_INIT(temp, ctx);
    int status = GR_SUCCESS;
    slong sz = ctx->sizeof_elem;
    slong N = 5;
    ulong cols[5] = {0, 2, 4, 6, 8};
    status |= gr_vec_init(dvec, N, ctx);
    status |= _gr_vec_randtest(GR_VEC_ENTRY(dvec, 0, sz), state, N, ctx);
    status |= gr_sparse_vec_init(vec, 2*N, ctx);
    status |= gr_sparse_vec_set_from_entries_sorted_deduped(vec, cols, GR_VEC_ENTRY(dvec, 0, sz), N, ctx);
    for (i = 0; i < 2*N; i++)
    {
        status |= gr_sparse_vec_find_entry(temp, vec, i, ctx);
        if ((i & 1) == 1 && (T_TRUE != gr_is_zero(temp, ctx))) /* Odd cols should all be zero */
            return GR_TEST_FAIL;
        if ((i & 1) == 0 && (T_TRUE != gr_equal(temp, GR_VEC_ENTRY(dvec, i>>1, ctx))))
            return GR_TEST_FAIL;
    }
    GR_TMP_CLEAR(temp, ctx);
    status |= gr_vec_clear(dvec, ctx);
    status |= gr_sparse_vec_clear(vec, ctx);
    return status;
}


TEST_FUNCTION_START(gr_sparse_vec_creation, state)
{   
    gr_ctx_t ctx;
    gr_ctx_init_random(ctx, state);
    CHECK_TEST(test_init(ctx), "Init");
    CHECK_TEST(test_init_from_entries(state, ctx), "Init from entries");
    gr_ctx_clear(ctx);
    TEST_FUNCTION_END(state);
}
