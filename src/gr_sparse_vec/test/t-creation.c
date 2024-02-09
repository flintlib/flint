#include "test_helpers.h"
#include "gr_sparse_vec.h"

#define CHECK_TEST(x, name) { if (GR_SUCCESS != (x)) { flint_printf("FAIL %s\n", (name)); flint_abort(); } }

int
test_init(gr_ctx_t ctx)
{
    gr_sparse_vec_t vec;
    int status = GR_SUCCESS;
    status |= gr_sparse_vec_init(vec, 5, ctx);
    status |= gr_sparse_vec_clear(vec);
    return status;
}

int
test_init_from_entries(gr_ctx_t ctx)
{
    int status = GR_SUCCESS;
    return status;
}


TEST_FUNCTION_START(gr_sparse_vec_creation, state)
{   
    gr_ctx_t ctx;
    gr_ctx_init_random(ctx, state);
    CHECK_TEST(test_init(ctx), "Init");
    CHECK_TEST(test_init_from_entries(ctx), "Init from entries");
    gr_ctx_clear(ctx);
    TEST_FUNCTION_END(state);
}