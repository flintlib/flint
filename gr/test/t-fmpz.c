#include "gr.h"

int main()
{
    gr_ctx_t ZZ;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("fmpz....");
    fflush(stdout);

    gr_ctx_init_fmpz(ZZ);
    ZZ->size_limit = 1000;
    gr_test_ring(ZZ, 10000, flags);

    {
        fmpz * a, *b;
        a = gr_heap_init(ZZ);
        b = gr_heap_init(ZZ);

        fmpz_set_str(a, "1000000000000000000000", 10);
        fmpz_set_str(b, "1000000000000000000001", 10);

        if (gr_sub(b, b, a, ZZ) != GR_SUCCESS || gr_is_one(b, ZZ) != T_TRUE)
            flint_abort();

        gr_heap_clear(a, ZZ);
        gr_heap_clear(b, ZZ);
    }

    gr_ctx_clear(ZZ);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
