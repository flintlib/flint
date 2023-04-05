#include "gr.h"

int main(void)
{
    gr_ctx_t PSL2Z;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("psl2z....");
    fflush(stdout);

    gr_ctx_init_psl2z(PSL2Z);
    gr_test_multiplicative_group(PSL2Z, 10000, flags);
    gr_ctx_clear(PSL2Z);

    flint_cleanup();
    flint_printf(" PASS\n");
    return 0;
}
