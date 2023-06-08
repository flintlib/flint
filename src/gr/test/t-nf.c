#include "long_extras.h"
#include "fmpq_poly.h"
#include "gr.h"
#include "nf.h"

int main(void)
{
    gr_ctx_t QQa;
    fmpq_poly_t f;
    fmpz_poly_t g;
    slong iter;
    flint_rand_t state;
    int flags = GR_TEST_ALWAYS_ABLE;

    flint_printf("nf....");
    fflush(stdout);

    flint_randinit(state);

    fmpq_poly_init(f);
    fmpz_poly_init(g);

    for (iter = 0; iter < 30; iter++)
    {
        do
        {
            fmpz_poly_randtest_irreducible(g, state, 2 + n_randint(state, 5), 1 + n_randint(state, 10));
        } while (g->length < 2);

        fmpq_poly_set_fmpz_poly(f, g);
        fmpq_poly_scalar_div_ui(f, f, 1 + n_randtest(state) % 256);

        gr_ctx_init_nf(QQa, f);
        gr_test_ring(QQa, 100, flags);
        gr_ctx_clear(QQa);
    }

    fmpq_poly_clear(f);
    fmpz_poly_clear(g);

    flint_randclear(state);

    flint_cleanup_master();
    flint_printf(" PASS\n");
    return 0;
}
