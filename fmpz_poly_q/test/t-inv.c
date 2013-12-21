#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("inv... ");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest_not_zero(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_inv(b, a);
        fmpz_poly_q_inv(a, a);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    /* Check a * (1/a) == 1 */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest_not_zero(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_inv(b, a);
        fmpz_poly_q_mul(b, a, b);

        result = fmpz_poly_q_is_one(b) && fmpz_poly_q_is_canonical(b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            fmpz_poly_q_print(a), flint_printf("\n\n");
            fmpz_poly_q_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
