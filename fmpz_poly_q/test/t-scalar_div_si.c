#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "long_extras.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("scalar_div_si... ");
    fflush(stdout);

    

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        slong x;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        x = z_randtest_not_zero(state);

        fmpz_poly_q_scalar_div_si(a, b, x);
        fmpz_poly_q_scalar_div_si(b, b, x);

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

    /* Check that (a + b) / x == a / x + b / x */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d;
        slong x;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);

        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, 
                                       n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, 
                                       n_randint(state, 50), 50);

        x = z_randtest_not_zero(state);

        fmpz_poly_q_scalar_div_si(c, a, x);
        fmpz_poly_q_scalar_div_si(d, b, x);
        fmpz_poly_q_add(d, c, d);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_scalar_div_si(c, c, x);

        result = fmpz_poly_q_equal(c, d) && fmpz_poly_q_is_canonical(c);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("a = "), fmpz_poly_q_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_q_print(b), flint_printf("\n\n");
            flint_printf("c = "), fmpz_poly_q_print(c), flint_printf("\n\n");
            flint_printf("d = "), fmpz_poly_q_print(d), flint_printf("\n\n");
            flint_printf("x = %wd\n\n", x);
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
