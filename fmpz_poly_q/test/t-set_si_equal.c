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

    flint_printf("set_si_equal... ");
    fflush(stdout);

    

    /* Equal polynomials */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        slong n;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        n = z_randtest(state);

        fmpz_poly_q_set_si(a, n);
        fmpz_poly_q_set(b, a);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("n = %wd\n\n", n);
            flint_printf("a = "), fmpz_poly_q_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_q_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        slong m, n;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);

        m = z_randtest(state);
        n = z_randtest(state);
        while (m == n)
            n = z_randtest(state);
        fmpz_poly_q_set_si(a, m);
        fmpz_poly_q_set_si(b, n);

        result = !fmpz_poly_q_equal(a, b);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("m = %wd\n\n", m);
            flint_printf("n = %wd\n\n", n);
            flint_printf("a = "), fmpz_poly_q_print(a), flint_printf("\n\n");
            flint_printf("b = "), fmpz_poly_q_print(b), flint_printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
