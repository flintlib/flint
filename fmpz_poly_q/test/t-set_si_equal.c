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
    flint_rand_t state;

    printf("set_si_equal... ");
    fflush(stdout);

    flint_randinit(state);

    /* Equal polynomials */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        len_t n;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        n = z_randtest(state);

        fmpz_poly_q_set_si(a, n);
        fmpz_poly_q_set(b, a);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            printf("n = %ld\n\n", n);
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_q_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        len_t m, n;

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
            printf("FAIL:\n");
            printf("m = %ld\n\n", m);
            printf("n = %ld\n\n", n);
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_q_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
