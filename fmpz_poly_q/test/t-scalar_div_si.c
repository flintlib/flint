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

    printf("scalar_div_si... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        len_t x;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        x = z_randtest_not_zero(state);

        fmpz_poly_q_scalar_div_si(a, b, x);
        fmpz_poly_q_scalar_div_si(b, b, x);

        result = fmpz_poly_q_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    /* Check that (a + b) / x == a / x + b / x */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d;
        len_t x;

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
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_q_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_q_print(c), printf("\n\n");
            printf("d = "), fmpz_poly_q_print(d), printf("\n\n");
            printf("x = %ld\n\n", x);
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
