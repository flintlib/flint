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

    printf("derivative... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_derivative(b, a);
        fmpz_poly_q_derivative(a, a);

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

    /* Check constants have derivative zero */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);

        fmpz_poly_q_set_si(a, z_randtest(state));

        fmpz_poly_q_derivative(b, a);

        result = fmpz_poly_q_is_zero(b);
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

    /* Check (f g)' = f' g + f g' */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d, lhs, rhs;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_poly_q_init(lhs);
        fmpz_poly_q_init(rhs);
        fmpz_poly_q_randtest(a, state, n_randint(state, 20), 20, n_randint(state, 20), 20);
        fmpz_poly_q_randtest(b, state, n_randint(state, 20), 20, n_randint(state, 20), 20);

        fmpz_poly_q_mul(lhs, a, b);
        fmpz_poly_q_derivative(lhs, lhs);
        fmpz_poly_q_derivative(c, a);
        fmpz_poly_q_derivative(d, b);
        fmpz_poly_q_mul(c, c, b);
        fmpz_poly_q_mul(d, a, d);
        fmpz_poly_q_add(rhs, c, d);

        result = fmpz_poly_q_equal(lhs, rhs) && fmpz_poly_q_is_canonical(lhs) 
                                             && fmpz_poly_q_is_canonical(rhs);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            fmpz_poly_q_print(lhs), printf("\n\n");
            fmpz_poly_q_print(rhs), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
        fmpz_poly_q_clear(lhs);
        fmpz_poly_q_clear(rhs);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
