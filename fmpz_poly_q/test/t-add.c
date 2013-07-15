#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("add... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and c */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_add(a, a, b);

        result = (fmpz_poly_q_equal(a, c));
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            fmpz_poly_q_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check aliasing of b and c */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_add(b, a, b);

        result = fmpz_poly_q_equal(b, c);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            fmpz_poly_q_print(c), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
    }

    /* Check commutativity */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_add(d, b, a);

        result = fmpz_poly_q_equal(c, d) && fmpz_poly_q_is_canonical(c) 
                                         && fmpz_poly_q_is_canonical(d);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            fmpz_poly_q_print(c), printf("\n\n");
            fmpz_poly_q_print(d), printf("\n\n");
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
