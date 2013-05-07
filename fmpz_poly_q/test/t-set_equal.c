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

    printf("set/equal... ");
    fflush(stdout);

    flint_randinit(state);

    /* Equal polynomials */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_set(b, a);

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

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        long coeff = n_randint(state, 50);
        fmpz_t x1, x2;

        fmpz_init(x1);
        fmpz_init(x2);
        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);

        fmpz_poly_q_set(b, a);

        fmpz_poly_get_coeff_fmpz(x2, fmpz_poly_q_numref(b), coeff);
        do
            fmpz_randtest(x1, state, 50);
        while (fmpz_equal(x1, x2));
        fmpz_poly_set_coeff_fmpz(fmpz_poly_q_numref(b), coeff, x1);
        fmpz_poly_q_canonicalise(b);

        result = !fmpz_poly_q_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            fmpz_poly_q_print(a), printf("\n\n");
            fmpz_poly_q_print(b), printf("\n\n");
            abort();
        }

        fmpz_clear(x1);
        fmpz_clear(x2);
        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
