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

    printf("scalar_div_mpq... ");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing of a and b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b;
        fmpz_t x1, x2;
        mpq_t y;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_init(x1);
        fmpz_init(x2);
        mpq_init(y);

        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_randtest_not_zero(x1, state, 50);
        fmpz_randtest_not_zero(x2, state, 50);
        fmpz_get_mpz(mpq_numref(y), x1);
        fmpz_get_mpz(mpq_denref(y), x2);
        mpq_canonicalize(y);

        fmpz_poly_q_scalar_div_mpq(a, b, y);
        fmpz_poly_q_scalar_div_mpq(b, b, y);

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
        fmpz_clear(x1);
        fmpz_clear(x2);
        mpq_clear(y);
    }

    /* Check that x (a + b) == x * a + x * b */
    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a, b, c, d;
        fmpz_t x1, x2;
        mpq_t y;

        fmpz_poly_q_init(a);
        fmpz_poly_q_init(b);
        fmpz_poly_q_init(c);
        fmpz_poly_q_init(d);
        fmpz_init(x1);
        fmpz_init(x2);
        mpq_init(y);

        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_randtest(b, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_randtest_not_zero(x1, state, 50);
        fmpz_randtest_not_zero(x2, state, 50);
        fmpz_get_mpz(mpq_numref(y), x1);
        fmpz_get_mpz(mpq_denref(y), x2);
        mpq_canonicalize(y);

        fmpz_poly_q_scalar_div_mpq(c, a, y);
        fmpz_poly_q_scalar_div_mpq(d, b, y);
        fmpz_poly_q_add(d, c, d);

        fmpz_poly_q_add(c, a, b);
        fmpz_poly_q_scalar_div_mpq(c, c, y);

        result = fmpz_poly_q_equal(c, d) && fmpz_poly_q_is_canonical(c);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            printf("b = "), fmpz_poly_q_print(b), printf("\n\n");
            printf("c = "), fmpz_poly_q_print(c), printf("\n\n");
            printf("d = "), fmpz_poly_q_print(d), printf("\n\n");
            gmp_printf("y = %Qd\n\n", y);
            abort();
        }

        fmpz_poly_q_clear(a);
        fmpz_poly_q_clear(b);
        fmpz_poly_q_clear(c);
        fmpz_poly_q_clear(d);
        fmpz_clear(x1);
        fmpz_clear(x2);
        mpq_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
