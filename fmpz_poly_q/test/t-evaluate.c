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

    printf("evaluate...");
    fflush(stdout);

    flint_randinit(state);

    /* Check aliasing */
    for (i = 0; i < 100; i++)
    {
        int ans1, ans2;
        mpq_t a, b;
        fmpz_t num, den;
        fmpz_poly_q_t f;

        mpq_init(a);
        mpq_init(b);
        fmpz_init(num);
        fmpz_init(den);
        fmpz_poly_q_init(f);
        fmpz_poly_q_randtest(f, state, n_randint(state, 10), 10, n_randint(state, 10), 10);

        fmpz_randtest(num, state, 50);
        fmpz_randtest_not_zero(den, state, 50);
        fmpz_get_mpz(mpq_numref(a), num);
        fmpz_get_mpz(mpq_denref(a), den);
        mpq_canonicalize(a);

        ans1 = fmpz_poly_q_evaluate(b, f, a);
        ans2 = fmpz_poly_q_evaluate(a, f, a);

        result = (ans1 == ans2) && mpq_equal(a, b);
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpz_poly_q_print(f), printf("\n");
            printf("num = "), fmpz_print(num), printf("\n");
            printf("den = "), fmpz_print(den), printf("\n");
            gmp_printf("a = %Qd\n", a);
            gmp_printf("b = %Qd\n", b);
            printf("ans1 = %d\n", ans1);
            printf("ans2 = %d\n", ans2);
            abort();
        }

        mpq_clear(a);
        mpq_clear(b);
        fmpz_clear(num);
        fmpz_clear(den);
        fmpz_poly_q_clear(f);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
