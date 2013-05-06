#include <stdlib.h>
#include <stdio.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("zero... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a;

        fmpz_poly_q_init(a);
        fmpz_poly_q_randtest(a, state, n_randint(state, 50), 50, n_randint(state, 50), 50);
        fmpz_poly_q_zero(a);

        result = fmpz_poly_q_is_zero(a) && fmpz_poly_q_is_canonical(a);
        if (!result)
        {
            printf("FAIL:\n");
            printf("a = "), fmpz_poly_q_print(a), printf("\n\n");
            abort();
        }

        fmpz_poly_q_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
