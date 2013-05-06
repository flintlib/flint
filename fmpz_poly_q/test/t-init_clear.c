/******************************************************************************

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <gmp.h>

#include "fmpz_poly_q.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;

    printf("init/clear... ");
    fflush(stdout);

    flint_randinit(state);

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a;
        long len1 = n_randint(state, 50);
        long len2 = n_randint(state, 50);
        mp_bitcnt_t bits1 = n_randint(state, 50);
        mp_bitcnt_t bits2 = n_randint(state, 50);

        fmpz_poly_q_init(a);
        fmpz_poly_q_randtest(a, state, len1, bits1, len2, bits2);
        fmpz_poly_q_clear(a);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return EXIT_SUCCESS;
}
