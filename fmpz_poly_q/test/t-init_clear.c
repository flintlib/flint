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
    FLINT_TEST_INIT(state);

    flint_printf("init/clear... ");
    fflush(stdout);

    

    for (i = 0; i < 100; i++)
    {
        fmpz_poly_q_t a;
        slong len1 = n_randint(state, 50);
        slong len2 = n_randint(state, 50);
        mp_bitcnt_t bits1 = n_randint(state, 50);
        mp_bitcnt_t bits2 = n_randint(state, 50);

        fmpz_poly_q_init(a);
        fmpz_poly_q_randtest(a, state, len1, bits1, len2, bits2);
        fmpz_poly_q_clear(a);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return EXIT_SUCCESS;
}
