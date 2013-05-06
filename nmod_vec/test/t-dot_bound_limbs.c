/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("dot_bound_limbs....");
    fflush(stdout);

    for (i = 0; i < 1000000; i++)
    {
        long len;
        nmod_t mod;
        mp_limb_t m;
        int limbs1, limbs2;
        mpz_t t;

        len = n_randint(state, 10000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        mpz_init2(t, 4*FLINT_BITS);
        mpz_set_ui(t, m-1);
        mpz_mul(t, t, t);
        mpz_mul_ui(t, t, len);
        limbs2 = mpz_size(t);

        if (limbs1 != limbs2)
        {
            printf("FAIL:\n");
            printf("m = %lu\n", m);
            printf("len = %ld\n", len);
            printf("limbs1 = %d\n", limbs1);
            printf("limbs2 = %d\n", limbs2);
            gmp_printf("bound: %Zd\n", t);
            abort();
        }

        mpz_clear(t);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
