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

    printf("dot....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        long len;
        nmod_t mod;
        mp_limb_t m, res;
        mp_ptr x, y;
        int limbs1;
        mpz_t s, t;
        long j;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        res = _nmod_vec_dot(x, y, len, mod, limbs1);

        mpz_init(s);
        mpz_init(t);

        for (j = 0; j < len; j++)
        {
            mpz_set_ui(t, x[j]);
            mpz_addmul_ui(s, t, y[j]);
        }

        mpz_mod_ui(s, s, m);

        if (mpz_get_ui(s) != res)
        {
            printf("FAIL:\n");
            printf("m = %lu\n", m);
            printf("len = %ld\n", len);
            printf("limbs1 = %d\n", limbs1);
            abort();
        }

        mpz_clear(s);
        mpz_clear(t);

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
