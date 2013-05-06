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

    Copyright (C) 2009 William Hart

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
    int i, result;
    flint_rand_t state;
    flint_randinit(state);

    printf("reduce....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        long j, len = n_randint(state, 100) + 1;
        mp_ptr vec = _nmod_vec_init(len);
        mp_ptr vec2 = _nmod_vec_init(len);

        mp_limb_t n = n_randtest_not_zero(state);
        nmod_t mod;
        nmod_init(&mod, n);

        for (j = 0; j < len; j++)
        {
            vec[j] = n_randtest(state);
            vec2[j] = vec[j];
        }

        _nmod_vec_reduce(vec, vec, len, mod);
        for (j = 0; j < len; j++)
            vec2[j] = n_mod2_preinv(vec2[j], mod.n, mod.ninv);

        result = _nmod_vec_equal(vec, vec2, len);
        if (!result)
        {
            printf("FAIL:\n");
            printf("len = %ld, n = %ld\n", len, n);
            abort();
        }

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
