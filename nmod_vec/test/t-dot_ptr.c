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

    printf("dot_ptr....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        long len;
        nmod_t mod;
        mp_limb_t m, res, res2;
        mp_ptr x, y;
        mp_ptr * z;
        int limbs1;
        long j, offset;

        len = n_randint(state, 1000) + 1;
        m = n_randtest_not_zero(state);
        offset = n_randint(state, 10);

        nmod_init(&mod, m);

        x = _nmod_vec_init(len);
        y = _nmod_vec_init(len);
        z = flint_malloc(sizeof(mp_ptr) * len);

        _nmod_vec_randtest(x, state, len, mod);
        _nmod_vec_randtest(y, state, len, mod);

        for (j = 0; j < len; j++)
            z[j] = &y[j] + offset;

        limbs1 = _nmod_vec_dot_bound_limbs(len, mod);

        res = _nmod_vec_dot_ptr(x, z, -offset, len, mod, limbs1);
        res2 = _nmod_vec_dot(x, y, len, mod, limbs1);

        if (res != res2)
        {
            printf("FAIL:\n");
            printf("m = %lu\n", m);
            printf("len = %ld\n", len);
            printf("limbs1 = %d\n", limbs1);
            abort();
        }

        _nmod_vec_clear(x);
        _nmod_vec_clear(y);
        flint_free(z);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
