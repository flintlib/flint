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

    Copyright (C) 2015 William Hart
    Copyright (C) 2015 Vladimir Glazachev

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
    FLINT_TEST_INIT(state);
    

    flint_printf("scalar_mul_nmod_shoup....");
    fflush(stdout);

    /* Check (a + b)*c == a*c + b*c */
    for (i = 0; i < 1000000; i++)
    {
        slong len = n_randint(state, 100) + 1;
        mp_limb_t n = n_randtest_not_zero(state) / 2 + 1;
        mp_limb_t c = n_randint(state, n) / 2;
        nmod_t mod;

        mp_ptr vec = _nmod_vec_init(len);
        mp_ptr vec2 = _nmod_vec_init(len);
        mp_ptr vec3 = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _nmod_vec_randtest(vec, state, len, mod);
        _nmod_vec_randtest(vec2, state, len, mod);

        _nmod_vec_add(vec3, vec, vec2, len, mod);
        _nmod_vec_scalar_mul_nmod_shoup(vec3, vec3, len, c, mod);

        _nmod_vec_scalar_mul_nmod_shoup(vec, vec, len, c, mod);
        _nmod_vec_scalar_mul_nmod_shoup(vec2, vec2, len, c, mod);
        _nmod_vec_add(vec, vec, vec2, len, mod);

        result = _nmod_vec_equal(vec, vec3, len);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("len = %wd, n = %wd\n", len, n);
            abort();
        }

        _nmod_vec_clear(vec);
        _nmod_vec_clear(vec2);
        _nmod_vec_clear(vec3);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
