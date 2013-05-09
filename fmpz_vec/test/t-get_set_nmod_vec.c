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

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"
#include "nmod_vec.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("get/set_nmod_vec....");
    fflush(stdout);

    flint_randinit(state);

    /* Check conversion to and from nmod_vec */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        fmpz *a, *b;
        mp_ptr c;
        nmod_t mod;
        len_t i;
        mp_limb_t t;

        len_t len = n_randint(state, 100);
        mp_limb_t n = n_randtest_not_zero(state);

        a = _fmpz_vec_init(len);
        b = _fmpz_vec_init(len);
        c = _nmod_vec_init(len);

        nmod_init(&mod, n);

        _fmpz_vec_randtest(a, state, len, 200);

        _fmpz_vec_get_nmod_vec(c, a, len, mod);
        _fmpz_vec_set_nmod_vec(b, c, len, mod);

        for (i = 0; i < len; i++)
        {
            fmpz_mod_ui(a + i, a + i, n);
            t = fmpz_get_ui(a + i);
            if (t > n / 2)
                fmpz_sub_ui(a + i, a + i, n);
        }

        result = (_fmpz_vec_equal(a, b, len));
        if (!result)
        {
            printf("FAIL:\n");
            _fmpz_vec_print(a, len), printf("\n\n");
            _fmpz_vec_print(b, len), printf("\n\n");
            abort();
        }

        _fmpz_vec_clear(a, len);
        _fmpz_vec_clear(b, len);
        _nmod_vec_clear(c);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
