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
#include "arith.h"
#include "nmod_vec.h"
#include "ulong_extras.h"

int main(void)
{
    flint_rand_t state;
    mp_ptr b1, b2;
    len_t n;

    const len_t maxn = 3000;

    printf("bell_number_nmod_vec....");
    fflush(stdout);

    flint_randinit(state);

    b1 = _nmod_vec_init(maxn);
    b2 = _nmod_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        nmod_t mod;
        mp_limb_t p;

        do {
            p = n_randtest_prime(state, 0);
        } while (p < n);

        nmod_init(&mod, p);

        arith_bell_number_nmod_vec_recursive(b1, n, mod);
        arith_bell_number_nmod_vec_series(b2, n, mod);

        if (!_nmod_vec_equal(b1, b2, n))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            abort();
        }
    }

    _nmod_vec_clear(b1);
    _nmod_vec_clear(b2);

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
