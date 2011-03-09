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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_mat.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_mat.h"

int
main(void)
{
    long m, mod, rep;
    long i;
    flint_rand_t state;
    flint_randinit(state);

    printf("det....");
    fflush(stdout);

    for (rep = 0; rep < 10000; rep++)
    {
        nmod_mat_t A;
        fmpz_mat_t B;
        mp_limb_t Adet;
        fmpz_t Bdet;

        m = n_randint(state, 20);
        mod = n_randtest_prime(state, 0);

        nmod_mat_init(A, m, m, mod);
        fmpz_mat_init(B, m, m);

        nmod_mat_randtest(A, state);

        for (i = 0; i < m*m; i++)
            fmpz_set_ui(&B->entries[i], A->entries[i]);

        Adet = nmod_mat_det(A);

        fmpz_init(Bdet);
        fmpz_mat_det(Bdet, B);
        fmpz_mod_ui(Bdet, Bdet, mod);

        if (Adet != fmpz_get_ui(Bdet))
        {
            printf("FAIL\n");
            abort();
        }

        nmod_mat_clear(A);
        fmpz_mat_clear(B);
        fmpz_clear(Bdet);
    }

    flint_randclear(state);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
