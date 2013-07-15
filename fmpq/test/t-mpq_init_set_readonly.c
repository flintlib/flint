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

    Copyright (C) 2010, 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "long_extras.h"
#include "fmpz.h"
#include "fmpq.h"

int main(void)
{
    int i;
    flint_rand_t state;
    
    printf("mpq_init_set_readonly....");
    fflush(stdout);
    
    flint_randinit(state);

    /* Create some small fmpq rationals, clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, FLINT_BITS - 2);

        flint_mpq_init_set_readonly(z, f);
        flint_mpq_clear_readonly(z);
    }

    /* Create some large fmpq rationals, do not clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, 2 * FLINT_BITS);

        if (COEFF_IS_MPZ(*fmpq_numref(f)) && COEFF_IS_MPZ(*(fmpq_denref(f))))
        {
            flint_mpq_init_set_readonly(z, f);
        }

        fmpq_clear(f);
    }

    /* Create some more fmpq rationals */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f, g;
        mpq_t z;

        fmpq_init(f);
        fmpq_init(g);
        fmpq_randtest(f, state, 2 * FLINT_BITS);

        flint_mpq_init_set_readonly(z, f);
        fmpq_set_mpq(g, z);

        if (!fmpq_equal(f, g))
        {
            printf("FAIL:\n\n");
            printf("f = "), fmpq_print(f), printf("\n");
            printf("g = "), fmpq_print(g), printf("\n");
            gmp_printf("z = %Qd\n", z);
        }

        fmpq_clear(f);
        fmpq_clear(g);
        flint_mpq_clear_readonly(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
