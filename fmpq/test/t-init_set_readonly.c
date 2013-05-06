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
    
    printf("init_set_readonly....");
    fflush(stdout);
    
    flint_randinit(state);

    /* Create some small fmpq rationals, clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, FLINT_BITS - 2);

        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g;

            fmpq_init_set_readonly(g, z);
            fmpq_clear_readonly(g);
        }

        mpq_clear(z);
    }

    /* Create some small fmpq ratioals, do *not* clear the mpq_t */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, FLINT_BITS - 2);

        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g;

            fmpq_init_set_readonly(g, z);
        }

        mpq_clear(z);
    }

    /* Create some more fmpq rationals */
    for (i = 0; i < 100000; i++)
    {
        fmpq_t f;
        mpq_t z;

        fmpq_init(f);
        fmpq_randtest(f, state, 2 * FLINT_BITS);
        mpq_init(z);
        fmpq_get_mpq(z, f);

        {
            fmpq_t g, h;

            fmpq_init_set_readonly(g, z);
            fmpq_init(h);
            fmpq_set_mpq(h, z);

            if (!fmpq_equal(g, h))
            {
                printf("FAIL:\n\n");
                printf("g = "), fmpq_print(g), printf("\n");
                printf("h = "), fmpq_print(h), printf("\n");
                gmp_printf("z = %Qd\n", z);
            }

            fmpq_clear_readonly(g);
            fmpq_clear(h);
        }

        fmpq_clear(f);
        mpq_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}

