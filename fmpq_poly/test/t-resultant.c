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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;
    
    printf("resultant....");
    fflush(stdout);

    flint_randinit(state);

    /* Check res(f, g) == (-1)^(deg f deg g) res(g, f) */
    for (i = 0; i < 500; i++)
    {
        fmpq_poly_t f, g;
        fmpq_t x, y;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_init(x);
        fmpq_init(y);

        fmpq_poly_randtest(f, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(g, state, n_randint(state, 60), 60);

        fmpq_poly_resultant(x, f, g);
        fmpq_poly_resultant(y, g, f);
        if ((fmpq_poly_degree(f) * fmpq_poly_degree(g)) % 2)
            fmpq_neg(y, y);

        result = fmpq_equal(x, y);
        if (!result)
        {
            printf("FAIL (res(f,g) == (-1)^(m * n) res(g, f)):\n");
            printf("f = "), fmpq_poly_print(f), printf("\n\n");
            printf("g = "), fmpq_poly_print(g), printf("\n\n");
            printf("x = "), fmpq_print(x), printf("\n\n");
            printf("y = "), fmpq_print(y), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
    }
    
    /* Check res(f h, g) == res(f, g) res(h, g) */
    for (i = 0; i < 500; i++)
    {
        fmpq_poly_t f, g, h;
        fmpq_t x, y, z;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_init(h);
        fmpq_init(x);
        fmpq_init(y);
        fmpq_init(z);

        fmpq_poly_randtest(f, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(g, state, n_randint(state, 60), 60);
        fmpq_poly_randtest(h, state, n_randint(state, 60), 60);

        fmpq_poly_resultant(y, f, g);
        fmpq_poly_resultant(z, h, g);
        fmpq_mul(y, y, z);
        fmpq_poly_mul(f, f, h);
        fmpq_poly_resultant(x, f, g);

        result = fmpq_equal(x, y);
        if (!result)
        {
            printf("FAIL (res(f h, g) == res(f, g) res(h, g)):\n");
            printf("f = "), fmpq_poly_print(f), printf("\n\n");
            printf("g = "), fmpq_poly_print(g), printf("\n\n");
            printf("h = "), fmpq_poly_print(h), printf("\n\n");
            printf("x = "), fmpq_print(x), printf("\n\n");
            printf("y = "), fmpq_print(y), printf("\n\n");
            abort();
        }
        
        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_poly_clear(h);
        fmpq_clear(x);
        fmpq_clear(y);
        fmpq_clear(z);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
