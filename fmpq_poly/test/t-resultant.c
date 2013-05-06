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
#include <gmp.h>
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
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
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
    for (i = 0; i < 50 * flint_test_multiplier(); i++)
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

    /* fredrik's test case */
    {
        fmpq_poly_t f, g;
        fmpq_t x, y;
        int result;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_init(x);
        fmpq_init(y);

        fmpq_poly_set_str(f, "49  16702090503 -23810415210 7561766512 801950253"
             " 56796743 40735271 -15934 820601 -2712604160 -1577466 0 0 -7967 0"
             " 0 0 -14491973 0 6566138489 -55769 0 130523361 4071137 15934"
             " -501921 -59067338 63860755253 23924901 -15934 -262911 -7967"
             " -4389817 0 185876611072 58470888545 130523361 -63736 -130618965"
             " -39835 0 7967 0 55769 -7967 103571 111298990 47802 -3808226"
             " -3800259");

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
        fmpq_poly_set_str(g, "59  -458395/219902324736 151585/4581298432"
           " 112595/219902324736 -2016245/54975581184 0 35/73300774912 0"
           " -234880919/219902324736 7/219902324736 -7/1278501888"
           " -6055/109951162368 7/27487790592 -504623/73300774912"
           " 53673977/219902324736 0 611667/73300774912 -497/13743895296"
           " 0 -6265/219902324736 2446675/73300774912 2345/219902324736"
           " -371/73300774912 -427/6871947648 -3758096377/219902324736"
           " 20595995/109951162368 -256459/73300774912 0 33690223/73300774912"
           " -229369/219902324736 93205/219902324736 -7/107374182"
           " -133/219902324736 -665/13743895296 -146503/219902324736 0"
           " 7/219902324736 66633/73300774912 -855190385/219902324736"
           " 229355/219902324736 0 161/219902324736 887299/219902324736"
           " -427/7582838784 -611667/18325193728 -7/5114007552 833/54975581184"
           " -7/109951162368 -5402264413/219902324736 7/5114007552 35/9162596864"
           " 1133545/219902324736 -151319/73300774912 0 7/219902324736"
           " 7/54975581184 0 -10367/109951162368 7/54975581184 -161/109951162368");
#pragma GCC diagnostic pop

        fmpq_poly_resultant(x, f, g);
        fmpq_poly_resultant(y, g, f);

        if ((fmpq_poly_degree(f) * fmpq_poly_degree(g)) % 2)
            fmpq_neg(y, y);

        result = fmpq_equal(x, y);
        if (!result)
        {
            printf("FAIL (res(f,g) == (-1)^(m * n) res(g, f)):\n");
            printf("x = "), fmpq_print(x), printf("\n\n");
            printf("y = "), fmpq_print(y), printf("\n\n");
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        fmpq_clear(x);
        fmpq_clear(y);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
