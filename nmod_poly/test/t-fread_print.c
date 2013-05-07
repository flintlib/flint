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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

#if !defined (__WIN32) && !defined (__CYGWIN__)

int
main(void)
{
    int i, result, r1;
    flint_rand_t state;
    flint_randinit(state);

    printf("fread_print....");
    fflush(stdout);

    /* Check reading and writing to a file */
    for (i = 0; i < 1000 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b;
        mp_limb_t n = n_randtest_not_zero(state);
        FILE * f = fopen("nmod_poly_test", "w+");

        if (!f)
        {
            printf("Error: unable to open file for writing.\n");
            abort();
        }

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_randtest(a, state, n_randint(state, 100));
        
        nmod_poly_fprint(f, a);
        fflush(f);
        fclose(f);
        f = fopen("nmod_poly_test", "r");
        r1 = nmod_poly_fread(f, b);
        
        result = (r1 && nmod_poly_equal(a, b));
        if (!result)
        {
            printf("FAIL:\n");
            printf("r1 = %d, n = %lu\n", r1, a->mod.n);
            nmod_poly_print(a), printf("\n\n");
            nmod_poly_print(b), printf("\n\n");
            fclose(f);
            remove("nmod_poly_test");
            abort();
        }

        fclose(f);
        if (remove("nmod_poly_test"))
        {
            printf("Error, unable to delete file nmod_poly_test\n");
            abort();
        }
        nmod_poly_clear(a);
        nmod_poly_clear(b);
    }

    flint_randclear(state);

    printf("PASS\n");
    return 0;
}

#else

int main(void)
{
    printf("print/ read....");
    fflush(stdout);
    printf("SKIPPED\n");
    return EXIT_SUCCESS;
}

#endif
