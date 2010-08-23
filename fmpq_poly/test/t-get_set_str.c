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

    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

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

    printf("get_set_str....");
    fflush(stdout);

    fmpq_poly_randinit();

    for (i = 0; i < 10000; i++)
    {
        int ans;
        char * str;
        fmpq_poly_t f, g;

        fmpq_poly_init(f);
        fmpq_poly_init(g);
        fmpq_poly_randtest(f, n_randint(100), n_randint(200));
        
        str = fmpq_poly_get_str(f);
        ans = fmpq_poly_set_str(g, str);

        result = (ans == 0 && fmpq_poly_equal(f, g));
        if (!result)
        {
            printf("FAIL:\n");
            printf("f = "), fmpq_poly_print(f), printf("\n\n");
            printf("g = "), fmpq_poly_print(g), printf("\n\n");
            printf("ans = %d\n\n", ans);
            abort();
        }

        fmpq_poly_clear(f);
        fmpq_poly_clear(g);
        free(str);
    }
    
    fmpq_poly_randclear();
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
