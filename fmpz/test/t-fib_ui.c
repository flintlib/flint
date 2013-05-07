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

    Copyright (C) 2009 William Hart
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
main(void)
{
    long i, n;
    fmpz_t x, y, z, w;

    printf("fib_ui....");
    fflush(stdout);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(z);
    fmpz_init(w);

    /* Twice to check demotion */
    for (n = 0; n < 2; n++)
    {
        for (i = 0; i < 200; i++)
        {
            fmpz_fib_ui(x, i);
            fmpz_fib_ui(y, i+1);
            fmpz_fib_ui(z, i+2);
            fmpz_add(w, x, y);

            if (!fmpz_equal(w, z))
            {
                printf("FAIL: %ld\n", i);
                fmpz_print(x);
                printf("\n");
                fmpz_print(y);
                printf("\n");
                fmpz_print(z);
                printf("\n");
                abort();
            }
        }
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(z);
    fmpz_clear(w);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
