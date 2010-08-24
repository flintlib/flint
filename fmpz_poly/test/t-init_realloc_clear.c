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

******************************************************************************/

#include <stdio.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    fmpz_randstate_t state;

    printf("init/init2/realloc/clear....");
    fflush(stdout);

    fmpz_poly_randinit(state);

    for (i = 0; i < 10000; i++)
    {
        fmpz_poly_t a;

        fmpz_poly_init2(a, n_randint(100));
        fmpz_poly_clear(a);
    }

    for (i = 0; i < 10000; i++)
    {
        fmpz_poly_t a;

        fmpz_poly_init2(a, n_randint(100));
        fmpz_poly_realloc(a, n_randint(100));
        fmpz_poly_realloc(a, n_randint(100));
        fmpz_poly_clear(a);
    }

    for (i = 0; i < 10000; i++)
    {
        fmpz_poly_t a;

        fmpz_poly_init(a);
        fmpz_poly_randtest(a, state, n_randint(100), n_randint(200));
        fmpz_poly_clear(a);
    }

    fmpz_poly_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
