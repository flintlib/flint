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
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"


int main(void)
{
    ulong k;
    fmpz_t x;
    fmpz_t y;

    printf("primorial....");
    fflush(stdout);

    fmpz_init(x);
    fmpz_init(y);
    fmpz_set_ui(y, 1);

    for (k = 0; k < 10000; k++)
    {
       arith_primorial(x, k);
       if (n_is_prime(k))
          fmpz_mul_ui(y, y, k);
       if (!fmpz_equal(x, y))
       {
          printf("FAIL:\n");
          printf("primorial of %lu disagrees with direct product\n", k); 
          fmpz_print(x);
          printf("\n");
          abort();
       }
    }

    fmpz_clear(x);
    fmpz_clear(y);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
