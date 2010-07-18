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
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int result;
    printf("max_limbs....");
    fflush(stdout);

    _fmpz_vec_randinit();

    for (ulong i = 0; i < 10000UL; i++)
    {
        fmpz * a;
        ulong length, bits, limbs, limbs2;

        length = n_randint(100);

        a = _fmpz_vec_init(length);
        bits  = n_randint(200);
        limbs = (bits + (ulong)(FLINT_BITS - 1)) / FLINT_BITS;
        _fmpz_vec_randtest(a, length, bits);

        limbs2 = _fmpz_vec_max_limbs(a, length);

        result = (limbs >= limbs2);
        if (!result)
        {
            printf("FAIL:\n");
            printf("limbs = %ld, limbs2 = %ld\n", limbs, limbs2);
            abort();
        }

        _fmpz_vec_clear(a, length);
    }

    _fmpz_vec_randclear();

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
