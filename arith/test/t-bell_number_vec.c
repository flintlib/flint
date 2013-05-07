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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int main(void)
{
    fmpz * b1;
    fmpz * b2;
    long n;

    const long maxn = 1000;

    printf("bell_number_vec....");
    fflush(stdout);

    b1 = _fmpz_vec_init(maxn);
    b2 = _fmpz_vec_init(maxn);

    for (n = 0; n < maxn; n += (n < 50) ? + 1 : n/4)
    {
        arith_bell_number_vec_recursive(b1, n);
        arith_bell_number_vec_multi_mod(b2, n);

        if (!_fmpz_vec_equal(b1, b2, n))
        {
            printf("FAIL:\n");
            printf("n = %ld\n", n);
            abort();
        }
    }

    _fmpz_vec_clear(b1, maxn);
    _fmpz_vec_clear(b2, maxn);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
