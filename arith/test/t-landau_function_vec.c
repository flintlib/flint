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

static const mp_limb_t known[] = {
    1, 1, 2, 3, 4, 6, 6, 12, 15, 20, 30, 30, 60, 60, 84, 105, 140, 210,
    210, 420, 420, 420, 420, 840, 840, 1260, 1260, 1540, 2310, 2520,
    4620, 4620, 5460, 5460, 9240, 9240, 13860, 13860, 16380, 16380,
    27720, 30030, 32760, 60060, 60060, 60060, 60060, 120120
};

int main(void)
{
    fmpz * res;
    len_t k, n;

    printf("landau_function_vec....");
    fflush(stdout);

    n = 45;
    res = _fmpz_vec_init(n);
    arith_landau_function_vec(res, n);

    for (k = 0; k < n; k++)
    {
        if (fmpz_cmp_ui(res + k, known[k]))
        {
            printf("FAIL:\n");
            printf("k = %ld, res[k] = %ld, expected: %ld\n",
                k, fmpz_get_si(res + k), known[k]);
            abort();
        }
    }

    _fmpz_vec_clear(res, n);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
