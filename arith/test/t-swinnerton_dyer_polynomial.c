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
#include "profiler.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "fmpq_poly.h"

static const mp_limb_t known_values[] =
{
    2147483629UL,
    1073742093UL,
    1342248677UL,
    3319936736UL,
    2947821228UL,
    1019513834UL,
    3324951530UL,
    1995039408UL,
    3505683295UL,
    3567639420UL,
    394942914UL
};

int main()
{
    fmpz_poly_t S;
    mp_limb_t r;
    long n;

    printf("swinnerton_dyer_polynomial....");
    fflush(stdout);

    for (n = 0; n <= 10; n++)
    {
        fmpz_poly_init(S);
        arith_swinnerton_dyer_polynomial(S, n);
        r = fmpz_poly_evaluate_mod(S, 2147483629UL, 4294967291UL);

        if (r != known_values[n])
        {
            printf("ERROR: wrong evaluation of S_%ld\n", n);
            abort();
        }

        fmpz_poly_clear(S);
    }

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
