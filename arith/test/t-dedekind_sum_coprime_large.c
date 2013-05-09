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
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arith.h"
#include "ulong_extras.h"
#include "math.h"

int main(void)
{
    fmpz_t hh, kk;
    fmpq_t s1, s2;
    len_t h, k;

    printf("dedekind_sum_coprime_large....");
    fflush(stdout);

    fmpz_init(hh);
    fmpz_init(kk);
    fmpq_init(s1);
    fmpq_init(s2);

    for (k = 0; k < 300; k++)
    {
        for (h = 0; h <= k; h++)
        {
            if (n_gcd(k, h) == 1)
            {
                fmpz_set_si(hh, h);
                fmpz_set_si(kk, k);

                arith_dedekind_sum_coprime_large(s1, hh, kk);
                arith_dedekind_sum_naive(s2, hh, kk);

                if (!fmpq_equal(s1, s2))
                {
                    printf("FAIL:\n");
                    printf("s(%ld,%ld)\n", h, k);
                    printf("s1: "); fmpq_print(s1); printf("\n");
                    printf("s2: "); fmpq_print(s2); printf("\n");
                    abort();
                }
            }
        }
    }

    fmpz_clear(hh);
    fmpz_clear(kk);
    fmpq_clear(s1);
    fmpq_clear(s2);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
