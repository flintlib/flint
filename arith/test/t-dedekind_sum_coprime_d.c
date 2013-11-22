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
    double s1, s2f;
    fmpz_t hh, kk;
    fmpq_t s2;
    slong h, k;

    FLINT_TEST_INIT(state);

    flint_printf("dedekind_sum_coprime_d....");
    fflush(stdout);

    fmpz_init(hh);
    fmpz_init(kk);
    fmpq_init(s2);

    for (k = 0; k < 300; k++)
    {
        for (h = 0; h <= k; h++)
        {
            if (n_gcd(k, h) == 1)
            {
                fmpz_set_si(hh, h);
                fmpz_set_si(kk, k);

                s1 = arith_dedekind_sum_coprime_d(h, k);
                arith_dedekind_sum_naive(s2, hh, kk);

                s2f = ((double)fmpz_get_si(fmpq_numref(s2))) /
                        fmpz_get_si(fmpq_denref(s2));

                if (fabs(s1 - s2f) > 1e-10)
                {
                    flint_printf("FAIL:\n");
                    flint_printf("s(%wd,%wd)\n", h, k);
                    flint_printf("s1: %.20f\n", s1);
                    flint_printf("s2: %.20f\n", s2f);
                    flint_printf("Exact: "); fmpq_print(s2); flint_printf("\n");
                    abort();
                }
            }
        }
    }

    fmpz_clear(hh);
    fmpz_clear(kk);
    fmpq_clear(s2);

    FLINT_TEST_CLEANUP(state);
    flint_printf("PASS\n");
    return 0;
}
