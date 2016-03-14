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

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "aprcl.h"

int main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);
   
    flint_printf("unity_zpq_gauss_sum....");
    fflush(stdout);
    
    for (i = 0; i < 100; i++)
    {
        int result;
        ulong p, q, pnum, ppow;
        fmpz_t n;
        unity_zpq gausssigma, gauss, gausspower;
        n_factor_t factors;

        n_factor_init(&factors);

        q = n_randprime(state, 6, 0);
        if (q == 2)
            q = 7;

        n_factor(&factors, q - 1, 0);

        pnum = n_randint(state, factors.num);
        p = factors.p[pnum];
        ppow = n_randint(state, factors.exp[pnum]);
        if (ppow == 0)
            ppow = 1;

        p = n_pow(p, ppow);

        fmpz_init_set_ui(n, n_randprime(state, 16, 0));

        unity_zpq_init(gausssigma, q, p, n);
        unity_zpq_init(gauss, q, p, n);
        unity_zpq_init(gausspower, q, p, n);

        unity_zpq_gauss_sum(gauss, q, p); 
        unity_zpq_gauss_sum_sigma_pow(gausssigma, q, p);

        unity_zpq_pow(gausspower, gauss, n);

        result = 0;
        for (j = 0; j < p; j++)
        {
            unity_zpq_mul_unity_p_pow(gauss, gausspower, j);
            if (unity_zpq_equal(gauss, gausssigma))
            { 
                result = 1;
                break;
            }
        }
        
        if (result == 0)
        {
            flint_printf("FAIL\n");
            abort();
        }

        unity_zpq_clear(gausssigma);
        unity_zpq_clear(gauss);
        unity_zpq_clear(gausspower);
        fmpz_clear(n);

    }

    FLINT_TEST_CLEANUP(state);

    flint_printf("PASS\n");
    return 0;
}

