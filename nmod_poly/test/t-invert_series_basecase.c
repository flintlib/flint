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

    Copyright (C) 2011 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    printf("invert_series_basecase....");
    fflush(stdout);

    /* Check result of div against divrem */
    for (i = 0; i < 1000; i++)
    {
        nmod_poly_t q, qinv, prod;
        long m;

        mp_limb_t n;
        do n = n_randtest_not_zero();
        while (!n_is_probabprime(n));

        nmod_poly_init(prod, n);
        nmod_poly_init(qinv, n);
        nmod_poly_init(q, n);
        
        do nmod_poly_randtest(q, n_randint(2000));
        while (q->length == 0 || q->coeffs[0] == 0);

        m = n_randint(q->length) + 1;

        nmod_poly_invert_series_basecase(qinv, q, m);
        
        nmod_poly_mul(prod, q, qinv);
        nmod_poly_truncate(prod, m);

        result = (prod->length == 1 && prod->coeffs[0] == 1);
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(qinv), printf("\n\n");
            nmod_poly_print(prod), printf("\n\n");
            printf("n = %ld\n", n);
            abort();
        }
        
        nmod_poly_clear(q);
        nmod_poly_clear(qinv);
        nmod_poly_clear(prod);
    }

    /* Check aliasing of q and qinv */
    for (i = 0; i < 1000; i++)
    {
        nmod_poly_t q, qinv;
        long m;

        mp_limb_t n;
        do n = n_randtest();
        while (!n_is_probabprime(n));

        nmod_poly_init(q, n);
        nmod_poly_init(qinv, n);
        do nmod_poly_randtest(q, n_randint(1000));
        while (q->length == 0 || q->coeffs[0] == 0);

        m = n_randint(q->length) + 1;

        nmod_poly_invert_series_basecase(qinv, q, m);
        nmod_poly_invert_series_basecase(q, q, m);
        
        result = (nmod_poly_equal(q, qinv));
        if (!result)
        {
            printf("FAIL:\n");
            nmod_poly_print(q), printf("\n\n");
            nmod_poly_print(qinv), printf("\n\n");
            nmod_poly_print(q), printf("\n\n");
            printf("n = %ld, m = %ld\n", n, m);
            abort();
        }

        nmod_poly_clear(q);
        nmod_poly_clear(qinv);
    }

    printf("PASS\n");
    return 0;
}
