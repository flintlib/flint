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

    Copyright (C) 2011, 2010 Sebastian Pancratz
    Copyright (C) 2009 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mod_poly.h"
#include "ulong_extras.h"

static int _check(fmpz_mod_poly_struct **B, 
                  const fmpz_mod_poly_t F, const fmpz_mod_poly_t R)
{
    const slong lenF = F->length;
    const slong lenR = R->length;
    const slong N = (lenF - 1) / (lenR - 1);

    slong i;
    int result;

    fmpz_mod_poly_t S;

    fmpz_mod_poly_init(S, &(R->p));
    fmpz_mod_poly_set(S, B[N]);
    for (i = N; i > 0; i--)
    {
        fmpz_mod_poly_mul(S, S, R);
        fmpz_mod_poly_add(S, S, B[i - 1]);
    }
    result = fmpz_mod_poly_equal(F, S);
    fmpz_mod_poly_clear(S);
    return result;
}

int main(void)
{
    int i, result;
    FLINT_TEST_INIT(state);

    flint_printf("radix....");
    fflush(stdout);

    

    for (i = 0; i < 500; i++)
    {
        fmpz_t p;
        fmpz_mod_poly_t f, r;
        fmpz_mod_poly_struct **b;
        fmpz_mod_poly_radix_t D;
        slong j, N;

        fmpz_init(p);
        fmpz_randtest_unsigned(p, state, 2 * FLINT_BITS);
        fmpz_add_ui(p, p, 2);

        fmpz_mod_poly_init(f, p);
        fmpz_mod_poly_init(r, p);
        fmpz_mod_poly_randtest(f, state, n_randint(state, 500));
        do
            fmpz_mod_poly_randtest_not_zero(r, state, n_randint(state, 20) + 2);
        while (r->length < 2);

        N = FLINT_MAX(0, fmpz_mod_poly_degree(f) / fmpz_mod_poly_degree(r));
        b = flint_malloc((N + 1) * sizeof(fmpz_mod_poly_struct *));
        for (j = 0; j <= N; j++)
        {
            b[j] = flint_malloc(sizeof(fmpz_mod_poly_struct));
            fmpz_mod_poly_init(b[j], p);
        }

        /* Ensure that lead(r) is a unit mod p */
        {
            fmpz_t d;
            fmpz *leadR = fmpz_mod_poly_lead(r);

            fmpz_init(d);
            fmpz_gcd(d, p, leadR);
            while (!fmpz_is_one(d))
            {
                fmpz_divexact(leadR, leadR, d);
                fmpz_gcd(d, p, leadR);
            }
            fmpz_clear(d);
        }

        fmpz_mod_poly_radix_init(D, r, f->length - 1 + n_randint(state, 50));

        fmpz_mod_poly_radix(b, f, D);

        result = _check(b, f, r);
        if (!result)
        {
            flint_printf("FAIL:\n");
            flint_printf("result = %d\n", result);
            flint_printf("p = "), fmpz_print(p), flint_printf("\n\n");
            flint_printf("f = "), fmpz_mod_poly_print(f), flint_printf("\n\n");
            flint_printf("r = "), fmpz_mod_poly_print(r), flint_printf("\n\n");
            flint_printf("N = %wd\n\n", N);
            for (j = 0; j <= N; j++)
            {
                flint_printf("b[%wd] = ", j), fmpz_mod_poly_print(b[j]), flint_printf("\n\n");
            }
            abort();
        }

        fmpz_mod_poly_clear(f);
        fmpz_mod_poly_clear(r);
        fmpz_mod_poly_radix_clear(D);
        for (j = 0; j <= N; j++)
        {
            fmpz_mod_poly_clear(b[j]);
            flint_free(b[j]);
        }
        flint_free(b);
        fmpz_clear(p);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

