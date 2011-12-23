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
    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "nmod_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    flint_rand_t state;

    printf("start_continue_hensel_lift....");
    fflush(stdout);

    flint_randinit(state);

    /* We check that lifting local factors of F_poly yields factors */
    for (i = 0; i < 1000; i++)
    {
        fmpz_poly_t F_poly, F_poly2, F_poly3, R;
        nmod_poly_t fac;
        nmod_poly_factor_t f_fac;
        fmpz_poly_factor_t F_fac;
        long bits, length, nbits, n, exp, j, part_exp;

        long r;
        fmpz_poly_t *v, *w;
        long *link;
        long prev_exp;

        bits = n_randint(state, 200) + 1;
        nbits = n_randint(state, FLINT_BITS - 6) + 6;

        fmpz_poly_init(F_poly);
        fmpz_poly_init(F_poly2);
        fmpz_poly_init(F_poly3);
        fmpz_poly_init(R);

        nmod_poly_factor_init(f_fac);
        fmpz_poly_factor_init(F_fac);

        n = n_randprime(state, nbits, 0); 
        exp = bits/(FLINT_BIT_COUNT(n) - 1) + 1;
        part_exp = n_randint(state, exp);

        nmod_poly_init(fac, n);

        do {
            length = n_randint(state, 200) + 2;  
            do { fmpz_poly_randtest(F_poly2, state, length, bits); } 
            while (F_poly2->length < 2);

            fmpz_set_ui(F_poly2->coeffs, n_randbits(state, FLINT_MIN(bits, FLINT_BITS - 2)));
            fmpz_set_ui(F_poly2->coeffs + F_poly2->length - 1, 1);

            length = n_randint(state, 200) + 2;  
            do { fmpz_poly_randtest(F_poly3, state, length, bits); } 
            while (F_poly3->length < 2);

            fmpz_set_ui(F_poly3->coeffs, n_randbits(state, FLINT_MIN(bits, FLINT_BITS - 2)));
            fmpz_set_ui(F_poly3->coeffs + F_poly3->length - 1, 1);

            fmpz_poly_mul(F_poly, F_poly2, F_poly3);

            fmpz_poly_get_nmod_poly(fac, F_poly);
        } while (!nmod_poly_is_squarefree(fac));

        fmpz_poly_get_nmod_poly(fac, F_poly2);
        nmod_poly_factor_insert(f_fac, fac, 1);

        fmpz_poly_get_nmod_poly(fac, F_poly3);
        nmod_poly_factor_insert(f_fac, fac, 1);
        nmod_poly_clear(fac);

        r = f_fac->num;
        v = malloc((2*r - 2)*sizeof(fmpz_poly_t));
        w = malloc((2*r - 2)*sizeof(fmpz_poly_t));
        link = malloc((2*r - 2)*sizeof(long));

        for(j = 0; j < 2*r - 2; j++)
        {
            fmpz_poly_init(v[j]);
            fmpz_poly_init(w[j]);
        }

        if (part_exp < 1)
        {
            _fmpz_poly_hensel_start_lift(F_fac, link, v, w, F_poly, f_fac, exp);
        }
        else
        {
            fmpz_t nn;

            fmpz_init_set_ui(nn, n);

            prev_exp = _fmpz_poly_hensel_start_lift(F_fac, link, v, w, 
                F_poly, f_fac, part_exp);
            _fmpz_poly_hensel_continue_lift(F_fac, link, v, w, 
                F_poly, prev_exp, part_exp, exp, nn);

            fmpz_clear(nn);
        }

        result = 1;
        for (j = 0; j < F_fac->num; j++)
        {
            fmpz_poly_rem(R, F_poly, F_fac->p + j);
            result &= (R->length == 0);
        }

        for (j = 0; j < 2*r - 2; j++)
        {
            fmpz_poly_clear(v[j]);
            fmpz_poly_clear(w[j]);
        }

        free(link);
        free(v);
        free(w);

        if (!result) 
        {
            printf("FAIL:\n");
            printf("length = %ld, bits = %ld, n = %ld, exp = %ld\n", length, bits, n, exp);
            fmpz_poly_print(F_poly); printf("\n\n");
            fmpz_poly_print(F_poly2); printf("\n\n");
            fmpz_poly_print(F_poly3); printf("\n\n");
            fmpz_poly_factor_print(F_fac); printf("\n\n");
            abort();
        } 

        nmod_poly_factor_clear(f_fac);
        fmpz_poly_factor_clear(F_fac);

        fmpz_poly_clear(R);

        fmpz_poly_clear(F_poly3);
        fmpz_poly_clear(F_poly2);
        fmpz_poly_clear(F_poly);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
