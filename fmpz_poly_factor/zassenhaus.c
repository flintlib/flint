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

    Copyright (C) 2011 Andy Novocin
    Copyright (C) 2011 Sebastian Pancratz
   
******************************************************************************/

#include <stdlib.h>
#include "fmpz_poly_factor.h"

void _fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t final_fac, 
                                  long exp, const fmpz_poly_t f, long cutoff)
{
    const long lenF = f->length;

    const long M_bits = fmpz_bits(fmpz_poly_lead(f)) 
                        + FLINT_ABS(fmpz_poly_max_bits(f)) 
                        + lenF + FLINT_CLOG2(lenF);
    long a;

    long i, num_primes;

    mp_limb_t p = 2;
    mp_limb_t min_p = p; 
    nmod_poly_factor_t fac;
    fmpz_poly_factor_t lifted_fac;
    long min_r;

    if (lenF <= 1)
    {
        return;
    }
    if (lenF == 2)
    {
        fmpz_poly_factor_insert(final_fac, f, exp);
        return;
    }

    nmod_poly_factor_init(fac);
    min_r = lenF;
    i = 0;

    for (num_primes = 1; num_primes < 3; num_primes++)
    {
        int check = 1;

        for ( ; i < 200 && check == 1; i++, p = n_nextprime(p, 0))
	    {
            nmod_poly_t F_d, F_sbo, F_tmp;

            nmod_poly_init(F_tmp, p);
            nmod_poly_init(F_d, p);
            nmod_poly_init(F_sbo, p);

            fmpz_poly_get_nmod_poly(F_tmp, f);
            if (F_tmp->length == lenF)
		    {
                nmod_poly_derivative(F_d, F_tmp);
                nmod_poly_gcd(F_sbo, F_tmp, F_d);

                if (nmod_poly_is_one(F_sbo))
                {
                    nmod_poly_factor_t temp_fac;

                    nmod_poly_factor_init(temp_fac);
                    nmod_poly_factor(temp_fac, F_tmp);

                    if (temp_fac->num <= min_r)
                    {
                        min_r = temp_fac->num;
                        min_p = p;
                        nmod_poly_factor_clear(fac);
                        nmod_poly_factor_init(fac);
                        nmod_poly_factor_concat(fac, temp_fac);
                    }
                    nmod_poly_factor_clear(temp_fac);
                    check = 0;
                }
            }
            nmod_poly_clear(F_d);
            nmod_poly_clear(F_sbo);
            nmod_poly_clear(F_tmp);
        }
   
	    if (i == 200)
	    {
            printf("Warning (fmpz_poly_factor_zassenhaus): \n");
            printf("Polynomial was not square_free after 200 primes, \n");
            printf("maybe an error?\n");

            nmod_poly_factor_clear(fac);
		    return;
        }
    }
    p = min_p;

    /* Check various abort conditions */
    {
        const long r = fac->num;

        if (r == 0)
        {
            printf("Exception (fmpz_poly_factor_zassenhaus): r = 0\n");
            nmod_poly_factor_clear(fac);
            abort();
        }
        if (r == 1)
        {
            fmpz_poly_factor_insert(final_fac, f, exp);
            nmod_poly_factor_clear(fac);
            return;
        }
        if (r > cutoff)
        {
            printf("Exception (fmpz_poly_factor_zassenhaus): r > cutoff\n");
            nmod_poly_factor_clear(fac);
            abort();
        }
    }

    fmpz_poly_factor_init(lifted_fac);

    a = M_bits / FLINT_CLOG2(p);

    /*
        Begin Hensel Lifting phase, make the tree in v, w, and link.
        Later we'll do a check if use_Hoeij_Novocin (try for smaller a).
     */
    {
        const long r = fac->num;
        long *link;
        fmpz_poly_t *v, *w;

        link = malloc((2*r-2) * sizeof(long));
        v    = malloc((2*r-2) * sizeof(fmpz_poly_t));
        w    = malloc((2*r-2) * sizeof(fmpz_poly_t));

        for (i = 0; i < 2*r - 2; i++)
        {
            fmpz_poly_init(v[i]);
            fmpz_poly_init(w[i]);
        }

        printf("Going with p = %lu to the a = %ld, r = %ld\n", p, a, r);
        _fmpz_poly_hensel_start_lift(lifted_fac, link, v, w, f, fac, a);

        for (i = 0; i < 2*r - 2; i++)
        {
            fmpz_poly_clear(v[i]);
            fmpz_poly_clear(w[i]);
        }
        free(link);
        free(v);
        free(w);
    }

    /* Recombination */
    {
        fmpz_t P;
        fmpz_init(P);
        fmpz_set_ui(P, p);
        fmpz_pow_ui(P, P, a);

        fmpz_poly_factor_zassenhaus_recombination(final_fac, lifted_fac, f, P, exp);

        fmpz_clear(P);
    }

    fmpz_poly_factor_clear(lifted_fac);
    nmod_poly_factor_clear(fac);
}

void fmpz_poly_factor_zassenhaus(fmpz_poly_factor_t fac, const fmpz_poly_t G)
{
    const long lenG = G->length;
    fmpz_poly_t g;

    if (lenG == 0)
    {
        fmpz_set_ui(&fac->c, 0);
        return;
    }

    if (lenG == 1)
    {
        fmpz_set(&fac->c, G->coeffs);
        return;
    }

    fmpz_poly_init(g);

    if (lenG == 2)
    {
        fmpz_poly_content(&fac->c, G);
        fmpz_poly_scalar_divexact_fmpz(g, G, &fac->c);
        fmpz_poly_factor_insert(fac, g, 1);
    }
    else
    {
        long j, k;
        fmpz_poly_factor_t sq_fr_fac;

        /* Does a presearch for a factor of form x^k */
        for (k = 0; fmpz_is_zero(G->coeffs + k); k++) ;

        if (k != 0)
        {
            fmpz_poly_t t;

            fmpz_poly_init(t);
            fmpz_poly_set_coeff_ui(t, 1, 1);
            fmpz_poly_factor_insert(fac, t, k);
            fmpz_poly_clear(t);
        }

        fmpz_poly_shift_right(g, G, k);

        /* Could make other tests for x-1 or simple things 
           maybe take advantage of the composition algorithm */
        fmpz_poly_factor_init(sq_fr_fac);
        fmpz_poly_factor_squarefree(sq_fr_fac, g);

        /* Factor each square-free part */
        for (j = 0; j < sq_fr_fac->num; j++)
            _fmpz_poly_factor_zassenhaus(fac, sq_fr_fac->exp[j], sq_fr_fac->p + j, 10);

        fmpz_poly_factor_clear(sq_fr_fac);
    }
    fmpz_poly_clear(g);
}

