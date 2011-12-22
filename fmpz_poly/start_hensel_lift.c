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

#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"

#define TRACE 1

long _fmpz_poly_start_hensel_lift(fmpz_poly_factor_t lifted_fac, long *link, 
    fmpz_poly_t *v, fmpz_poly_t *w, fmpz_poly_t f, 
    const nmod_poly_factor_t local_fac, long target_exp)
{
    const long r = local_fac->num_factors;

    long i;
    fmpz_t p, P, big_P;

    long current_exp = 1;

    fmpz_poly_t monic_f;


    ulong prev_exp;


    fmpz_init(p);
    fmpz_init(P);
    fmpz_init(big_P);

    fmpz_poly_init(monic_f);

    /*
        Set P := p, big_P, monic_f := monic(f).
     */
    fmpz_set_ui(p, (local_fac->factors[0])->mod.n);
    fmpz_set(P, p);
    fmpz_pow_ui(big_P, P, target_exp);

    if (fmpz_is_one(fmpz_poly_lead(f)))
    {
        fmpz_poly_set(monic_f, f);
    }
    else if (fmpz_cmp_si(fmpz_poly_lead(f), -1) == 0)
    {
        fmpz_poly_neg(monic_f, f);
    }
    else
    {
        fmpz_t t;

        fmpz_init(t);
        fmpz_mod(t, fmpz_poly_lead(f), big_P);

        if (fmpz_invmod(t, t, big_P))
        {
            printf("Exception in fmpz_poly_start_hensel_lift.\n");
            abort();
        }

        fmpz_poly_scalar_mul_fmpz(monic_f, f, t);
        fmpz_poly_scalar_mod_fmpz(monic_f, monic_f, big_P);
        fmpz_clear(t);
    }

    #if TRACE
    printf("building hensel tree\n");
    #endif

    fmpz_poly_build_hensel_tree(link, v, w, local_fac);

    /* i = floor(log2(target_exp)) */
    for (i = 0; (1L << (i + 1)) < target_exp; i++) ;

    {
        long num_steps = FLINT_MIN(target_exp, 5 + i);
        long *copy_exponents;
        long *exponents;
        long pow = target_exp;
        long max_steps;

        exponents = malloc(num_steps * sizeof(long));
        copy_exponents = malloc(num_steps * sizeof(long));

        for (i = 0; i < num_steps && pow > 1; i++)
        { 
            copy_exponents[i] = pow;
            pow = (pow + 1)/2;
        }

        max_steps = i;
        exponents[0] = 1;
        for (i = 1; i <= max_steps; i++)
            exponents[i] = copy_exponents[max_steps - i];

        num_steps = max_steps + 1;

        /*
    here num_steps is actually the number of meaningful numbers in the array 
    exponent so num_steps - 2 means that the final time in the loop has 1 and 
    one last time outside of the loop with 0
        */

        #if TRACE
        printf("doing lift\n");
        #endif

        for (i = 0; i < num_steps - 2; i++)
            fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, 
                p, exponents[i], exponents[i + 1], 1);

        /* Last run doesn't calculate the inverses */
        if (num_steps > 1)
            fmpz_poly_tree_hensel_lift(link, v, w, P, monic_f, r, 
                p, exponents[i], exponents[i + 1], 0);

        #if TRACE
        printf("done lifts\n");
        #endif

        prev_exp = exponents[i];

        free(copy_exponents);
        free(exponents);
    }

    fmpz_poly_clear(monic_f);

    /*
        Now everything is lifted to p^target_exp, we just need to 
        insert the factors into their correct places in lifted_fac.
     */
    fmpz_poly_factor_fit_length(lifted_fac, r);

    for (i = 0; i < 2*r - 2; i++)
    { 
        if (link[i] < 0)
        {
            fmpz_poly_scalar_mod_fmpz(lifted_fac->factors + (- link[i] - 1), v[i], P);
            lifted_fac->exponents[- link[i] - 1] = 1L; 
        }
    }
    lifted_fac->length = r;

    fmpz_clear(big_P);
    fmpz_clear(P);
    fmpz_clear(p);

    return prev_exp;
}

