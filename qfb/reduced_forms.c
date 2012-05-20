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

    Copyright (C) 2012 William Hart

******************************************************************************/

#undef ulong /* prevent clash with stdlib */
#include <stdlib.h>
#define ulong unsigned long
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

long compute_c(long a, long b, long d)
{
    fmpz_t fb;

    fmpz_init(fb);

    fmpz_set_si(fb, b);

    fmpz_mul_si(fb, fb, b);
    
    if (d > 0) 
        fmpz_sub_ui(fb, fb, d);
    else
        fmpz_add_ui(fb, fb, -d);
    
    fmpz_fdiv_q_ui(fb, fb, a);
    fmpz_fdiv_q_2exp(fb, fb, 2);
    
    b = fmpz_get_si(fb);
    
    fmpz_clear(fb);

    return b;
}

long qfb_reduced_forms(qfb ** forms, long d)
{
    long a, j, p, alim, alloc, num, roots, sqrt, i, prod, prime_i;
    mp_limb_t a2;
    n_factor_t * fac;

    if (d >= 0)
    {
       printf("Exception: qfb_reduced_forms not implemented for positive discriminant.\n");
       abort();
    }

    alim = n_sqrt(-d/3); /* maximum a value to check */

    (*forms) = NULL; 
    alloc = 0;
    
    if ((-d & 4) == 2 || (-d & 4) == 3) /* ensure d is 0, 1 mod 4 */
        return 0;

    fac = flint_calloc(alim + 1, sizeof(n_factor_t));

    for (a = 2; a <= alim; a += 2) /* find powers of 2 dividing 4a values */
    {
        a2 = a;
        fac[a].exp[0] = n_remove(&a2, 2) + 2;
        fac[a].p[0] = 2;
        fac[a].num = 1;
    }

    for (a = 1; a <= alim; a += 2)
    {
        fac[a].exp[0] = 2;
        fac[a].p[0] = 2;
        fac[a].num = 1;
    }
    
    sqrt = n_sqrt(alim);
    n_compute_primes(FLINT_MAX(sqrt, 10000));

    prime_i = 1;
    while ((p = flint_primes[prime_i]) <= sqrt) /* sieve for factors of 4a values */
    {
        for (a = p; a <= alim; a+= p)
        {
            a2 = a;
            num = fac[a].num;
            fac[a].exp[num] = n_remove2_precomp(&a2, p, flint_prime_inverses[prime_i]);
            fac[a].p[num] = p;
            fac[a].num++;
        }
        prime_i++;
    }
    
    for (a = 1; a <= alim; a++) /* write any remaining prime factor of each 4a value */
    {
        prod = 1;
        for (i = 0; i < fac[a].num; i++)
            prod *= n_pow(fac[a].p[i], fac[a].exp[i]);
        p = (4*a)/prod;
        if (p != 1)
        {
            num = fac[a].num;
            fac[a].exp[num] = 1;
            fac[a].p[num] = p;
            fac[a].num++;
        }
    }

    num = 0;

    for (a = 1; a <= alim; a++) /* loop through possible a's */
    {
        mp_limb_t * s;

        roots = n_sqrtmodn(&s, n_negmod((-d)%(4*a), 4*a), fac + a);

        for (j = 0; j < roots; j++) /* loop through all square roots of d mod 4a */
        {
           mp_limb_signed_t b = s[j];
           
           if (b > 2*a) b -= 4*a;
           
           if (-a < b && b <= a) /* we may have a form */
           {
               mp_limb_t c = compute_c(a, b, d);
               
               if (c >= a && (b >= 0 || a != c)) /* we have a form */
               {
                   mp_limb_t g;
                   
                   if (b)
                   {   
                      g= n_gcd(c, FLINT_ABS(b));
                      g = n_gcd(a, g);
                   } else
                      g = n_gcd(c, a);

                   if (g == 1) /* we have a primitive form */
                   {
                       if (num == alloc) /* realloc if necessary */
                       {
                           (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(alim, 100))*sizeof(qfb));
                           alloc += FLINT_MIN(alim, 100);
                       }

                       (*forms)[num].a = a; /* record form */
                       (*forms)[num].b = b;
                       (*forms)[num++].c = c;
                   }
               }
           }
        }

        free(s);
    }

    flint_free(fac);

    return num;
}
