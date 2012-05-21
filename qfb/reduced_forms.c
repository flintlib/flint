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

int pow_incr(int * pows, int * exp, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        pows[i]++;
        if (pows[i] > exp[i])
            pows[i] = 0;
        else
            return 1;
    }

    return 0;
}

long qfb_reduced_forms(qfb ** forms, long d)
{
    long a, b, c, p, blim, alloc, num, sqrt, i, prod, prime_i;
    mp_limb_t b2, ppow, exp;
    n_factor_t * fac;
    mp_limb_t * s;

    if (d >= 0)
    {
       printf("Exception: qfb_reduced_forms not implemented for positive discriminant.\n");
       abort();
    }

    blim = n_sqrt(-d/3); /* maximum a value to check */

    (*forms) = NULL; 
    alloc = 0;
    
    if ((-d & 3) == 2 || (-d & 3) == 1) /* ensure d is 0, 1 mod 4 */
        return 0;

    fac = flint_calloc(blim + 1, sizeof(n_factor_t));

    sqrt = n_sqrt(blim*blim - d);
    n_compute_primes(FLINT_MAX(sqrt, 10000));

    prime_i = 1;
    while ((p = flint_primes[prime_i]) <= sqrt) /* sieve for factors of p^exp */
    {
        ppow = p;
        exp = 1;
        while (ppow <= blim*blim - d) /* loop through powers of p */
        {
            num = n_sqrtmod_primepow(&s, n_negmod((-d) % ppow, ppow), p, exp);

            if (exp == 1)
            {
                for (i = 0; i < num; i++) /* sieve with each sqrt mod p */
                {
                     long off = s[i];
                     while (off <= blim)
                     {
                         fac[off].p[fac[off].num] = p;
                         fac[off].exp[fac[off].num] = 1;
                         fac[off].num++;
                         off += ppow;
                     }
                }
            } else /* exp > 1 */
            {
                for (i = 0; i < num; i++) /* sieve with each sqrt mod p^exp */
                {
                     long off = s[i];
                     while (off <= blim)
                     {
                         fac[off].exp[fac[off].num - 1]++;
                         off += ppow;
                     }
                }
            }

            free(s);

            exp++;
            ppow *= p;
        }

        prime_i++;
    }

    for (b = (d & 1); b <= blim; b += 2) /* write any remaining factors, including 2^exp */
    {
        b2 = (b*b - d)/4;

        count_trailing_zeros(exp, b2); /* powers of 2 */
        if (exp)
        {
            fac[b].p[fac[b].num] = 2;
            fac[b].exp[fac[b].num] = exp;
            fac[b].num++;
        }
        
        prod = 1;
        for (i = 0; i < fac[b].num; i++)
            prod *= n_pow(fac[b].p[i], fac[b].exp[i]);

        b2 /= prod;
        
        if (b2 != 1)
        {
            fac[b].p[fac[b].num] = b2;
            fac[b].exp[fac[b].num] = 1;
            fac[b].num++;     
        }
    }

    num = 0;
    for (b = (d & 1); b <= blim; b += 2) /* compute forms for each b */
    {
         int pows[FLINT_MAX_FACTORS_IN_LIMB];
         int n = fac[b].num;

         b2 = (b*b - d)/4;

         for (i = 0; i < n; i++)
             pows[i] = 0;

         do
         {
             a = 1;
             
             for (i = 0; i < n; i++)
                 a *= n_pow(fac[b].p[i], pows[i]);
             
             c = b2 / a;
             
             if (a <= c && b <= a) /* we have a form */
             {
                 mp_limb_t g;
                 
                 if (b)
                 {   
                     g= n_gcd(c, b);
                     g = n_gcd(a, g);
                 } else
                     g = n_gcd(c, a);

                 if (g == 1) /* primitive form */
                 {
                     if (num == alloc) /* realloc if necessary */
                     {
                         (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(blim, 100))*sizeof(qfb));
                         alloc += FLINT_MIN(blim, 100);
                     }
                    
                     (*forms)[num].a = a; /* record form */
                     (*forms)[num].b = b;
                     (*forms)[num++].c = c;
                     
                     if (b && a != c && b != a)
                     {
                         if (num == alloc) /* realloc if necessary */
                         {
                             (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(blim, 100))*sizeof(qfb));
                             alloc += FLINT_MIN(blim, 100);
                         }
                    
                        (*forms)[num].a = a; /* record opposite form */
                        (*forms)[num].b = -b;
                        (*forms)[num++].c = c;
                     }
                 }
             }
         } while (pow_incr(pows, fac[b].exp, n));
    }

    flint_free(fac);

    return num;
}
