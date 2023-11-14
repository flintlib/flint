/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "qfb.h"

/*
   Iterate through all factors of a number given factorisation
   into n prime powers whose maximum values are stored in exp,
   storing the values at the current iteration in pows.
*/
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

slong qfb_reduced_forms_large(qfb ** forms, slong d)
{
    slong a, j, k, p, alim, alloc, num, roots, sqrt, i, prod, prime_i;
    mp_srcptr primes;
    const double * prime_inverses;
    mp_limb_t a2;
    n_factor_t * fac;

    if (d >= 0)
        flint_throw(FLINT_ERROR, "%s not implemented for positive discriminant\n", __func__);

    alim = n_sqrt(-d/3); /* maximum a value to check */

    (*forms) = NULL;
    alloc = 0;

    if ((-d & 3) == 2 || (-d & 3) == 1) /* ensure d is 0, 1 mod 4 */
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
    primes = n_primes_arr_readonly(FLINT_MAX(sqrt, 10000));
    prime_inverses = n_prime_inverses_arr_readonly(FLINT_MAX(sqrt, 10000));

    prime_i = 1;
    while ((p = primes[prime_i]) <= sqrt) /* sieve for factors of 4a values */
    {
        for (a = p; a <= alim; a+= p)
        {
            a2 = a;
            num = fac[a].num;
            fac[a].exp[num] = n_remove2_precomp(&a2, p, prime_inverses[prime_i]);
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
               /*
                  let B = FLINT_BITS
                  -sqrt(2^(B-1)/3) < b < sqrt(2^(B-1)/3)
                  0 < -d < 2^(B-1)
               */
               mp_limb_t c = ((mp_limb_t) (b*b) + (mp_limb_t) (-d))/(4*(mp_limb_t) a);

               if (c >= (mp_limb_t) a && (b >= 0 || a != c)) /* we have a form */
               {
                   mp_limb_t g;

                   if (b)
                   {
                      g = n_gcd(c, FLINT_ABS(b));
                      g = n_gcd(a, g);
                   } else
                      g = n_gcd(c, a);

                   if (g == 1) /* we have a primitive form */
                   {
                       if (num == alloc) /* realloc if necessary */
                       {
                           (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(alim, 100))*sizeof(qfb));
                           alloc += FLINT_MIN(alim, 100);
                           for (k = num; k < alloc; k++)
                              qfb_init((*forms) + k);
                       }

                       fmpz_set_si((*forms)[num].a, a); /* record form */
                       fmpz_set_si((*forms)[num].b, b);
                       fmpz_set_si((*forms)[num++].c, c);
                   }
               }
           }
        }

        flint_free(s);
    }

    flint_free(fac);

    return num;
}

slong qfb_reduced_forms(qfb ** forms, slong d)
{
    slong a, b, k, c, p, blim, alloc, num, sqrt, i, prod, prime_i;
    mp_srcptr primes;
    const double * prime_inverses;
    mp_limb_t b2, exp, primes_cutoff = 0;
    n_factor_t * fac;
    mp_limb_t * s;

    if (d >= 0)
        flint_throw(FLINT_ERROR, "%s not implemented for positive discriminant\n", __func__);

    blim = n_sqrt(-d/3); /* maximum a value to check */

    (*forms) = NULL;
    alloc = 0;

    if ((-d & 3) == 2 || (-d & 3) == 1) /* ensure d is 0, 1 mod 4 */
        return 0;

    sqrt = n_sqrt(blim*blim - d);
    n_nth_prime_bounds(&primes_cutoff, &primes_cutoff, sqrt);
    if (primes_cutoff > FLINT_PRIMES_SMALL_CUTOFF*FLINT_PRIMES_SMALL_CUTOFF)
       return qfb_reduced_forms_large(forms, d);

    primes = n_primes_arr_readonly(FLINT_MAX(sqrt, 10000));
    prime_inverses = n_prime_inverses_arr_readonly(FLINT_MAX(sqrt, 10000));

    fac = flint_calloc(blim + 1, sizeof(n_factor_t));

    prime_i = 1;
    while ((p = primes[prime_i]) <= sqrt) /* sieve for factors of p^exp */
    {
        num = n_sqrtmod_primepow(&s, n_negmod((-d) % p, p), p, 1);

        for (i = 0; i < num; i++) /* sieve with each sqrt mod p */
        {
            mp_limb_t off = s[i];
            while (off <= blim)
            {
                b2 = (off*off - (mp_limb_t) d)/4;

                fac[off].p[fac[off].num] = p;
                fac[off].exp[fac[off].num] = n_remove2_precomp(&b2, p, prime_inverses[prime_i]);
                fac[off].num++;

                off += p;
            }
        }

        prime_i++;

        flint_free(s);
    }

    for (b = (d & 1); b <= blim; b += 2) /* write any remaining factors, including 2^exp */
    {
        b2 = ((mp_limb_t)(b*b - d))/4;

        exp = flint_ctz(b2); /* powers of 2 */
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

         b2 = ((mp_limb_t)(b*b - d))/4;

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
                     g = n_gcd(c, b);
                     g = n_gcd(a, g);
                 } else
                     g = n_gcd(c, a);

                 if (g == 1) /* primitive form */
                 {
                     if (num == alloc) /* realloc if necessary */
                     {
                         (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(blim, 100))*sizeof(qfb));
                         alloc += FLINT_MIN(blim, 100);
                         for (k = num; k < alloc; k++)
                            qfb_init((*forms) + k);
                    }

                     fmpz_set_si((*forms)[num].a, a); /* record form */
                     fmpz_set_si((*forms)[num].b, b);
                     fmpz_set_si((*forms)[num++].c, c);

                     if (b && a != c && b != a)
                     {
                         if (num == alloc) /* realloc if necessary */
                         {
                             (*forms) = flint_realloc(*forms, (alloc + FLINT_MIN(blim, 100))*sizeof(qfb));
                             alloc += FLINT_MIN(blim, 100);
                             for (k = num; k < alloc; k++)
                                qfb_init((*forms) + k);
                         }

                        fmpz_set_si((*forms)[num].a, a); /* record opposite form */
                        fmpz_set_si((*forms)[num].b, -b);
                        fmpz_set_si((*forms)[num++].c, c);
                     }
                 }
             }
         } while (pow_incr(pows, fac[b].exp, n));
    }

    flint_free(fac);

    return num;
}
