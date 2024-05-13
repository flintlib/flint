/*
    Copyright (C) 2008 Peter Shrimpton
    Copyright (C) 2009 William Hart
    Copyright (C) 2014, 2015 Dana Jacobsen
    Copyright (C) 2015 Kushagra Singh

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"

int n_is_prime(ulong n)
{
    /* flint's "BPSW" checked against Feitsma and Galway's database [1, 2]
       up to 2^64 by Dana Jacobsen.
       [1]  http://www.janfeitsma.nl/math/psp2/database
       [2]  http://www.cecm.sfu.ca/Pseudoprimes/index-2-to-64.html
    */

    if (n < 11) {
        if (n == 2 || n == 3 || n == 5 || n == 7)   return 1;
        else                                        return 0;
    }
    if (!(n%2) || !(n%3) || !(n%5) || !(n%7))       return 0;
    if (n <  121) /* 11*11 */                       return 1;
    if (!(n%11) || !(n%13) || !(n%17) || !(n%19) ||
        !(n%23) || !(n%29) || !(n%31) || !(n%37) ||
        !(n%41) || !(n%43) || !(n%47) || !(n%53))   return 0;
    if (n < 3481) /* 59*59 */                       return 1;
    if (n > 1000000 &&
        (!(n% 59) || !(n% 61) || !(n% 67) || !(n% 71) || !(n% 73) ||
         !(n% 79) || !(n% 83) || !(n% 89) || !(n% 97) || !(n%101) ||
         !(n%103) || !(n%107) || !(n%109) || !(n%113) || !(n%127) ||
         !(n%131) || !(n%137) || !(n%139) || !(n%149)))  return 0;

    return n_is_probabprime(n);
}

int
n_is_prime_pocklington(ulong n, ulong iterations)
{
    int pass;
    slong i;
    ulong j;
    ulong n1, cofactor, b, c, ninv, limit, F, Fsq, det, rootn, val, c1, c2, upper_limit;
    n_factor_t factors;
    c = 0;

#if FLINT64
    upper_limit = 2642246;                  /* 2642246^3 is approximately 2^64 */
#else
    upper_limit = 1626;                     /* 1626^3 is approximately 2^32 */
#endif

    if (n == 1)
        return 0;
    else if (n % 2 == 0)
        return (n == UWORD(2));

    rootn = n_sqrt(n);                      /* floor(sqrt(n)) */

    if (n == rootn*rootn)
        return 0;

    n1 = n - 1;
    n_factor_init(&factors);
    limit = (ulong) pow((double)n1, 1.0/3);

    val = n_pow(limit, 3);

    while (val < n1 && limit < upper_limit)    /* ensuring that limit >= n1^(1/3) */
    {
        limit++;
        val = n_pow(limit, 3);
    }


    cofactor = n_factor_partial(&factors, n1, limit, 1);

    if (cofactor != 1)                      /* check that cofactor is coprime to factors found */
    {
        for (i = 0; i < factors.num; i++)
        {
            if (factors.p[i] > FLINT_FACTOR_TRIAL_PRIMES_PRIME)
            {
                while (cofactor >= factors.p[i] && (cofactor % factors.p[i]) == 0)
                {
                    factors.exp[i]++;
                    cofactor /= factors.p[i];
                }
            }
        }
    }
    F = n1/cofactor;                        /* n1 = F*cofactor */
    Fsq = F*F;

    if (F <= rootn)                         /* cube root method applicable only if n^1/3 <= F < n^1/2 */
    {
        c2 = n1/(Fsq);                      /* expressing n as c2*F^2 + c1*F + 1  */
        c1 = (n1 - c2*Fsq )/F;
        det = c1*c1 - 4*c2;
        if (n_is_square(det))               /* BSL's test for (n^1/3 <= F < n^1/2) */
            return 0;
    }
    ninv = n_preinvert_limb(n);
    c = 1;
    for (i = factors.num - 1; i >= 0; i--)
    {
        ulong exp = n1 / factors.p[i];
        pass = 0;

        for (j = 2; j < iterations && pass == 0; j++)
        {
            b = n_powmod2_preinv(j, exp, n, ninv);
            if (n_powmod2_ui_preinv(b, factors.p[i], n, ninv) != UWORD(1))
                return 0;

            b = n_submod(b, UWORD(1), n);
            if (b != UWORD(0))
            {
                c = n_mulmod2_preinv(c, b, n, ninv);
                pass = 1;
            }
            else if (c == 0)
                return 0;
        }
        if (j == iterations)
            return -1;
    }
    return (n_gcd(n, c) == UWORD(1));
}

ulong flint_pseudosquares[] = {17, 73, 241, 1009, 2641, 8089, 18001,
          53881, 87481, 117049, 515761, 1083289, 3206641, 3818929, 9257329,
          22000801, 48473881, 48473881, 175244281, 427733329, 427733329,
          898716289u, 2805544681u, 2805544681u, 2805544681u
#ifndef FLINT64
          };
#else
          , 10310263441u, 23616331489u, 85157610409u, 85157610409u,
          196265095009u, 196265095009u, 2871842842801u, 2871842842801u,
          2871842842801u, 26250887023729u, 26250887023729u, 112434732901969u,
          112434732901969u, 112434732901969u, 178936222537081u,
          178936222537081u, 696161110209049u, 696161110209049u,
          2854909648103881u, 6450045516630769u, 6450045516630769u,
          11641399247947921u, 11641399247947921u, 190621428905186449u,
          196640248121928601u, 712624335095093521u, 1773855791877850321u };
#endif

#if FLINT64
#define FLINT_NUM_PSEUDOSQUARES 52
#else
#define FLINT_NUM_PSEUDOSQUARES 25
#endif

int n_is_prime_pseudosquare(ulong n)
{
    unsigned int i, j, m1;
    ulong p, B, NB, exp, mod8;
    const ulong * primes;
    const double * inverses;

    if (n < UWORD(2))
        return 0;
    else if ((n & UWORD(1)) == UWORD(0))
        return (n == UWORD(2));

    primes = n_primes_arr_readonly(FLINT_PSEUDOSQUARES_CUTOFF+1);
    inverses = n_prime_inverses_arr_readonly(FLINT_PSEUDOSQUARES_CUTOFF+1);

    for (i = 0; i < FLINT_PSEUDOSQUARES_CUTOFF; i++)
    {
        double ppre;
        p = primes[i];
        if (p*p > n)
            return 1;
        ppre = inverses[i];
        if (!n_mod2_precomp(n, p, ppre))
            return 0;
    }

    B  = primes[FLINT_PSEUDOSQUARES_CUTOFF];
    NB = (n - 1)/B + 1;
    m1 = 0;

    for (i = 0; i < FLINT_NUM_PSEUDOSQUARES; i++)
        if (flint_pseudosquares[i] > NB)
            break;

    exp = (n - 1)/2;

    for (j = 0; j <= i; j++)
    {
        ulong mod = n_powmod2(primes[j], exp, n);
        if ((mod != UWORD(1)) && (mod != n - 1))
            return 0;
        else if (mod == n - 1)
            m1 = 1;
    }

    mod8 = n % 8;

    if ((mod8 == 3) || (mod8 == 7))
        return 1;
    else if (mod8 == 5)
    {
        ulong mod = n_powmod2(UWORD(2), exp, n);
        if (mod == n - 1)
            return 1;
        else
            flint_throw(FLINT_ERROR, "Whoah, %wu is a probable prime, but not prime, please report!!\n", n);
    }
    else
    {
        if (m1) return 1;
        for (j = i + 1; j < FLINT_NUM_PSEUDOSQUARES + 1; j++)
        {
            ulong mod = n_powmod2(primes[j], exp, n);
            if (mod == n - 1)
                return 1;
            else if (mod != 1)
                flint_throw(FLINT_ERROR, "Whoah, %wu is a probable prime, but not prime, please report!!\n", n);
        }
        flint_throw(FLINT_ERROR, "Whoah, %wu is a probable prime, but not prime, please report!!\n", n);
    }
}
