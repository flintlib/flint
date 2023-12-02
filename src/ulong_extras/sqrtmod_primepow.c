/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "flint.h"
#include "ulong_extras.h"

slong n_sqrtmod_2pow(mp_limb_t ** sqrt, mp_limb_t a, slong exp)
{
    mp_limb_t r = (a & 1);
    mp_limb_t * s;

    if (exp == 0) /* special case for sqrt of 0 mod 1 */
    {
        *sqrt = flint_malloc(sizeof(mp_limb_t));
        (*sqrt)[0] = 0;

        return 1;
    }

    if (exp == 1) /* special case mod 2 */
    {
        *sqrt = flint_malloc(sizeof(mp_limb_t));

        if (r) (*sqrt)[0] = 1;
        else (*sqrt)[0] = 0;

        return 1;
    }

    if (exp == 2) /* special case mod 4 */
    {
        r = (a & 3);

        if (r < 2) /* 0, 1 mod 4 */
        {
           *sqrt = flint_malloc(sizeof(mp_limb_t)*2);

           (*sqrt)[0] = r;
           (*sqrt)[1] = r + 2;

           return 2;
        } else /* 2, 3 mod 4 */
        {
           *sqrt = NULL;
           return 0;
        }
    }

    if (r) /* a is odd */
    {
        mp_limb_t roots[2];
        slong i, ex, pow;

        if ((a & 7) != 1) /* check square root exists */
        {
            *sqrt = NULL;
            return 0;
        }

        roots[0] = 1; /* one of each pair of square roots mod 8 */
        roots[1] = 3;

        pow = 8;
        for (ex = 3; ex < exp; ex++, pow *= 2) /* lift roots */
        {
            i = 0;

            r = roots[0];
            if (((r*r) & (2*pow - 1)) == (a & (2*pow - 1)))
               roots[i++] = r;

            r = pow - r;
            if (((r*r) & (2*pow - 1)) == (a & (2*pow - 1)))
            {
               roots[i++] = r;
               if (i == 2) continue;
            }

            r = roots[1];
            if (((r*r) & (2*pow - 1)) == (a & (2*pow - 1)))
            {
               roots[i++] = r;
               if (i == 2) continue;
            }

            r = pow - r;
            roots[i] = r;
        }

        *sqrt = flint_malloc(sizeof(mp_limb_t)*4);

        (*sqrt)[0] = roots[0]; /* write out both pairs of roots */
        (*sqrt)[1] = pow - roots[0];
        (*sqrt)[2] = roots[1];
        (*sqrt)[3] = pow - roots[1];

        return 4;
    } else /* a is even */
    {
        slong i, k, num, pow;

        for (k = 2; k <= exp; k++) /* find highest power of 2 dividing a */
        {
            if (a & ((UWORD(1)<<k) - 1))
                break;
        }
        k--;

        if (a == 0)
        {
           a = (UWORD(1)<<(exp - k/2));
           num = (UWORD(1)<<(k/2));
           s = flint_malloc(num*sizeof(mp_limb_t));
           for (i = 0; i < num; i++)
               s[i] = i*a;

           *sqrt = s;
           return num;
        }

        if (k & 1) /* not a square */
        {
            *sqrt = NULL;
            return 0;
        }

        pow = (UWORD(1)<<k);

        exp -= k;
        a /= pow;

        num = n_sqrtmod_2pow(&s, a, exp); /* divide through by 2^k and recurse */

        a = (UWORD(1)<<(k/2));
        r = a*(UWORD(1)<<exp);

        if (num == 0) /* check that roots were actually returned */
        {
            *sqrt = NULL;
            return 0;
        }

        for (i = 0; i < num; i++) /* multiply roots by 2^(k/2) */
            s[i] *= a;

        if (num == 1) /* one root */
        {
           s = flint_realloc(s, a*sizeof(mp_limb_t));

           for (i = 1; i < a; i++)
              s[i] = s[i - 1] + r;

        } else if (num == 2) /* two roots */
        {
           s = flint_realloc(s, 2*a*sizeof(mp_limb_t));

           for (i = 1; i < a; i++)
           {
              s[2*i] = s[2*i - 2] + r;
              s[2*i + 1] = s[2*i - 1] + r;
           }
        } else /* num == 4, i.e. four roots */
        {
           s = flint_realloc(s, 4*a*sizeof(mp_limb_t));

           for (i = 1; i < a; i++)
           {
              s[4*i] = s[4*i - 4] + r;
              s[4*i + 1] = s[4*i - 3] + r;
              s[4*i + 2] = s[4*i - 2] + r;
              s[4*i + 3] = s[4*i - 1] + r;
           }
        }

        *sqrt = s;

        return num*a;
    }
}

slong n_sqrtmod_primepow(mp_limb_t ** sqrt, mp_limb_t a, mp_limb_t p, slong exp)
{
    mp_limb_t r, ex, pow, k, a1, pinv, powinv;
    mp_limb_t * s;
    slong i, num;

    if (exp < 0)
    {
        flint_throw(FLINT_ERROR, "Exception (n_sqrtmod_primepow). exp must be non-negative.\n");
    }

    if (exp == 0) /* special case, sqrt of 0 mod 1 */
    {
        *sqrt = flint_malloc(sizeof(mp_limb_t));
        (*sqrt)[0] = 0;

        return 1;
    }

    if (p == 2) /* deal with p = 2 specially */
       return n_sqrtmod_2pow(sqrt, a, exp);

    if (exp == 1) /* special case, roots mod p */
    {
        r = n_sqrtmod(a, p);

        if (r == 0 && a != 0)
        {
            *sqrt = NULL;
            return 0;
        }

        *sqrt = flint_malloc(sizeof(mp_limb_t)*(1 + (r != 0)));
        (*sqrt)[0] = r;
        if (r) (*sqrt)[1] = p - r;

        return 1 + (r != 0);
    }

    pinv = n_preinvert_limb(p);
    a1 = n_mod2_preinv(a, p, pinv);
    r = n_sqrtmod(a1, p);

    if (r == 0 && a1 != 0)
    {
        *sqrt = NULL;
        return 0;
    }

    if (r) /* gcd(a, p) = 1, p is odd, lift r and p - r */
    {
        for (ex = 1, pow = p; ex < exp; ex++, pow *= p) /* lift root */
        {
            /* set k = ((r^2 - a) mod (p^(ex + 1))) / p^ex */
            powinv = n_preinvert_limb(pow*p);
            a1 = n_mulmod2_preinv(r, r, pow*p, powinv);
            k = (a < a1 ? a1 - a : a - a1);
            k = n_mod2_preinv(k, pow*p, powinv);
            k /= pow;
            if (a < a1)
                k = n_negmod(k, p);

            /* set k = k / 2r mod p */
            a1 = n_mulmod2_preinv(2, r, p, pinv);
            k = n_mulmod2_preinv(n_invmod(a1, p), k, p, pinv);

            /* set r = r + k*p^ex */
            r += k*pow;
        }

        *sqrt = flint_malloc(sizeof(mp_limb_t)*2);
        (*sqrt)[0] = r;
        (*sqrt)[1] = pow - r;

        return 2;
    } else /* special case, one root lifts to p roots */
    {
        for (k = 1, pow = p; k < exp; k++) /* find highest power of p dividing a */
        {
            mp_limb_t pow2 = pow * p;

            if (a % pow2 != 0)
                break;

            pow = pow2;
        }

        if (a == 0) /* special case, a == 0 */
        {
           a = n_pow(p, exp - k/2);
           num = n_pow(p, k/2);
           s = flint_malloc(num*sizeof(mp_limb_t));
           for (i = 0; i < num; i++)
               s[i] = i*a;

           *sqrt = s;
           return num;
        }

        if (k & 1) /* not a square */
        {
            *sqrt = NULL;
            return 0;
        }

        exp -= k;
        a /= pow;

        num = n_sqrtmod_primepow(&s, a, p, exp); /* divide through by p^k and recurse */

        if (num == 0)
        {
            *sqrt = NULL;
            return 0;
        }

        a = n_pow(p, k/2);
        r = a*n_pow(p, exp);

        s[0] *= a; /* multiply roots by p^(k/2) */
        s[1] *= a;

        s = flint_realloc(s, 2*a*sizeof(mp_limb_t));

        for (i = 1; i < a; i++)
        {
            s[2*i] = s[2*i - 2] + r;
            s[2*i + 1] = s[2*i - 1] + r;
        }

        *sqrt = s;

        return 2*a;
    }
}

