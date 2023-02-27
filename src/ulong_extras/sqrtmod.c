/*
    Copyright (C) 2009 William Hart
    Copyright (C) 2011 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"

mp_limb_t n_sqrtmod(mp_limb_t a, mp_limb_t p) 
{
    slong i, r, m, iter;
    mp_limb_t p1, k, b, g, bpow, gpow, res;
    mp_limb_t pinv;

    if (a <= 1)
    {
        return a;
    }

    /* just do a brute force search */
    if (p < 600)
    {
        mp_limb_t t, t2;

        if (p > 50 && n_jacobi_unsigned(a, p) == -1)
            return 0;

        t = t2 = 0;

        while (t < (p - 1) / 2)
        {
            t2 = n_addmod(t2, 2*t + 1, p);
            t++;
            
            if (t2 == a)
                return t;
        }

        return 0;
    }

    if (n_is_square(p)) /* modulus is a square */
       return 0;
   
    if ((p & 1) == 0) /* modulus is even */
       return 0;

    pinv = n_preinvert_limb(p);

    if (n_jacobi_unsigned(a, p) == -1)
        return 0;

    if ((p & UWORD(3)) == 3)
    {
        return n_powmod2_ui_preinv(a, (p + 1)/4, p, pinv); /* p == 2^B - 1 isn't prime */
    }

    if ((p & UWORD(7)) == 5)
    {
       b = n_powmod2_ui_preinv(a, (p + 3)/8, p, pinv); /* p == 2^B - 3 isn't prime */
       g = n_mulmod2_preinv(b, b, p, pinv);

       if (g == a)
          return b;

       g = n_powmod2_ui_preinv(2, (p - 1)/4, p, pinv);
       return n_mulmod2_preinv(g, b, p, pinv);
    }

    r = 0;
    p1 = p - 1;

    do {
        p1 >>= UWORD(1); 
        r++;
    } while ((p1 & UWORD(1)) == 0);

    b = n_powmod2_ui_preinv(a, p1, p, pinv);

    for (k = 3; ; k+=2) /* 2 is a quadratic residue mod p = 8k + 1 */
    {
        if (n_jacobi_unsigned(k, p) == -1) break;
    }

    g = n_powmod2_ui_preinv(k, p1, p, pinv);
    res = n_powmod2_ui_preinv(a, (p1 + 1) / 2, p, pinv);

    iter = r - 1; /* maximum number of iterations possible if p is prime */
    
    while (b != 1)
    {
        bpow = b;
        m = 0;
        do
        {
            bpow = n_mulmod2_preinv(bpow, bpow, p, pinv);
            m++;
        } while (m < r && bpow != 1);
        gpow = g;
        for (i = 1; i < r - m; i++)
        {
            gpow = n_mulmod2_preinv(gpow, gpow, p, pinv);
        }
        res = n_mulmod2_preinv(res, gpow, p, pinv);
        g = n_mulmod2_preinv(gpow, gpow, p, pinv);
        b = n_mulmod2_preinv(b, g, p, pinv);
        r = m;
        if (iter-- == 0)
            return 0; /* too many iterations, p is not prime, prevents hang */
    }

    return res;
}

