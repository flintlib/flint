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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"

int
mp_limb_rootrem(mp_limb_t* base, mp_limb_t* remainder, mp_limb_t n, mp_limb_t root)
{
    mp_limb_t a, u, t, s, q, max_pow;

    if (n<1 || root<1)
        return 0;

    if (root == 1)
    {
        *base = n;
        *remainder = 0;
        return 1;
    }
        
    if (n <= root)
    {
        *base = 1;
        *remainder = n - 1;
        return 1;
    }

    if (root == 2)
    {
        *base = n_sqrt(n);      /* floor of sqrt(n) */
        a = n_pow(*base, 2);

        if (a>n)                /* doc says n_sqrt may return value rounded up */
        {
            *base = *base - 1;
            a = n_pow(*base, 2);
        }

        *remainder = n - a;
        return 0;
    }

    max_pow = n_clog(n, 2);      /* 2^(max_pow) > n */

    if (n == n_pow(2, max_pow)  && (root == max_pow))
    {
        *base = 2;
        *remainder = 0;
        return 1;
    }

    if (root >= max_pow)         /* if root > max_pow, base has to be 1 (we know its not 2) */
    {

        *base = 1;
        *remainder = n - 1;
        return 1;
    }


    u = n_sqrt(n) + 1;           /* u can be any value greater than n^1/k, we know k > 2 */
    s = n;

    do{
        s = u;
        q = n_clog(n, s);
        t = (root-1)*s;
        if (q > (root-1))
            t+= (n/n_pow(s, root-1));
        u = t/root;
    }while (u<s);

    *base = s;

    a = n_pow(s, root);
    *remainder = n - a;

    return 1;
}
