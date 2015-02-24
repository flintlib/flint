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


#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#define ulong mp_limb_t
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"

int
n_rootrem_newton_iteration(mp_limb_t* base, mp_limb_t* remainder, mp_limb_t n, mp_limb_t root)
{
    mp_limb_t max_pow;
    double x, dx, a, b, small_float;

    if (root > 10)
        small_float = 0.001;
    else if (root >= 7)
        small_float = 0.00001;
    else 
        small_float = 0.000001;

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

    /* Newton iteration begins */

    x = n/root;                 

    dx = 0;
    a = root;
    b = root - 1;

    while (1)
    {
        dx = (n / pow(x, b));   // dx = n / x^n-1
        dx = dx - x;            // dx = ((n / x^n-1) - x)
        dx = dx / n;            // dx = ((n / x^n-1) - x) / n

        x = x + dx;

        if (absolute(dx) < absolute(x)*(small_float))
            break;
    }

    *remainder = x;
    *base = *remainder;
    *remainder = pow(*remainder, root);
    *remainder = n - *remainder;

    return 1;
}
