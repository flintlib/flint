/*
    Copyright (C) 2009 Thomas Boothby
    Copyright (C) 2009, 2017 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#define ulong ulongxx /* interferes with system includes */
#include <math.h>
#undef ulong
#include "flint.h"
#include "ulong_extras.h"

int n_is_perfect_power(ulong * root, ulong n)
{
    static unsigned char mod63[63] = {7,7,4,0,5,4,0,5,6,5,4,4,0,4,4,0,5,4,5,4,
                              4,0,5,4,0,5,4,6,7,4,0,4,4,0,4,6,7,5,4,0,4,4,0,5,
                              4,4,5,4,0,5,4,0,4,4,4,6,4,0,5,4,0,4,6};
    static unsigned char mod61[61] = {7,7,0,3,1,1,0,0,2,3,0,6,1,5,5,1,1,0,0,1,
                              3,4,1,2,2,1,0,3,2,4,0,0,4,2,3,0,1,2,2,1,4,3,1,0,
                              0,1,1,5,5,1,6,0,3,2,0,0,1,1,3,0,7};
    static unsigned char mod44[44] = {7,7,0,2,3,3,0,2,2,3,0,6,7,2,0,2,3,2,0,2,
                              3,6,0,6,2,3,0,2,2,2,0,2,6,7,0,2,3,3,0,2,2,2,0,6};
    static unsigned char mod31[31] = {7,7,3,0,3,5,4,1,3,1,1,0,0,0,1,2,3,0,1,1,
                              1,0,0,2,0,5,4,2,1,2,6};
    static unsigned char mod72[72] = {7,7,0,0,0,7,0,7,7,7,0,7,0,7,0,0,7,7,0,7,
                              0,0,0,7,0,7,0,7,0,7,0,7,7,0,0,7,0,7,0,0,7,7,0,7,
                              0,7,0,7,0,7,0,0,0,7,0,7,7,0,0,7,0,7,0,7,7,7,0,7,
                              0,0,0,7};
    static unsigned char mod49[49] = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,
                              0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,0,1};
    static unsigned char mod67[67] = {2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,2,2,0,0,0,0,0,0,2,2,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2};
    static unsigned char mod79[79] = {4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,4,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,0,0,4,4,0,0,0,0,0,0,0,0,0,0,0,
                              0,0,0,0,0,0,0,0,0,0,4};

    unsigned char t;
    ulong count, exp, r;
    
    /* check for powers 2, 3, 5 */
    t = mod31[n%31];
    t &= mod44[n%44];
    t &= mod61[n%61];
    t &= mod63[n%63];
    
    if (t & 1) 
    {
        ulong y = n_sqrtrem(&r, n);
        if (r == 0)
        {
           *root = y;
           return 2;
        }
    }

    if (t & 2) 
    {
        ulong y = n_cbrtrem(&r, n);
        if (r == 0)
        if (n == n_pow(y, 3))
        {
           *root = y;
           return 3;
        }
    }

    if (t & 4) 
    {
        ulong y = n_rootrem(&r, n, 5);
        if (r == 0)
        {
           *root = y;
           return 5;
        }
    }

    /* check for power 7, 11, 13 */
    t = mod49[n%49];
    t |= mod67[n%67];
    t |= mod79[n%79];
    t &= mod72[n%72];
    
    if (t & 1) 
    {
        ulong y = n_rootrem(&r, n, 7);
        if (r == 0)
        {
           *root = y;
           return 7;
        }
    }

    if (t & 2) 
    {
        ulong y = n_rootrem(&r, n, 11);
        if (r == 0)
        {
           *root = y;
           return 11;
        }
    }

    if (t & 13) 
    {
        ulong y = n_rootrem(&r, n, 13);
        if (r == 0)
        {
           *root = y;
           return 13;
        }
    }

    /* highest power of 2 */
    count_trailing_zeros(count, n);
    n >>= count;

    if (n == 1)
    {
       if (count == 1)
          return 0;
       *root = 2;
       return count;
    }

    /* check other powers (exp >= 17, root <= 13 and odd) */
    exp = 0;
    while ((n % 3) == 0)
    {
       n /= 3;
       exp += 1;
    }
    if (exp > 0)
    {
       if (n == 1 && exp > 1)
       {
          if (count == 0)
          {
             *root = 3;
             return exp;
          } else if (count == exp)
          {
             *root = 6;
             return exp;
          } else if (count == 2*exp)
          {
             *root = 12;
             return exp;
          }
       }
       return 0;
    }

#if FLINT64

    exp = 0;
    while ((n % 5) == 0)
    {
       n /= 5;
       exp += 1;
    }
    if (exp > 0)
    {
       if (n == 1 && exp > 1)
       {
          if (count == 0)
          {
             *root = 5;
             return exp;
          } else if (count == exp)
          {
             *root = 10;
             return exp;
          }
       }
       return 0;
    }

    if (count > 0)
       return 0;

    exp = 0;
    while ((n % 7) == 0)
    {
       n /= 7;
       exp += 1;
    }
    if (exp > 0)
    {
       if (n == 1 && exp > 1)
       {
          *root = 7;
          return exp;
       } 
       return 0;
    }

    exp = 0;
    while ((n % 11) == 0)
    {
       n /= 11;
       exp += 1;
    }
    if (exp > 0)
    {
       if (n == 1 && exp > 1)
       {
          *root = 11;
          return exp;
       } 
       return 0;
    }

    exp = 0;
    while ((n % 13) == 0)
    {
       n /= 13;
       exp += 1;
    }
    if (exp > 0)
    {
       if (n == 1 && exp > 1)
       {
          *root = 13;
          return exp;
       } 
       return 0;
    }

#endif

    return 0;
}

