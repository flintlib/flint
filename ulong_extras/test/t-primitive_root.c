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
    Copyright (C) 2015 Nitin Kumar
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

int main()
{
   int i, b, j, k;
   mp_limb_t p, n, root, phi;
   double ninv;
   n_factor_t factors;
   FLINT_TEST_INIT(state);

   flint_printf("n_primitive_root....");
   fflush(stdout);

   /* check for p ^ k, where p is and odd prime */

   for (i = 0; i < 100; i++)
   {
       p = n_randtest_prime(state, 1);

       if (p == 2) continue;

       b = FLINT_BIT_COUNT(p);

       for (j = 2; j <= FLINT_BITS / b; j++)
       {
           n = n_pow(p, j);
           phi = (n / p) * (p - 1);
           root = n_primitive_root(n);
           ninv = n_precompute_inverse(n);
           n_factor_init(&factors);
           n_factor(&factors, p - 1, 1);
           factors.p[factors.num] = p;
           factors.exp[factors.num] = 1;
           factors.num++;

           for (k = 0; k < factors.num; k++)
           {
               if (n_powmod_precomp(root, phi / factors.p[k], n, ninv) == 1)
               {
                   flint_printf("FAIL:\n");
                   flint_printf("%wu ** (%wu / %wu) == 1 mod %wu\n", root, phi, factors.p[k], n);
                   abort();
               }
           }
       }

   }

   /* check for 2 * (p ^ k), where p is an odd prime*/
   for (i = 0; i < 100; i++)
   {
       p = n_randtest_prime(state, 1);

       if (p == 2) continue;

       b = FLINT_BIT_COUNT(p);

       for (j = 1; j <= (FLINT_BITS - 2) / b; j++)
       {
           n = 2 * n_pow(p, j);
           phi = (n / p) * (p - 1);
           phi /= 2;
           root = n_primitive_root(n);
           ninv = n_precompute_inverse(n);
           n_factor_init(&factors);

           if (j > 1)
           {
               n_factor(&factors, p - 1, 1);
               factors.p[factors.num] = p;
               factors.exp[factors.num] = 1;
               factors.num++;
           }

           for (k = 0; k < factors.num; k++)
           {
               if (n_powmod_precomp(root, phi / factors.p[k], n, ninv) == 1)
               {
                   flint_printf("FAIL:\n");
                   flint_printf("%wu ** (%wu / %wu) == 1 mod %wu\n", root, phi, factors.p[k], n);
                   abort();
               }
           }
       }

   }

   FLINT_TEST_CLEANUP(state);

   flint_printf("PASS\n");
   return 0;
}
