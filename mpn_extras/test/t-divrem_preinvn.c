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

    Copyright (C) 2013 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "longlong.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    mpz_t a, d, q1, q2, r1;
    gmp_randstate_t st;
    flint_rand_t state;
    mp_ptr dinv;
    mp_size_t size, size2;
    mp_bitcnt_t norm;
    
    flint_printf("divrem_preinvn....");
    fflush(stdout);

    mpz_init(a);
    mpz_init(d);
    mpz_init(q1);
    mpz_init(r1);
    /* don't init q2 */

    gmp_randinit_default(st);
    flint_randinit(state);

    /* test flint_mpn_divrem_n_preinvn */
    for (i = 0; i < 10000; i++)
    {
       size = n_randint(state, 200) + 1;
       
       do {
          mpz_rrandomb(a, st, 2*size*FLINT_BITS);
          do {
             mpz_rrandomb(d, st, size*FLINT_BITS);
          } while (mpz_sgn(d) == 0);
       
          /* normalise */
          count_leading_zeros(norm, d->_mp_d[d->_mp_size - 1]);
          mpz_mul_2exp(d, d, norm);
          mpz_mul_2exp(a, a, norm);
          if (a->_mp_size > 2*size)
             a->_mp_size = 2*size;
          if (mpn_cmp(a->_mp_d + size, d->_mp_d, size) >= 0)
             mpn_sub_n(a->_mp_d + size, a->_mp_d + size, d->_mp_d, size);
       
          size2 = 2*size;
          while (size2 && a->_mp_d[size2 - 1] == 0) size2--;
          if (size2 < size)
          {
             a->_mp_d[size - 1] = 1;
             size2 = size;
          }
          a->_mp_size = size2;
       } while (mpn_cmp(a->_mp_d + size2 - size, d->_mp_d, size) >= 0);

       /* reduce a mod d */
       mpz_fdiv_qr(q1, r1, a, d);
       
       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       q2->_mp_d = flint_malloc(size*sizeof(mp_limb_t));
       
       flint_mpn_divrem_n_preinvn(q2->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv); 

       /* normalise */
       size2 -= size;
       while (size && a->_mp_d[size - 1] == 0) size--;
       a->_mp_size = size;
       a->_mp_alloc = size;

       while (size2 && q2->_mp_d[size2 - 1] == 0) size2--;
       q2->_mp_size = size2;
       q2->_mp_alloc = size2;

       result = (mpz_cmp(r1, a) == 0 && mpz_cmp(q1, q2) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", q1);
          gmp_printf("%Zd\n", q2);
          gmp_printf("%Zd\n", r1);
          flint_printf("size = %wd\n", size);
          abort();
       }

       flint_free(q2->_mp_d);
       flint_free(dinv);
    }

    /* test flint_mpn_divrem_preinvn */
    for (i = 0; i < 10000; i++)
    {
       size = n_randint(state, 200) + 1;
       size2 = n_randint(state, 200) + size;
       
       mpz_rrandomb(a, st, size2*FLINT_BITS);
       do {
          mpz_rrandomb(d, st, size*FLINT_BITS);
       } while (mpz_sgn(d) == 0);
       
       /* normalise */
       count_leading_zeros(norm, d->_mp_d[d->_mp_size - 1]);
       mpz_mul_2exp(d, d, norm);
       mpz_mul_2exp(a, a, norm);
       size2 = a->_mp_size;

       /* allocate space */
       q2->_mp_size = size2 - size + 1;  
       q2->_mp_d = flint_malloc(q2->_mp_size*sizeof(mp_limb_t));
       
       /* reduce a mod d */
       mpz_fdiv_qr(q1, r1, a, d);
       
       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       q2->_mp_d[q2->_mp_size - 1] = flint_mpn_divrem_preinvn(q2->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv); 

       /* normalise */
       while (size && a->_mp_d[size - 1] == 0) size--;
       a->_mp_size = size;
       
       size2 = q2->_mp_size;
       while (size2 && q2->_mp_d[size2 - 1] == 0) size2--;
       q2->_mp_size = size2;
       
       result = (mpz_cmp(r1, a) == 0 && mpz_cmp(q1, q2) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          gmp_printf("%Zd\n", q1);
          gmp_printf("%Zd\n", q2);
          flint_printf("size = %wd\n", size);
          flint_printf("size2 = %wd\n", size2);
          abort();
       }

       flint_free(dinv);
       flint_free(q2->_mp_d);
    }

    mpz_clear(a);
    mpz_clear(d);
    mpz_clear(q1);
    mpz_clear(r1);
    /* don't init q2 */

    gmp_randclear(st);
    flint_randclear(state);

    flint_printf("PASS\n");
    return 0;
}
