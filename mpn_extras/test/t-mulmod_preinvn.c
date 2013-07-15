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
    mpz_t a, b, d, r1, r2;
    gmp_randstate_t st;
    flint_rand_t state;
    mp_ptr dinv;
    mp_size_t size;
    mp_bitcnt_t norm;
    
    printf("mulmod_preinvn....");
    fflush(stdout);

    mpz_init(a);
    mpz_init(b);
    mpz_init(d);
    mpz_init(r1);
    /* don't init r2 */

    gmp_randinit_default(st);
    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
       size = n_randint(state, 200) + 2;
       
       mpz_rrandomb(a, st, size*FLINT_BITS);
       mpz_rrandomb(b, st, size*FLINT_BITS);
       do {
          mpz_rrandomb(d, st, size*FLINT_BITS);
       } while (mpz_sgn(d) == 0);
       
       /* reduce a, b mod d */
       mpz_fdiv_r(a, a, d);
       mpz_fdiv_r(b, b, d);

       mpz_mul(r1, a, b);
       mpz_fdiv_r(r1, r1, d);
       
       /* normalise */
       count_leading_zeros(norm, d->_mp_d[d->_mp_size - 1]);
       mpz_mul_2exp(a, a, norm);
       mpz_mul_2exp(b, b, norm);
       mpz_mul_2exp(d, d, norm);

       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       r2->_mp_d = flint_malloc(size*sizeof(mp_limb_t));
       
       flint_mpn_mulmod_preinvn(r2->_mp_d, a->_mp_d, b->_mp_d, size, d->_mp_d, dinv, norm); 

       /* normalise */
       while (size && r2->_mp_d[size - 1] == 0) size--;
       r2->_mp_size = size;
       r2->_mp_alloc = size;

       result = (mpz_cmp(r1, r2) == 0);
       if (!result)
       {
          printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", b);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          gmp_printf("%Zd\n", r2);
          printf("size = %ld\n", size);
          abort();
       }

       flint_free(r2->_mp_d);
       flint_free(dinv);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(d);
    mpz_clear(r1);
    /* don't init r2 */

    gmp_randclear(st);
    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
