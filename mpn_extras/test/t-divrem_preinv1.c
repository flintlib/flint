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
#include "mpn_extras.h"
#include "ulong_extras.h"

int main(void)
{
    int i, result;
    mpz_t a, a2, b, q, r, q2;
    gmp_randstate_t st;
    flint_rand_t state;
    mp_limb_t d1, d2, inv;
    len_t s1, s2;
    
    printf("divrem_preinv1....");
    fflush(stdout);

    mpz_init(a);
    mpz_init(a2);
    mpz_init(b);
    mpz_init(q);
    mpz_init(r);
    
    gmp_randinit_default(st);
    flint_randinit(state);

    for (i = 0; i < 10000; i++)
    {
       do {
          mpz_rrandomb(a, st, n_randint(state, 200));
          do {
             mpz_rrandomb(b, st, n_randint(state, 200));
          } while (mpz_sgn(b) == 0);

          s1 = a->_mp_size;
          s2 = b->_mp_size;
       } while (s1 < s2 || s2 < 2);
       
       mpz_set(a2, a);
       
       /* normalise b */
       b->_mp_d[b->_mp_size - 1] |= ((mp_limb_t) 1 << (GMP_LIMB_BITS - 1));

       d1 = b->_mp_d[b->_mp_size - 1];
       d2 = b->_mp_d[b->_mp_size - 2];

       mpz_fdiv_qr(q, r, a, b);

       inv = flint_mpn_preinv1(d1, d2);

       q2->_mp_d = flint_malloc((s1 - s2 + 1)*sizeof(mp_limb_t));
       
       q2->_mp_d[s1 - s2] = flint_mpn_divrem_preinv1(q2->_mp_d, a2->_mp_d, a2->_mp_size, b->_mp_d, b->_mp_size, inv); 

       /* normalise */
       s1 -= (s2 - 1);
       while (s1 && q2->_mp_d[s1 - 1] == 0) s1--;
       q2->_mp_size = s1;
       q2->_mp_alloc = s1;

       while (s2 && a2->_mp_d[s2 - 1] == 0) s2--;
       a2->_mp_size = s2;

       result = (mpz_cmp(q, q2) == 0 && mpz_cmp(a2, r) == 0);
       if (!result)
       {
          printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", b);
          gmp_printf("%Zd\n", q);
          gmp_printf("%Zd\n", r);
          gmp_printf("%Zd\n", q2);
          gmp_printf("%Zd\n", a2);
          abort();
       }

       flint_free(q2->_mp_d);
    }

    mpz_clear(a);
    mpz_clear(a2);
    mpz_clear(b);
    mpz_clear(q);
    mpz_clear(r);
    /* don't clear g */
    gmp_randclear(st);
    flint_randclear(state);

    printf("PASS\n");
    return 0;
}
