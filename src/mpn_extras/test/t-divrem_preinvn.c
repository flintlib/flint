/*
    Copyright (C) 2013 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_divrem_preinvn, state)
{
    int i, result;
    mpz_t a, d, q1, q2, r1, r2;
    gmp_randstate_t st;
    mp_ptr dinv;
    mp_size_t size, size2;
    flint_bitcnt_t norm;

    mpz_init(a);
    mpz_init(d);
    mpz_init(q1);
    mpz_init(r1);
    /* don't init r2, q2 */

    gmp_randinit_default(st);

    /* test flint_mpn_divrem_preinvn alias r and a */
    for (i = 0; i < 10000; i++)
    {
       size = n_randint(state, 200) + 1;
       size2 = n_randint(state, 200) + size;

       mpz_rrandomb(a, st, size2*FLINT_BITS);
       do {
          mpz_rrandomb(d, st, size*FLINT_BITS);
       } while (mpz_sgn(d) == 0);

       /* normalise */
       norm = flint_clz(d->_mp_d[d->_mp_size - 1]);
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

       q2->_mp_d[q2->_mp_size - 1] = flint_mpn_divrem_preinvn(q2->_mp_d, a->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv);

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
          fflush(stdout);
          flint_abort();
       }

       flint_free(dinv);
       flint_free(q2->_mp_d);
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
       norm = flint_clz(d->_mp_d[d->_mp_size - 1]);
       mpz_mul_2exp(d, d, norm);
       mpz_mul_2exp(a, a, norm);
       size2 = a->_mp_size;

       /* allocate space */
       q2->_mp_size = size2 - size + 1;
       q2->_mp_d = flint_malloc(q2->_mp_size*sizeof(mp_limb_t));

       r2->_mp_size = size2;
       r2->_mp_d = flint_malloc(r2->_mp_size*sizeof(mp_limb_t));

       /* reduce a mod d */
       mpz_fdiv_qr(q1, r1, a, d);

       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       q2->_mp_d[q2->_mp_size - 1] = flint_mpn_divrem_preinvn(q2->_mp_d, r2->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv);

       /* normalise */
       while (size && r2->_mp_d[size - 1] == 0) size--;
       r2->_mp_size = size;

       size2 = q2->_mp_size;
       while (size2 && q2->_mp_d[size2 - 1] == 0) size2--;
       q2->_mp_size = size2;

       result = (mpz_cmp(r1, r2) == 0 && mpz_cmp(q1, q2) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          gmp_printf("%Zd\n", r2);
          gmp_printf("%Zd\n", q1);
          gmp_printf("%Zd\n", q2);
          flint_printf("size = %wd\n", size);
          flint_printf("size2 = %wd\n", size2);
          fflush(stdout);
          flint_abort();
       }

       flint_free(dinv);
       flint_free(q2->_mp_d);
       flint_free(r2->_mp_d);
    }

    mpz_clear(a);
    mpz_clear(d);
    mpz_clear(q1);
    mpz_clear(r1);
    /* don't clear r2, q2 */
    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
