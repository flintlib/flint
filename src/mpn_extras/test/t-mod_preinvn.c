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

TEST_FUNCTION_START(flint_mpn_mod_preinvn, state)
{
    int i, result;
    mpz_t a, d, r1, r2;
    gmp_randstate_t st;
    mp_ptr dinv;
    mp_size_t size, size2;
    flint_bitcnt_t norm;

    mpz_init(a);
    mpz_init(d);
    mpz_init(r1);
    /* don't init r2 */

    gmp_randinit_default(st);

    /* test flint_mpn_mod_preinvn */
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

       /* make space for r */
       r2->_mp_size = size2;
       r2->_mp_d = flint_malloc(r2->_mp_size*sizeof(mp_limb_t));

       /* reduce a mod d */
       mpz_fdiv_r(r1, a, d);

       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       flint_mpn_mod_preinvn(r2->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv);

       /* normalise */
       while (size && r2->_mp_d[size - 1] == 0) size--;
       r2->_mp_size = size;

       result = (mpz_cmp(r1, r2) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          flint_printf("size = %wd\n", size);
          flint_printf("size2 = %wd\n", size2);
          fflush(stdout);
          flint_abort();
       }

       flint_free(dinv);
       flint_free(r2->_mp_d);
    }

    /* test flint_mpn_mod_preinvn alias r and a */
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

       /* reduce a mod d */
       mpz_fdiv_r(r1, a, d);

       dinv = flint_malloc(size*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size);

       flint_mpn_mod_preinvn(a->_mp_d, a->_mp_d, size2, d->_mp_d, size, dinv);

       /* normalise */
       while (size && a->_mp_d[size - 1] == 0) size--;
       a->_mp_size = size;

       result = (mpz_cmp(r1, a) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          flint_printf("size = %wd\n", size);
          flint_printf("size2 = %wd\n", size2);
          fflush(stdout);
          flint_abort();
       }

       flint_free(dinv);
    }

    mpz_clear(a);
    mpz_clear(d);
    mpz_clear(r1);
    /* don't clear r2 */

    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
