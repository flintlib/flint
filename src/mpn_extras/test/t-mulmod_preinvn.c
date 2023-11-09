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

TEST_FUNCTION_START(flint_mpn_mulmod_preinvn, state)
{
    int i, result;
    mpz_t a, b, d, r1, r2;
    gmp_randstate_t st;
    mp_ptr dinv;
    mp_size_t size, size_d;
    flint_bitcnt_t norm;

    mpz_init(a);
    mpz_init(b);
    mpz_init(d);
    mpz_init(r1);
    /* don't init r2 */

    gmp_randinit_default(st);

    for (i = 0; i < 10000; i++)
    {
       size = n_randint(state, 200) + 1;

       mpz_rrandomb(a, st, size);
       mpz_rrandomb(b, st, size);
       do {
          mpz_rrandomb(d, st, size);
       } while (mpz_sgn(d) == 0);

       size_d = d->_mp_size;
       /* reduce a, b mod d */
       mpz_fdiv_r(a, a, d);
       mpz_fdiv_r(b, b, d);

       mpz_mul(r1, a, b);
       mpz_fdiv_r(r1, r1, d);

       /* normalise */
       norm = flint_clz(d->_mp_d[size_d - 1]);
       mpz_mul_2exp(a, a, norm);
       mpz_mul_2exp(b, b, norm);
       mpz_mul_2exp(d, d, norm);

       dinv = flint_malloc(size_d*sizeof(mp_limb_t));
       flint_mpn_preinvn(dinv, d->_mp_d, size_d);

       r2->_mp_d = flint_malloc(size_d*sizeof(mp_limb_t));

       flint_mpn_mulmod_preinvn(r2->_mp_d, a->_mp_d, b->_mp_d, size_d, d->_mp_d, dinv, norm);

       /* normalise */
       r2->_mp_alloc = size_d;
       while (size_d && r2->_mp_d[size_d - 1] == 0) size_d--;
       r2->_mp_size = size_d;

       mpz_div_2exp(r2, r2, norm);
       mpz_div_2exp(a, a, norm);
       mpz_div_2exp(b, b, norm);
       mpz_div_2exp(d, d, norm);

       result = (mpz_cmp(r1, r2) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", b);
          gmp_printf("%Zd\n", d);
          gmp_printf("%Zd\n", r1);
          gmp_printf("%Zd\n", r2);
          flint_printf("size_d = %wd\n", size_d);
          fflush(stdout);
          flint_abort();
       }

       if (r2->_mp_alloc)
          flint_free(r2->_mp_d);
       flint_free(dinv);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(d);
    mpz_clear(r1);
    /* don't init r2 */

    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
