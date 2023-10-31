/*
    Copyright (C) 2011 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "mpn_extras.h"

TEST_FUNCTION_START(flint_mpn_divides, state)
{
    int i, result;
    mpz_t a, b, c, g, s;
    mp_ptr temp;
    gmp_randstate_t st;

    mpz_init(a);
    mpz_init(b);
    mpz_init(c);
    mpz_init(s);
    /* don't init g */
    gmp_randinit_default(st);

    /* check if b divides a*b */
    for (i = 0; i < 10000; i++)
    {
       do {
          mpz_urandomb(a, st, n_randint(state, 200));
       } while (mpz_sgn(a) == 0);
       do {
          mpz_urandomb(b, st, n_randint(state, 200));
       } while (mpz_sgn(b) == 0);
       do {
          mpz_urandomb(c, st, n_randint(state, 200));
       } while (mpz_sgn(c) == 0);

       mpz_mul(c, a, b);

       g->_mp_d = flint_malloc((c->_mp_size - b->_mp_size + 1)*sizeof(mp_limb_t));
       temp = flint_malloc(b->_mp_size * sizeof(mp_limb_t));

       result = flint_mpn_divides(g->_mp_d, c->_mp_d, c->_mp_size, b->_mp_d, b->_mp_size, temp);
       g->_mp_size = c->_mp_size - b->_mp_size + 1;
       g->_mp_size -= (g->_mp_d[g->_mp_size - 1] == 0);

       result &= (mpz_cmp(g, a) == 0);
       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", c);
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", b);
          fflush(stdout);
          flint_abort();
       }

       flint_free(g->_mp_d);
       flint_free(temp);
    }

    /* check b does not divide a*b + s for s < b */
    for (i = 0; i < 10000; i++)
    {
       do {
          mpz_urandomb(a, st, n_randint(state, 200));
       } while (mpz_sgn(a) == 0);
       do {
          mpz_urandomb(b, st, n_randint(state, 200) + 2);
       } while (mpz_sgn(b) == 0 || b->_mp_size == 1);
       do {
          mpz_urandomb(c, st, n_randint(state, 200));
       } while (mpz_sgn(c) == 0);
       do {
          mpz_urandomb(s, st, n_randint(state, b->_mp_size));
       } while (mpz_sgn(s) == 0);

       mpz_mul(c, a, b);
       mpz_add(c, c, s);

       g->_mp_d = flint_malloc((c->_mp_size - b->_mp_size + 1)*sizeof(mp_limb_t));
       temp = flint_malloc(b->_mp_size * sizeof(mp_limb_t));

       result = !flint_mpn_divides(g->_mp_d, c->_mp_d, c->_mp_size, b->_mp_d, b->_mp_size, temp);

       if (!result)
       {
          flint_printf("FAIL:\n");
          gmp_printf("%Zd\n", c);
          gmp_printf("%Zd\n", a);
          gmp_printf("%Zd\n", b);
          gmp_printf("%Zd\n", s);
          fflush(stdout);
          flint_abort();
       }

       flint_free(g->_mp_d);
       flint_free(temp);
    }

    mpz_clear(a);
    mpz_clear(b);
    mpz_clear(c);
    mpz_clear(s);
    /* don't clear g */
    gmp_randclear(st);

    TEST_FUNCTION_END(state);
}
