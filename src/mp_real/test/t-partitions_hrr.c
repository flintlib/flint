/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arb.h"
#include "arith.h"
#include "mp_real.h"

/* mp_real_partitions_hrr and _mp_real_partitions_hrr (every flag
   combination: the dfloat batches or not, the rarely taken paths, one
   thread or several) against the power series for n < 2500 and against
   values mod 10^9 (generated with Sage) from 10^5 to 10^8 */

/* values mod 10^9 generated with Sage (as in partitions/test) */
static const ulong testdata_mp_real_partitions[][2] =
{
  {100000, 421098519},
  {100001, 33940350},
  {100002, 579731933},
  {100003, 625213730},
  {100004, 539454200},
  {100005, 69672418},
  {100006, 865684292},
  {100007, 641916724},
  {100008, 36737908},
  {100009, 293498270},
  {100010, 177812057},
  {100011, 756857293},
  {100012, 950821113},
  {100013, 824882014},
  {100014, 533894560},
  {100015, 660734788},
  {100016, 835912257},
  {100017, 302982816},
  {100018, 468609888},
  {100019, 221646940},
  {1000000, 104673818},
  {1000001, 980212296},
  {1000002, 709795681},
  {1000003, 530913758},
  {1000004, 955452980},
  {1000005, 384388683},
  {1000006, 138665072},
  {1000007, 144832602},
  {1000008, 182646067},
  {1000009, 659145045},
  {1000010, 17911162},
  {1000011, 606326324},
  {1000012, 99495156},
  {1000013, 314860251},
  {1000014, 497563335},
  {1000015, 726842109},
  {1000016, 301469541},
  {1000017, 227491620},
  {1000018, 704160927},
  {1000019, 995311980},
  {10000000, 677288980},
  {10000001, 433805210},
  {10000002, 365406948},
  {10000003, 120899894},
  {10000004, 272822040},
  {10000005, 71938624},
  {10000006, 637670808},
  {10000007, 766947591},
  {10000008, 980210244},
  {10000009, 965734705},
  {10000010, 187411691},
  {10000011, 485652153},
  {10000012, 825498761},
  {10000013, 895802660},
  {10000014, 152775845},
  {10000015, 791493402},
  {10000016, 299640598},
  {10000017, 383615481},
  {10000018, 378922331},
  {10000019, 37059200},
  {100000000, 836637702},
  {100000001, 66421565},
  {100000002, 747849093},
  {100000003, 465329748},
  {100000004, 166747980},
  {0, 0},
};

static int
_mp_real_equal_fmpz(const mp_real_t x, const fmpz_t f)
{
    arb_t a;
    fmpz_t g;
    int r;
    arb_init(a);
    fmpz_init(g);
    mp_real_get_arb(a, x);
    r = arb_is_exact(a) && arb_get_unique_fmpz(g, a) && fmpz_equal(g, f)
        && (x->err == 0) && (x->size == 0 || x->d[0] != 0);
    arb_clear(a);
    fmpz_clear(g);
    return r;
}

#define NUM 2500

TEST_FUNCTION_START(mp_real_partitions_hrr, state)
{
    fmpz * v;
    fmpz_t n, f;
    mp_real_t x;
    slong i, iter;
    int flags;

    v = _fmpz_vec_init(NUM);
    fmpz_init(n);
    fmpz_init(f);
    mp_real_init(x);

    arith_number_of_partitions_vec(v, NUM);

    for (i = 0; i < NUM; i++)
    {
        mp_real_partitions_hrr(x, 0, i);
        if (!_mp_real_equal_fmpz(x, v + i))
            TEST_FUNCTION_FAIL("p(%wd)\n", i);

        flags = n_randint(state, 4);
        if (n_randint(state, 4) == 0)
            flint_set_num_threads(1 + n_randint(state, 4));
        fmpz_set_si(n, i);
        _mp_real_partitions_hrr(x, n, flags);
        flint_set_num_threads(1);
        if (!_mp_real_equal_fmpz(x, v + i))
            TEST_FUNCTION_FAIL("p(%wd) by the formula, flags %d\n", i, flags);
    }

    /* negative n */
    fmpz_set_si(n, -1 - (slong) n_randint(state, 1000));
    _mp_real_partitions_hrr(x, n, 0);
    if (x->size != 0)
        TEST_FUNCTION_FAIL("p(-n)\n");

    for (iter = 0; testdata_mp_real_partitions[iter][0] != 0; iter++)
    {
        ulong m = testdata_mp_real_partitions[iter][0];
        arb_t a;

        if (m > 10000000 && n_randint(state, 4) != 0)
            continue;

        flags = n_randint(state, 4);
        if (n_randint(state, 2))
            flint_set_num_threads(2 + n_randint(state, 3));
        if (n_randint(state, 2))
            mp_real_partitions_hrr(x, 0, m);
        else
        {
            fmpz_set_ui(n, m);
            _mp_real_partitions_hrr(x, n, flags);
        }
        flint_set_num_threads(1);

        arb_init(a);
        mp_real_get_arb(a, x);
        if (!arb_is_exact(a) || !arb_get_unique_fmpz(f, a)
            || fmpz_fdiv_ui(f, 1000000000) != testdata_mp_real_partitions[iter][1])
            TEST_FUNCTION_FAIL("p(%wu) mod 10^9, flags %d\n", m, flags);
        arb_clear(a);
    }

    _fmpz_vec_clear(v, NUM);
    fmpz_clear(n);
    fmpz_clear(f);
    mp_real_clear(x);

    TEST_FUNCTION_END(state);
}
