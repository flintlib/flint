/*
    Copyright (C) 2011-2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "partitions.h"

/* Defined in t-partitions_fmpz_ui.c and t-partitions_fmpz_ui_threaded.c */
#define testdata testdata_fmpz_ui_threaded
/* Values mod 10^9 */
static const ulong testdata[][2] =
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
  {500000000, 143535392},
  {1000000000, 685688339},
  {UWORD(4000000000), 525944299},
  {0, 0},
};

#define NUM 5000

TEST_FUNCTION_START(partitions_fmpz_ui_threaded, state)
{
    slong i;

    flint_set_num_threads(2);

    {
        fmpz_t p;
        fmpz * v;

        fmpz_init(p);
        v = _fmpz_vec_init(NUM);

        arith_number_of_partitions_vec(v, NUM);

        for (i = 0; i < NUM; i++)
        {
            partitions_fmpz_ui(p, i);
            if (!fmpz_equal(p, v + i))
            {
                flint_printf("FAIL:\n");
                flint_printf("p(%wd) does not agree with power series\n", i);
                flint_printf("Computed p(%wd): ", i); fmpz_print(p); flint_printf("\n");
                flint_printf("Expected: "); fmpz_print(v + i); flint_printf("\n");
                flint_abort();
            }
        }

        _fmpz_vec_clear(v, NUM);

        for (i = 0; testdata[i][0] != 0; i++)
        {
            partitions_fmpz_ui(p, testdata[i][0]);

            if (fmpz_fdiv_ui(p, 1000000000) != testdata[i][1])
            {
                flint_printf("FAIL:\n");
                flint_printf("p(%wd) does not agree with known value mod 10^9\n",
                    testdata[i][0]);
                flint_printf("Computed: %wu\n", fmpz_fdiv_ui(p, 1000000000));
                flint_printf("Expected: %wu\n", testdata[i][1]);
                flint_abort();
            }
        }

        fmpz_clear(p);
    }

    TEST_FUNCTION_END(state);
}
#undef testdata
