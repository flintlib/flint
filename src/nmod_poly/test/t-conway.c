/*
    Copyright (C) 2024 Albin Ahlbäck

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <string.h>
#include "test_helpers.h"
#include "nmod_poly.h"

typedef struct
{
    ulong prime;
    slong deg;
    ulong * coeffs;
}
conway_polynomial_t;

typedef struct
{
    ulong prime;
    slong degree;
}
prime_degree_t;

/* As the source code splits it into 7 different cases, we need to check every
   case against some variation against CPimport.txt by Frank Lübeck. We do this
   by checking every available degree of some primes for each case. We also
   check that some non-primes and unavailable degrees return a value indicating
   an error. */

TEST_FUNCTION_START(_nmod_poly_conway, state)
{
    ulong op[500];
    slong ix;

    /* As the source code for primes < 260 does not depend on specific cases, we
       here simply check that a couple entries agree with the data. */
    ulong cp_2_1[] = {1,1};
    ulong cp_2_42[] = {1,1,1,0,0,1,1,0,0,1,0,1,1,0,0,0,0,0,1,0,1,0,0,0,1,1,1,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1};
    ulong cp_3_263[] = {1,1,2,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1};
    ulong cp_11_1[] = {9,1};
    ulong cp_257_12[] = {3,20,148,249,173,215,225,13,2,0,0,0,1};

    /* For 260 < prime < 300, we have the cases deg = 1, 2, ..., 12 */
    ulong cp_263_1[] = {258,1};
    ulong cp_283_2[] = {3,282,1};
    ulong cp_293_3[] = {291,2,0,1};
    ulong cp_293_4[] = {2,166,3,0,1};
    ulong cp_277_5[] = {272,1,0,0,0,1};
    ulong cp_269_6[] = {2,206,101,120,1,0,1};
    ulong cp_263_7[] = {258,1,0,0,0,0,0,1};
    ulong cp_271_9[] = {265,186,266,10,0,0,0,0,0,1};
    ulong cp_271_10[] = {6,126,74,256,10,133,1,0,0,0,1};
    ulong cp_271_11[] = {265,10,0,0,0,0,0,0,0,0,0,1};
    ulong cp_281_12[] = {3,191,28,58,116,103,68,202,0,0,0,0,1};

    /* For 300 <= prime < 1000, we have the cases deg = 1, 2, ..., 9 */
    ulong cp_311_1[] = {294,1};
    ulong cp_311_2[] = {17,310,1};
    ulong cp_307_3[] = {302,7,0,1};
    ulong cp_409_4[] = {21,407,12,0,1};
    ulong cp_419_5[] = {417,4,0,0,0,1};
    ulong cp_443_6[] = {2,41,218,298,1,0,1};
    ulong cp_523_7[] = {521,13,0,0,0,0,0,1};
    ulong cp_787_8[] = {2,715,26,612,5,0,0,0,1};
    ulong cp_919_9[] = {912,623,410,7,0,0,0,0,0,1};

    /* For 1000 <= prime < 3371, we have the cases deg = 1, ..., 7, 9. However,
       there are seven primes that do not have degree 7 and 9. We cover them
       later on in the test. */
    ulong cp_1009_1[] = {998,1};
    ulong cp_3359_2[] = {11,3357,1};
    ulong cp_2311_3[] = {2308,3,0,1};
    ulong cp_2333_4[] = {2,1591,7,0,1};
    ulong cp_2293_5[] = {2291,2,0,0,0,1};
    ulong cp_2287_6[] = {19,233,630,1812,0,0,1};
    ulong cp_2273_7[] = {2270,20,0,0,0,0,0,1};
    ulong cp_3361_9[] = {3339,3147,1704,5,0,0,0,0,0,1};

    /* For 3371 <= prime < 11000, we have the cases deg = 1, 2, ..., 6 */
    ulong cp_3371_1[] = {3369,1};
    ulong cp_5011_2[] = {2,5010,1};
    ulong cp_6679_3[] = {6672,1,0,1};
    ulong cp_8191_4[] = {17,8183,3,0,1};
    ulong cp_9679_5[] = {9676,8,0,0,0,1};
    ulong cp_10993_6[] = {7,4015,6869,8311,0,0,1};

    /* For 11000 <= prime < 65536, we have the cases deg = 1, 2, ..., 4 */
    ulong cp_11003_1[] = {11001,1};
    ulong cp_36109_2[] = {2,36108,1};
    ulong cp_38867_3[] = {38865,2,0,1};
    ulong cp_65521_4[] = {17,42121,20,0,1};

    /* For 65536 <= prime < 109988, we have the case deg = 4 */
    ulong cp_65537_4[] = {3,48035,7,0,1};
    ulong cp_65647_4[] = {3,54892,2,0,1};
    ulong cp_65651_4[] = {2,58463,2,0,1};
    ulong cp_65657_4[] = {3,54810,0,0,1};
    ulong cp_65677_4[] = {2,52138,3,0,1};
    ulong cp_77029_4[] = {2,58146,0,0,1};
    ulong cp_77137_4[] = {13,71265,4,0,1};
    ulong cp_88771_4[] = {11,68103,5,0,1};
    ulong cp_100069_4[] = {10,65117,2,0,1};
    ulong cp_109943_4[] = {5,109939,3,0,1};
    ulong cp_109961_4[] = {3,109949,0,0,1};
    ulong cp_109987_4[] = {3,100525,3,0,1};

    /* Expected outputs */
    conway_polynomial_t cp[] =
    {
        {2, 1, cp_2_1},
        {2, 42, cp_2_42},
        {3, 263, cp_3_263},
        {11, 1, cp_11_1},
        {257, 12, cp_257_12},

        {263, 1, cp_263_1},
        {283, 2, cp_283_2},
        {293, 3, cp_293_3},
        {293, 4, cp_293_4},
        {277, 5, cp_277_5},
        {269, 6, cp_269_6},
        {263, 7, cp_263_7},
        {271, 9, cp_271_9},
        {271, 10, cp_271_10},
        {271, 11, cp_271_11},
        {281, 12, cp_281_12},

        {311, 1, cp_311_1},
        {311, 2, cp_311_2},
        {307, 3, cp_307_3},
        {409, 4, cp_409_4},
        {419, 5, cp_419_5},
        {443, 6, cp_443_6},
        {523, 7, cp_523_7},
        {787, 8, cp_787_8},
        {919, 9, cp_919_9},

        {1009, 1, cp_1009_1},
        {3359, 2, cp_3359_2},
        {2311, 3, cp_2311_3},
        {2333, 4, cp_2333_4},
        {2293, 5, cp_2293_5},
        {2287, 6, cp_2287_6},
        {2273, 7, cp_2273_7},
        {3361, 9, cp_3361_9},

        {3371, 1, cp_3371_1},
        {5011, 2, cp_5011_2},
        {6679, 3, cp_6679_3},
        {8191, 4, cp_8191_4},
        {9679, 5, cp_9679_5},
        {10993, 6, cp_10993_6},

        {11003, 1, cp_11003_1},
        {36109, 2, cp_36109_2},
        {38867, 3, cp_38867_3},
        {65521, 4, cp_65521_4},

        {65537, 4, cp_65537_4},
        {65647, 4, cp_65647_4},
        {65651, 4, cp_65651_4},
        {65657, 4, cp_65657_4},
        {65677, 4, cp_65677_4},
        {77029, 4, cp_77029_4},
        {77137, 4, cp_77137_4},
        {88771, 4, cp_88771_4},
        {100069, 4, cp_100069_4},
        {109943, 4, cp_109943_4},
        {109961, 4, cp_109961_4},
        {109987, 4, cp_109987_4},
    };

    /* Expected errors */
    prime_degree_t error_list[] =
    {
        {0, 1},
        {1, 1},
        {2, 0},
        {2, -1},
        {2, 408},
        {2, 410},
        {2, 490},
        {6, 4},
        {259, 4},
        {263, -1},
        {263, 0},
        {263, 13},
        {279, 4},
        {297, 4},
        {307, 0},
        {307, 10},
        {999, 0},
        {999, 10},
        {999, 10},
        {1013, 8},
        {1013, 10},
        {1031, 8},
        {1031, 10},
        {2689, 7},
        {2689, 9},
        {2797, 7},
        {2797, 9},
        {2833, 7},
        {2833, 9},
        {3019, 7},
        {3019, 9},
        {3163, 7},
        {3163, 9},
        {3209, 7},
        {3209, 9},
        {3331, 7},
        {3331, 9},
        {3367, 4},
        {3371, 0},
        {3371, 7},
        {10993, 7},
        {10997, 4},
        {11003, 0},
        {11003, 5},
        {11117, 5},
        {11121, 4},
        {65521, 5},
        {65521, 0},
        {65527, 4},
        {65541, 4},
        {109981, 4},
        {109987, 0},
        {109989, 4},
        {UWORD(1782381726837121), 4},
    };

    for (ix = 0; ix < sizeof(cp) / sizeof(conway_polynomial_t); ix++)
    {
        ulong prime;
        slong deg;
        ulong * coeffs;
        int result;

        prime = cp[ix].prime;
        deg = cp[ix].deg;
        coeffs = cp[ix].coeffs;

        result = _nmod_poly_conway(op, prime, deg);

        result = result && !memcmp(coeffs, op, sizeof(ulong) * (deg + 1));

        if (!result)
            flint_throw(FLINT_TEST_FAIL,
                    "prime = %wu\n"
                    "degree = %wd\n"
                    "Got:      %{ulong*}\n"
                    "Expected: %{ulong*}\n",
                    prime, deg, op, deg + 1, coeffs, deg + 1);
    }

    for (ix = 0; ix < sizeof(error_list) / sizeof(prime_degree_t); ix++)
    {
        ulong prime;
        slong deg;
        int result;

        prime = error_list[ix].prime;
        deg = error_list[ix].degree;

        result = _nmod_poly_conway(op, prime, deg);

        if (result)
            flint_throw(FLINT_TEST_FAIL,
                    "Exected return value 0 for prime = %wu and degree %wd.\n"
                    "Got return value %d.\n",
                    prime, deg, result);
    }

    TEST_FUNCTION_END(state);
}

TEST_FUNCTION_START(_nmod_poly_conway_rand, state)
{
    slong ix;
    int result;
    ulong op[500];

    for (ix = 0; ix < 10000 * flint_test_multiplier(); ix++)
    {
        ulong prime;
        slong degree;

        /* Test all primes and degrees */
        prime = _nmod_poly_conway_rand(&degree, state, 0);

        result = _nmod_poly_conway(op, prime, degree);

        if (!result)
            flint_throw(FLINT_TEST_FAIL,
                    "Prime %wu and degree %wd not present in database.\n",
                    prime, degree);
    }

    TEST_FUNCTION_END(state);
}
