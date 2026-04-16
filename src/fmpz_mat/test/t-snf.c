/*
    Copyright (C) 2026 Edgar Costa

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_mat.h"

/*
    Helper: compute SNF of A, check is_in_snf, rank, idempotency, aliasing.
*/
static void
_check_snf_rand(const fmpz_mat_t A, slong expected_rank)
{
    fmpz_mat_t S, T;
    slong i, m, n, snf_rank;

    m = fmpz_mat_nrows(A);
    n = fmpz_mat_ncols(A);

    fmpz_mat_init(S, m, n);
    fmpz_mat_init(T, m, n);

    fmpz_mat_snf(S, A);

    if (!fmpz_mat_is_in_snf(S))
    {
        flint_printf("FAIL: matrix not in snf!\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    snf_rank = 0;
    for (i = 0; i < FLINT_MIN(m, n); i++)
        if (!fmpz_is_zero(fmpz_mat_entry(S, i, i)))
            snf_rank++;

    if (snf_rank != expected_rank)
    {
        flint_printf("FAIL: snf rank %wd != expected %wd\n",
            snf_rank, expected_rank);
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    /* Idempotency */
    fmpz_mat_snf(T, S);

    if (!fmpz_mat_equal(S, T))
    {
        flint_printf("FAIL: snf not idempotent!\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fmpz_mat_print_pretty(T); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    /* Aliasing */
    fmpz_mat_set(T, A);
    fmpz_mat_snf(T, T);

    if (!fmpz_mat_equal(T, S))
    {
        flint_printf("FAIL: aliasing test failed!\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fmpz_mat_print_pretty(T); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mat_clear(T);
    fmpz_mat_clear(S);
}

/*
    Helper: build m x n matrix from slong array, compute SNF, and check
    that the diagonal matches the expected values (verified with Magma).
*/
static void
_test_snf_case(slong m, slong n, const slong * data,
        const slong * expected_snf)
{
    fmpz_mat_t A, S, S2;
    slong i, j, d;

    d = FLINT_MIN(m, n);

    fmpz_mat_init(A, m, n);
    fmpz_mat_init(S, m, n);
    fmpz_mat_init(S2, m, n);

    for (i = 0; i < m; i++)
        for (j = 0; j < n; j++)
            fmpz_set_si(fmpz_mat_entry(A, i, j), data[i * n + j]);

    fmpz_mat_snf(S, A);

    if (!fmpz_mat_is_in_snf(S))
    {
        flint_printf("FAIL:\n");
        flint_printf("matrix not in snf!\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    for (i = 0; i < d; i++)
    {
        if (!fmpz_equal_si(fmpz_mat_entry(S, i, i), expected_snf[i]))
        {
            flint_printf("FAIL:\n");
            flint_printf("diagonal entry %wd is %{fmpz}, expected %wd\n",
                    i, fmpz_mat_entry(S, i, i), expected_snf[i]);
            fmpz_mat_print_pretty(A); flint_printf("\n\n");
            fmpz_mat_print_pretty(S); flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }
    }

    fmpz_mat_snf(S2, S);
    if (!fmpz_mat_equal(S, S2))
    {
        flint_printf("FAIL:\n");
        flint_printf("snf of a matrix in snf should be the same!\n");
        fmpz_mat_print_pretty(A); flint_printf("\n\n");
        fmpz_mat_print_pretty(S); flint_printf("\n\n");
        fmpz_mat_print_pretty(S2); flint_printf("\n\n");
        fflush(stdout);
        flint_abort();
    }

    fmpz_mat_clear(S2);
    fmpz_mat_clear(S);
    fmpz_mat_clear(A);
}

TEST_FUNCTION_START(fmpz_mat_snf, state)
{
    slong iter;

    /* Hardcoded test cases verified with Magma */

    /* Issue #2592: 27x25 non-square matrix that previously hung */
    {
        slong data[] = {4,-8,9,-9,3,0,9,-2,-10,-8,-3,8,8,-10,-2,8,-9,-5,5,6,10,4,-2,-5,8,3,10,5,-8,5,1,3,0,0,-7,-5,0,3,5,-1,2,7,-9,4,-8,0,-2,0,-7,2,6,-10,7,4,3,-9,-4,6,1,9,5,10,4,-9,-4,-2,7,-6,-1,4,5,-7,-10,10,9,-3,-5,-1,7,-10,7,3,-8,-3,-7,4,-7,10,-6,5,-1,6,-2,3,5,5,-3,4,7,-6,2,-4,9,6,-6,-8,-2,3,0,6,-2,-10,-1,-1,8,8,5,-6,4,7,5,1,0,7,7,2,4,0,-4,-3,8,2,-3,3,-9,0,5,2,2,10,-6,5,-9,-6,6,8,0,-7,4,-7,6,4,-10,-6,3,10,-6,-8,5,-2,0,9,2,10,-8,0,7,2,0,10,5,7,-9,9,-8,-3,10,-1,-3,-8,3,-7,10,-7,4,-5,-1,-10,-9,0,-9,-1,1,1,3,-6,-3,6,3,8,-5,-5,-5,-8,9,2,9,-3,5,8,-6,-3,4,10,-2,4,-2,-10,4,-1,7,-5,-8,4,1,8,-1,10,3,-2,4,-1,-4,2,5,-7,-3,2,8,1,8,-1,-1,-10,2,-2,-10,8,-9,9,5,-1,-3,9,1,-3,10,-4,9,-2,-6,10,-7,10,10,-9,-1,4,-9,8,1,-6,-8,-1,0,3,-5,-4,-6,7,1,6,6,-2,-5,-2,5,-1,0,-7,4,-8,-6,-3,2,7,1,-8,2,-10,-2,7,-7,4,1,-2,8,2,10,1,-7,-3,5,-10,9,7,0,9,-3,10,-8,10,4,-1,10,3,-7,-6,-9,-9,-1,5,-7,-7,-3,7,-6,2,4,1,7,3,8,-6,3,10,-7,5,9,3,-2,-9,1,-4,4,4,-3,1,-7,1,7,10,1,-9,2,-2,-4,-7,4,-8,-4,10,10,9,-10,-9,0,-3,-6,8,-4,-8,7,-4,8,-4,-3,0,-6,9,-10,-2,-6,-6,7,-2,-5,-7,-10,-6,-10,1,-3,8,0,-10,-5,-2,-9,-6,3,6,-7,-8,5,4,1,6,8,-7,4,6,-3,9,-9,6,-1,4,10,-10,-9,5,2,3,-7,5,4,-8,-8,0,9,-6,-8,-6,-2,9,10,8,7,0,2,9,6,-1,4,6,9,3,-7,-7,10,10,7,-4,3,4,-3,3,0,4,2,3,-7,0,3,0,-2,1,-6,5,-8,-8,-8,-8,3,-7,1,-6,7,-9,8,7,7,0,-7,3,1,3,-9,-1,9,-1,1,-7,8,6,-4,-6,5,-3,-7,1,7,1,-7,-2,8,-3,3,7,9,9,10,7,-10,9,-2,-10,-5,-2,-1,0,1,-10,-5,-6,8,2,-8,-6,10,-10,-8,6,-4,2,3,4,0,-5,1,-1,0,8,9,-8,-9,-6,-5,9,-9,-8,-2,4,3,5,9,4,3,-2,-4,6,-7,1,3,-7,-1,8,5,6,-1,-9,-3,2,9,-9,-10,-4,-1,-4,-6,-2,-1,0,-7,-10,5,3,-5,-6,2,7,-3,6,7,1,-8,2,-9,3,-10,4,-8,0,8,3,8,2,10,3,-1,-7,2,-10,0,-5,9,4,1,-8,3,-7,-3,3,8,2,6,-8,2,-1,0,-3,0,-5,-8,6,10,-7,6,6,-4,1,1,10,-6,-3,-7,-6,-2,-4,-5,9,-6,10,-8,-5,10,5,4,8,8,4,8,10,10};
        slong snf[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
        _test_snf_case(27, 25, data, snf);
    }

    /* Square singular matrices (det = 0) */
    {
        slong data[] = {90,-65,170,83,-59,60,-85,-64,-152,-26,76,86,-142,-119,-73,79,139,72,102,-25,-64,161,105,-155,-5,0,93,27,-12,-27,-234,13,-9,-36,-94,-70,-1,112,1,67,80,109,20,87,12,101,-68,-73,1,25,3,-48,41,-9,-20,-69,7,-95,-9,10,-1,28,-129,-41,-60,-128,120,-194,-137,-42,31,139,112,59,-66,-190,181,-128,227,65,-140,16,-128,-48,47,81,16,-146,95,71,173,-107,-21,116,-121,-7,-133,-53,90,-2,57,-5,12,-52,55,-41,55,-85,79,-10,16,46,106,-54,68,-10,-24,-19,36,-107,129};
        slong snf[] = {1,1,1,1,1,1,1,1,2,0,0};
        _test_snf_case(11, 11, data, snf);
    }
    {
        slong data[] = {-46,-61,-154,168,22,146,-1,37,-191,-21,59,-272,-15,-145,-68,-11,150,25,54,-38,-111,-10,72,144,-96,-64,142,118,-89,48,-8,48,41,1,-31,91,5,-16,18,-166,-243,-17,-124,191,-8,-87,28,133,-22,223,4,-52,147,-159,3,80,-18,-8,183,-234,48,-200,-28,-107,146,59,-15,-97,-29,-39,-152,61,-105,21,-154,233,178,-80,111,-38,32,180,-23,47,-76,-7,-20,27,42,48,125,-39,72,11,93,-144,-152,52,-61,31,76,-4,-119,-24,108,76,-142,-4,-18,41,59,-15,-24,118,-175,226,73,13,113,-113,35,81,126,60,18,87,-1,101,5,24,113,60,38,-8,-189,-5,66,-20,-165,-24,-63,-190,81,127,-2,95,131,76,-231,-126,-62,-87,-82,-76,-173,105,261,-42,-39,-60,69,8,151,-47,-58,-45,123,-18,1,-44,-68,175,-30,108,-19,73,-41,-50,-61,-18,68,27,100,-82,-88,-151,86,-25,-37,-81,-131,147,58,129,67,75,-166,-41,25,13,-71,6,-54,8,4,-132,-94,87,-38,29,-27,5,131,51,-27,130,-78,15,134,-92,-140,12,-3,181,-108,-78,88,151,-101,81,51,-1,42,57,-21,143,-7,33,86,-237,-6,68,-65,-110,-80,-139,30,-131,224,123,-14,144,-90,160,297,-18,-207,53,116,-43,98,70,88,-89,-50,-133,19,-7,-67,-112,-128,81,165,35,27,-11,-198,-58,-16,109,37,-67,-115,-52,48,-76,-32,-49,-123,-42,-168,-32,105,-6,-24,165,-69,-59,200,-21,122,53,151,-103,-20,-9,20,-27,55,-72,-3,155,-36,134,-37,9,-40,-19,-124,65,-28,98,136,-81,-24,-176,-44,72,122,39,-31,189,162,41,20,-30,-175,-23,-53,28,-76,-165,49,38,68,108,108,-37,-42,-38,-104,-78,-46,-281,-33,-142,192,-81,-86,109,158,-47,90,-103,54,100,75,-113,27,-46,-60,145,-133,90,-110,111,85,-72,-61,16,74,-111,-150,-32,56,112,56,-170,-37,0,-78,-68,-156,-7,34,91,88,-27,133,-22};
        slong snf[] = {1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
        _test_snf_case(20, 20, data, snf);
    }
    {
        slong data[] = {14,-72,24,-10,0,2,-82,-36,18,74,-34,-28,60,50,-6,46,-23,70,-45,9,-26,-7,66,35,2,-77,54,46,-41,-71,-5,-31,25,-110,45,-15,10,5,-120,-55,20,115,-60,-50,85,85,-5,65,-15,118,-21,17,22,1,146,59,-46,-117,38,30,-113,-63,19,-87,34,-56,72,-6,64,14,-30,-28,-34,70,-78,-68,4,94,22,2,-6,-10,-15,-2,-22,-4,-23,-5,19,6,13,12,23,-12,-10,18,-20,-16,-48,-4,-64,-12,-52,-8,52,4,44,40,56,-44,-28,44,-8,30,-15,4,-6,-2,31,15,-3,-32,19,16,-21,-26,0,-16,14,32,36,6,56,10,66,16,-50,-22,-30,-28,-64,26,26,-50,17,-80,30,-11,4,3,-89,-40,17,83,-41,-34,64,59,-5,49,7,-10,15,-1,14,3,-4,-5,-8,13,-16,-14,-1,19,5,-1,0,26,3,4,14,2,37,13,-17,-24,1,0,-31,-6,8,-24,19,10,45,3,58,11,42,5,-46,1,-42,-38,-47,43,25,-37,25,-110,45,-15,10,5,-120,-55,20,115,-60,-50,85,85,-5,65,6,-94,3,-14,-34,-4,-125,-47,49,90,-17,-12,101,36,-22,78,-2,92,6,14,44,6,128,46,-56,-86,8,4,-106,-26,26,-82};
        slong snf[] = {1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        _test_snf_case(16, 16, data, snf);
    }
    {
        slong data[] = {-6,122,4,-1,40,-92,2,-66,-71,-22,-61,-120,40,23,-3,12,221,-141,41,-106,-28,48,20,-96,-27,47,-4,-44,77,-21,128,30,51,31,-39,146,-17,64,69,-70,3,-64,-89,-32,-159,155,-14,17,40,148,-130,-27,-146,-25,23,-69,45,-178,23,-80,-23,54,11,64,107,60,33,3,4,-26,-16,-25,6,19,-13,1,-27,-48,-92,74,153,144,-34,60,-166,143,-52,87,68,-4,9,-63,-44,-27,1,59,-32,4,-52,-84,-16,135,55,116,158,-109,-55,-58,-98,41,-84,-116,-32,166,80,-12,92,-132,210,20,146};
        slong snf[] = {1,1,1,1,0,0,0,0,0,0,0};
        _test_snf_case(11, 11, data, snf);
    }

    /* Square nonsingular matrices with nontrivial invariant factors */
    {
        slong data[] = {8051,3,350,729,-1060,-2,2199,-1866,177,-732,-41534,-2,0,-3776,3806,6,-11326,11298,-12,3776,6330,9,1057,567,-1560,0,1737,-738,516,-576,-2831,5,-477,-255,735,0,-785,288,-255,256,-2110,-3,-352,-189,550,0,-579,216,-172,192,-11880,0,0,-1080,1080,2,-3240,3240,-6,1080,-14684,-10,-460,-1330,1752,4,-4004,3588,-222,1336,24861,30,3171,2235,-5244,0,6810,-3786,1539,-2262,-10098,8,-356,-916,1216,4,-2768,2448,-214,916,-37620,0,0,-3420,3420,6,-10260,10260,-18,3420};
        slong snf[] = {1,1,1,1,1,2,6,30,60,180};
        _test_snf_case(10, 10, data, snf);
    }
    {
        slong data[] = {901,0,0,298,0,0,-900,0,0,300,-3,-298,0,1,0,0,0,0,0,0,0,0,0,0,4494,-45,7,1526,180,0,-4305,0,-67,1500,7114,-1526,-21,44894,12,-22461,9,2,48,0,-12,-4,-54067,-43,63,-134685,-36,67383,-26,-6,-143,0,36,12,162199,129,-6,-3,0,-2,0,1,6,0,0,-2,1,0,0,-9,0,0,0,0,3,0,0,0,-6,0,8088,0,0,2700,540,-36,-7560,12,-144,2700,-1140,-2700,-2700,0,0,-900,-180,0,2520,0,60,-900,360,900,3600,-44976,0,23700,171,0,-3429,0,-60,1200,67158,-1200,3,6,-3,-12,4,0,-2,0,3,0,-2999,12,0,-14992,0,7500,-3,0,-3,0,0,0,22506,0};
        slong snf[] = {1,1,1,1,1,1,3,12,60,300,1500,7500};
        _test_snf_case(12, 12, data, snf);
    }

    /* Non-square matrices with nontrivial invariant factors */
    {
        slong data[] = {6,0,0,0,0,0,0,18,0,0,6,0,0,0,0,0,12,0,-432,0,0,0,0,0,0,0,0,0,0,0,0,0,48,0,0,48,0,0,0,0,0,0,0,0,0,0,-24,0,1008,0,0,0,0,0,5184,0,0,-31104,0,0,0,0,0,0,15984,3456,0,0,1728,2592,0,0,0,7776,0,0,0,-48,-5184,-46656,-14304,2592,0,-7776,-7776,0,-1296,0,-23328,0,0,36,0,-3888,0,0,1296,0,0,0,10368,-15552,0,0,0,-2592,0,0,2592,1296,-5184,-1296,-6480,-2592,15552,-2592,0,-3888,0,0,0,0,0,0,15552,5184,0,0,2592,2592,0,0,0,7776,0,0,0,0,0,0,0,0,0,0,5184,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5184,0,0,0,0,0,-12,0,432,1296,-5184,0,0,-2592,0,-10368,15552,0,0,0,0,0,0,0,-1296,5184,0,0,2592,0,31104,-46656,15552,0,0};
        slong snf[] = {6,12,48,144,432,1296,1296,1296,2592,5184,5184,15552,15552};
        _test_snf_case(13, 15, data, snf);
    }
    {
        slong data[] = {25,-6350,75,-300,3200,0,0,0,0,-900,0,0,0,-600,0,0,0,10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,-780,40,-600,5200,2400,0,0,115200,-1800,0,0,0,-1650,0,0,-10,3180,-40,100,-1600,0,0,0,0,300,0,0,0,150,0,0,0,-800,0,0,400,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-300,4800,2400,0,0,57600,-900,0,0,0,-900,0,0,-10,-20,10,200,0,0,3600,0,0,-3000,0,0,-7200,600,0,10800,15,30,-15,-300,0,0,0,14400,0,-900,-43200,0,0,-900,0,0,25,-6350,75,-600,5600,1200,0,0,57600,-1800,0,0,0,-1500,0,0};
        slong snf[] = {5,10,50,100,400,1200,3600,14400,57600};
        _test_snf_case(9, 16, data, snf);
    }
    {
        slong data[] = {78,-48,120,0,0,0,-42,48,0,12,0,6,0,0,0,0,0,0,0,0,24,-24,60,0,-108,-972,24,24,0,216,-2760,132,2952,1308,72,648,504,1512,0,-144,-24,24,3180,0,36,324,-24,1596,0,-72,4128,-270,600,-1944,-216,-2484,-732,240,0,432,-1416,132,-300,648,252,2484,204,-120,0,-504,48,-48,-6360,0,-108,-972,48,-3192,0,216};
        slong snf[] = {6,6,12,12,36,108,324,1620};
        _test_snf_case(8, 10, data, snf);
    }
    {
        slong data[] = {6,-24,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0,0,-270000,0,0,0,0,0,-12,66,12,0,24,0,-288000,72,90000,24,0,0,0,24,0,-18,72,-120,60,0,0,-360,54000,120,0,180,0,0,0,-162000,0,0,0,0,300,0,0,0,0,810000,0,0,0,0,-900,6,-6,0,0,0,1500,6000,1500,0,0,0,0,0,0,0,0,0,0,0,0,0,6000,0,0,0,0,0,0,0,0,0,0,0,0,900,0,0,18000,0,2430000,0,0,0,0,-56700,6,-6,240,-120,0,4500,-269280,4500,89760,0,-360,0,0,0,0,12,-12,480,-240,0,9000,-538560,9000,179520,270000,-720,0,0,0,0};
        slong snf[] = {6,6,12,60,300,1500,6000,18000,90000,270000};
        _test_snf_case(10, 15, data, snf);
    }

    /* Edge cases */
    {
        /* 0x0 matrix */
        fmpz_mat_t S;
        fmpz_mat_init(S, 0, 0);
        fmpz_mat_snf(S, S);
        fmpz_mat_clear(S);
    }
    {
        /* 0x5 matrix */
        fmpz_mat_t A, S;
        fmpz_mat_init(A, 0, 5);
        fmpz_mat_init(S, 0, 5);
        fmpz_mat_snf(S, A);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        /* 3x0 matrix */
        fmpz_mat_t A, S;
        fmpz_mat_init(A, 3, 0);
        fmpz_mat_init(S, 3, 0);
        fmpz_mat_snf(S, A);
        fmpz_mat_clear(S);
        fmpz_mat_clear(A);
    }
    {
        slong data[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
        slong snf[] = {0,0,0,0,0};
        _test_snf_case(5, 8, data, snf);
    }
    {
        slong data[] = {-7};
        slong snf[] = {7};
        _test_snf_case(1, 1, data, snf);
    }
    {
        slong data[] = {6,-4,10,0,8,-2,14,0,-12,20};
        slong snf[] = {2};
        _test_snf_case(1, 10, data, snf);
    }
    {
        slong data[] = {12,-18,6,0,-24,30,0,42};
        slong snf[] = {6};
        _test_snf_case(8, 1, data, snf);
    }

    /* Randomized tests */
    for (iter = 0; iter < 2000 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong m, n, b, d, r;

        m = n_randint(state, 30);
        n = n_randint(state, 30);
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);

        b = 1 + n_randint(state, 10) * n_randint(state, 10);
        d = n_randint(state, 2 * m * n + 1);
        fmpz_mat_randrank(A, state, r, b);

        if (n_randint(state, 2))
            fmpz_mat_randops(A, state, d);

        _check_snf_rand(A, r);
        fmpz_mat_clear(A);
    }

    /* Non-square matrices of moderate size */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong m, n, b, d, r;

        m = 10 + n_randint(state, 40);
        n = 10 + n_randint(state, 40);
        if (m == n)
            m++;
        r = n_randint(state, FLINT_MIN(m, n) + 1);

        fmpz_mat_init(A, m, n);

        b = 1 + n_randint(state, 10);
        d = n_randint(state, 2 * m * n + 1);
        fmpz_mat_randrank(A, state, r, b);
        fmpz_mat_randops(A, state, d);

        _check_snf_rand(A, r);
        fmpz_mat_clear(A);
    }

    /* Square nonsingular matrices above cutoff */
    for (iter = 0; iter < 20 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong n, b, d;

        n = 10 + n_randint(state, 20);

        fmpz_mat_init(A, n, n);

        b = 1 + n_randint(state, 10);
        d = n_randint(state, 2 * n * n + 1);
        fmpz_mat_randrank(A, state, n, b);
        fmpz_mat_randops(A, state, d);

        _check_snf_rand(A, n);
        fmpz_mat_clear(A);
    }

    /* Square singular matrices (zero det) */
    for (iter = 0; iter < 100 * flint_test_multiplier(); iter++)
    {
        fmpz_mat_t A;
        slong n, b, d, r;

        n = 5 + n_randint(state, 25);
        r = n_randint(state, n);

        fmpz_mat_init(A, n, n);

        b = 1 + n_randint(state, 10);
        d = n_randint(state, 2 * n * n + 1);
        fmpz_mat_randrank(A, state, r, b);
        fmpz_mat_randops(A, state, d);

        _check_snf_rand(A, r);
        fmpz_mat_clear(A);
    }

    /*
        Cross-check all implementations against each other on small
        square nonsingular matrices:
        - fmpz_mat_snf (dispatcher)
        - fmpz_mat_snf_kannan_bachem
        - fmpz_mat_snf_iliopoulos
        - fmpz_mat_snf_transform (diagonal of S)
        - fmpz_mat_elementary_divisors

        Keep max size small since kannan_bachem can hang on larger
        matrices (see #2592).
    */
    {
        const char * method_names[] = {
            "kannan_bachem", "iliopoulos",
            "snf_transform", "elementary_divisors"
        };
        slong num_methods = 4;

        for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
        {
            fmpz_mat_t A, S_ref, S_other, U, V;
            fmpz_t det;
            fmpz * ed;
            slong n, b, d, i, method, rank;

            n = n_randint(state, 9);

            fmpz_mat_init(A, n, n);
            fmpz_mat_init(S_ref, n, n);
            fmpz_mat_init(S_other, n, n);
            fmpz_mat_init(U, n, n);
            fmpz_mat_init(V, n, n);
            fmpz_init(det);
            ed = _fmpz_vec_init(n);

            b = 1 + n_randint(state, 10) * n_randint(state, 10);
            d = n_randint(state, 2 * n * n + 1);
            fmpz_mat_randrank(A, state, n, b);

            if (n_randint(state, 2))
                fmpz_mat_randops(A, state, d);

            /* Reference: dispatcher */
            fmpz_mat_snf(S_ref, A);

            fmpz_mat_det(det, A);
            fmpz_abs(det, det);

            for (method = 0; method < num_methods; method++)
            {
                fmpz_mat_zero(S_other);

                switch (method)
                {
                    case 0:
                        fmpz_mat_snf_kannan_bachem(S_other, A);
                        break;
                    case 1:
                        fmpz_mat_snf_iliopoulos(S_other, A, det);
                        break;
                    case 2:
                        fmpz_mat_snf_transform(S_other, U, V, A);
                        break;
                    case 3:
                        rank = fmpz_mat_elementary_divisors(ed, A);
                        if (rank != n)
                        {
                            flint_printf("FAIL: %s rank %wd != %wd"
                                "\n", method_names[method], rank, n);
                            fflush(stdout);
                            flint_abort();
                        }
                        for (i = 0; i < n; i++)
                            fmpz_set(fmpz_mat_entry(S_other, i, i),
                                &ed[i]);
                        break;
                }

                if (!fmpz_mat_equal(S_ref, S_other))
                {
                    flint_printf("FAIL: %s disagrees with "
                        "dispatcher!\n", method_names[method]);
                    flint_printf("n=%wd b=%wd\n", n, b);
                    flint_printf("A:          ");
                    fmpz_mat_print_pretty(A);
                    flint_printf("\n");
                    flint_printf("dispatcher: ");
                    fmpz_mat_print_pretty(S_ref);
                    flint_printf("\n");
                    flint_printf("%s: ", method_names[method]);
                    fmpz_mat_print_pretty(S_other);
                    flint_printf("\n");
                    fflush(stdout);
                    flint_abort();
                }
            }

            _fmpz_vec_clear(ed, n);
            fmpz_clear(det);
            fmpz_mat_clear(V);
            fmpz_mat_clear(U);
            fmpz_mat_clear(S_other);
            fmpz_mat_clear(S_ref);
            fmpz_mat_clear(A);
        }
    }

    TEST_FUNCTION_END(state);
}
