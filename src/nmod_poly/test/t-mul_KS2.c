/*
    Copyright (C) 2010 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_poly.h"

#pragma GCC diagnostic ignored "-Woverlength-strings"

TEST_FUNCTION_START(nmod_poly_mul_KS2, state)
{
    int i, result;

    /* Check aliasing of a and b */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_KS2(a, b, c);
        nmod_poly_mul_KS2(b, b, c);

        result = (nmod_poly_equal(a, b));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(b), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Check aliasing of a and c */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_KS2(a, b, c);
        nmod_poly_mul_KS2(c, b, c);

        result = (nmod_poly_equal(a, c));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a), flint_printf("\n\n");
            nmod_poly_print(c), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Compare with mul_classical */
    for (i = 0; i < 200 * flint_test_multiplier(); i++)
    {
        nmod_poly_t a1, a2, b, c;
        mp_limb_t n = n_randtest_not_zero(state);

        nmod_poly_init(a1, n);
        nmod_poly_init(a2, n);
        nmod_poly_init(b, n);
        nmod_poly_init(c, n);
        nmod_poly_randtest(b, state, n_randint(state, 50));
        nmod_poly_randtest(c, state, n_randint(state, 50));

        nmod_poly_mul_classical(a1, b, c);
        nmod_poly_mul_KS2(a2, b, c);

        result = (nmod_poly_equal(a1, a2));
        if (!result)
        {
            flint_printf("FAIL:\n");
            nmod_poly_print(a1), flint_printf("\n\n");
            nmod_poly_print(a2), flint_printf("\n\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a1);
        nmod_poly_clear(a2);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
    }

    /* Test bug */
#if FLINT64
    {
        nmod_poly_t a, b, c, d;
        mp_limb_t mod = UWORD(2289285083314003039);

        nmod_poly_init(a, mod);
        nmod_poly_init(b, mod);
        nmod_poly_init(c, mod);
        nmod_poly_init(d, mod);

        nmod_poly_set_str(a, "33 2289285083314003039  1904046980750977501 45214318121847844 54776012950656710 1873599154826904295 870135259339956207 1979537227904177235 1077518167660740425 1914467488553877071 590032441505981152 615648453231634975 1569207985886566133 787136232386763586 263398180027397052 1072218041043012468 1506848477788670239 1400920451698857943 1647489479045838018 916805681536849287 418919486780459023 1905019227786610376 521214770020411309 1686949157600795332 1694792566051380615 859359964104912916 379633023194464188 1707896212900599917 2116886930226258819 1784312697572836983 1657809908840472396 187671865737908075 24295635882237532 1236324514297047805 1");

        nmod_poly_set_str(b, "33 2289285083314003039  1248232897348310716 2007649684622883320 1745423566053694824 377131386390909586 510000530211074748 2252357745328136405 121845726568848648 1776552397423576952 1087867787512095029 1542258909416820365 333466255708033066 1369137418799234940 1019757368632255267 1708944472854199555 424720899313432606 1832061349660947209 1680133579237359807 1011480105427615896 469487889066997171 250011917413341143 1528554661920601095 1534412092525363608 1513772824610977843 1546061195142719878 2202390411223099645 1876123527101934779 777483746883723994 1298099847659446396 1845415258502877349 82593368833294076 7103430832858844 1414221645839662896 1");

        nmod_poly_mul_classical(c, a, b);
        nmod_poly_mul_KS2(d, a, b);

        if (!nmod_poly_equal(c, d))
        {
            flint_printf("FAIL: NMOD_RED2 bug!\n");
            fflush(stdout);
            flint_abort();
        }

        nmod_poly_clear(a);
        nmod_poly_clear(b);
        nmod_poly_clear(c);
        nmod_poly_clear(d);
    }
#endif

    TEST_FUNCTION_END(state);
}
