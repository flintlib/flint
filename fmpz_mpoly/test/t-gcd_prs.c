/*
    Copyright (C) 2017 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mpoly.h"
#include "ulong_extras.h"
#include "profiler.h"

int
main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);

    flint_printf("gcd_prs....");
    fflush(stdout);


/*
    {
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c;
        fmpz_mpoly_univar_t ax, bx, cx;
        timeit_t time;
        const char* vars[] = {"x","y","z","t","a","b","c","d","e","f","g","h"};

        fmpz_mpoly_ctx_init(ctx, 4, ORD_LEX);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);

        fmpz_mpoly_set_str_pretty(a, "2*x^2078*y^9*z^53*t^15 + x^2100*y^6*z^40 + x^2003*y^7*z^40 + x^2000*y^6*z^40 - 28*x^1080*y^9*z^44*t^15 - 14*x^1102*y^6*z^31 - 4*x^1079*y^6*z^33*t^16 + 156*x^1078*y^6*z^33*t^15 - 2*x^1101*y^3*z^20*t + 78*x^1100*y^3*z^20 - 14*x^1005*y^7*z^31 - 14*x^1002*y^6*z^31 - 2*x^1004*y^4*z^20*t + 78*x^1003*y^4*z^20 - 2*x^1001*y^3*z^20*t + 78*x^1000*y^3*z^20 + 98*x^82*y^9*z^35*t^15 + 49*x^104*y^6*z^22 + 28*x^81*y^6*z^24*t^16 - 1092*x^80*y^6*z^24*t^15 + 14*x^103*y^3*z^11*t - 546*x^102*y^3*z^11 + 2*x^80*y^3*z^13*t^17 - 156*x^79*y^3*z^13*t^16 + 3042*x^78*y^3*z^13*t^15 + x^102*t^2 - 78*x^101*t + 1521*x^100 + 49*x^7*y^7*z^22 + 49*x^4*y^6*z^22 + 14*x^6*y^4*z^11*t - 546*x^5*y^4*z^11 + 14*x^3*y^3*z^11*t - 546*x^2*y^3*z^11 + x^5*y*t^2 - 78*x^4*y*t + 1521*x^3*y + x^2*t^2 - 78*x*t + 1521", vars, ctx);
        fmpz_mpoly_set_str_pretty(b, "4*x^1156*y^9*z^46*t^30 + 4*x^1178*y^6*z^33*t^15 + x^1200*y^3*z^20 + 4*x^1081*y^7*z^33*t^15 + 4*x^1078*y^6*z^33*t^15 + 2*x^1103*y^4*z^20 + 2*x^1100*y^3*z^20 + x^1006*y^5*z^20 + 2*x^1003*y^4*z^20 + x^1000*y^3*z^20 - 28*x^158*y^9*z^37*t^30 - 28*x^180*y^6*z^24*t^15 - 4*x^157*y^6*z^26*t^31 + 156*x^156*y^6*z^26*t^30 - 7*x^202*y^3*z^11 - 4*x^179*y^3*z^13*t^16 + 156*x^178*y^3*z^13*t^15 - x^201*t + 39*x^200 - 28*x^83*y^7*z^24*t^15 - 28*x^80*y^6*z^24*t^15 - 14*x^105*y^4*z^11 - 14*x^102*y^3*z^11 - 4*x^82*y^4*z^13*t^16 + 156*x^81*y^4*z^13*t^15 - 4*x^79*y^3*z^13*t^16 + 156*x^78*y^3*z^13*t^15 - 2*x^104*y*t + 78*x^103*y - 2*x^101*t + 78*x^100 - 7*x^8*y^5*z^11 - 14*x^5*y^4*z^11 - 7*x^2*y^3*z^11 - x^7*y^2*t + 39*x^6*y^2 - 2*x^4*y*t + 78*x^3*y - x*t + 39", vars, ctx);

        printf("A: "); fmpz_mpoly_print_pretty(a, vars, ctx); printf("\n");
        printf("B: "); fmpz_mpoly_print_pretty(b, vars, ctx); printf("\n");

        timeit_start(time);
        fmpz_mpoly_gcd_prs(c, a, b, ctx);
        timeit_stop(time);
        printf("C: "); fmpz_mpoly_print_pretty(c, vars, ctx); printf("\n");
        flint_printf("time: %wd\n",time->wall);

        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
    }
*/

    for (i = 0; i < 100 * flint_test_multiplier(); i++)
    {
        int ok;
        fmpz_t ac, bc, d;
        fmpz_mpoly_ctx_t ctx;
        fmpz_mpoly_t a, b, c, g, t;
        ordering_t ord;
        slong nvars, len1, len2, len3, exp_bound1, exp_bound2, exp_bound3;
        slong coeff_bits;

        ord = mpoly_ordering_randtest(state);
        nvars = n_randint(state, 3) + 1;

        fmpz_mpoly_ctx_init(ctx, nvars, ord);

        fmpz_init(ac);
        fmpz_init(bc);
        fmpz_init(d);
        fmpz_mpoly_init(a, ctx);
        fmpz_mpoly_init(b, ctx);
        fmpz_mpoly_init(c, ctx);
        fmpz_mpoly_init(g, ctx);
        fmpz_mpoly_init(t, ctx);

        for (j = 0; j < 1; j++)
        {
            len1 = n_randint(state, 18);
            len2 = n_randint(state, 18);
            len3 = n_randint(state, 14);
            exp_bound1 = n_randint(state, 4) + 2;
            exp_bound2 = n_randint(state, 4) + 2;
            exp_bound3 = n_randint(state, 3) + 2;
            coeff_bits = n_randint(state, 3) + 1;
            fmpz_mpoly_randtest(a, state, len1, exp_bound1, coeff_bits, ctx);
            fmpz_mpoly_randtest(b, state, len2, exp_bound2, coeff_bits, ctx);
            fmpz_mpoly_randtest(c, state, len3, exp_bound3, coeff_bits, ctx);

            fmpz_mpoly_mul_johnson(a, a, c, ctx);
            fmpz_mpoly_mul_johnson(b, b, c, ctx);
            fmpz_mpoly_gcd_prs(g, a, b, ctx);

            if (fmpz_mpoly_is_zero(g, ctx))
            {
                if (!fmpz_mpoly_is_zero(a, ctx) || !fmpz_mpoly_is_zero(b, ctx))
                {
                    flint_printf("FAIL\ngcd is zero but both inputs are not\n"
                                                     "i: %wd  j: %wd\n", i, j);
                    flint_abort();
                }
                continue;
            }

            ok = fmpz_mpoly_divides_monagan_pearce(a, a, g, ctx);
            ok = ok && fmpz_mpoly_divides_monagan_pearce(b, b, g, ctx);
            if (!ok)
            {
                flint_printf("FAIL\ngcd doesn't divide both inputs\n"
                                                     "i: %wd  j: %wd\n", i, j);
                flint_abort();
            }

            if (!fmpz_mpoly_gcd_is_unit(a, b, ctx))
                flint_printf("FAIL\ncofactors are not relatively prime\n"
                                                     "i: %wd  j: %wd\n", i, j);

        }

        fmpz_clear(ac);
        fmpz_clear(bc);
        fmpz_clear(d);
        fmpz_mpoly_clear(a, ctx);
        fmpz_mpoly_clear(b, ctx);
        fmpz_mpoly_clear(c, ctx);
        fmpz_mpoly_clear(g, ctx);       
        fmpz_mpoly_clear(t, ctx);
    }

    FLINT_TEST_CLEANUP(state);

    printf("PASS\n");
    return 0;
}

