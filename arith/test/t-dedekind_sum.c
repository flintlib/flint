/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "arith.h"
#include "ulong_extras.h"
#include "math.h"

/*
    The results in the following random test cases were computed with the
    naive implementation. Doing a live comparison with large values against
    the naive implementation would take too much time.
*/
static const slong testdata[][4] =
{
    /* h, k,  p/q */
    {WORD(20816815), WORD(29229), WORD(-10669), WORD(87687)},
    {WORD(-481962612), WORD(709105), WORD(-910639), WORD(141821)},
    {WORD(-70965), WORD(3384), WORD(1785), WORD(752)},
    {WORD(1899905), WORD(6657), WORD(-43795), WORD(5706)},
    {WORD(-1893), WORD(511167), WORD(-3411568), WORD(170389)},
    {WORD(1417295), WORD(10180), WORD(3543), WORD(4072)},
    {WORD(-1149), WORD(9350), WORD(6971), WORD(9350)},
    {WORD(-15520), WORD(22977640), WORD(-70331425), WORD(574441)},
    {WORD(3339), WORD(9873153), WORD(270746882), WORD(1097017)},
    {WORD(470645896), WORD(71754), WORD(-21713), WORD(107631)},
    {WORD(1153), WORD(1332403), WORD(258755243), WORD(2664806)},
    {WORD(-501576), WORD(292801), WORD(269095), WORD(292801)},
    {WORD(1861), WORD(34440), WORD(-723059), WORD(206640)},
    {WORD(-4278761), WORD(239321), WORD(791947), WORD(239321)},
    {WORD(9414763), WORD(30776409), WORD(-93285463), WORD(92329227)},
    {WORD(4872687), WORD(2199), WORD(146), WORD(733)},
    {WORD(-22349505), WORD(60581653), WORD(27694241), WORD(60581653)},
    {WORD(85739724), WORD(9289), WORD(961), WORD(2654)},
    {WORD(-5616), WORD(124023), WORD(-31447), WORD(41341)},
    {WORD(99382204), WORD(1378843), WORD(-2537405), WORD(2757686)},
    {WORD(1903), WORD(15842), WORD(102), WORD(89)},
    {WORD(-907226), WORD(5818), WORD(5608), WORD(2909)},
    {WORD(-948920), WORD(4768), WORD(-4815), WORD(1192)},
    {WORD(-352220914), WORD(15390287), WORD(-171358081), WORD(30780574)},
    {WORD(-159206), WORD(3028284), WORD(12921745), WORD(4542426)},
    {WORD(61951448), WORD(1624), WORD(-341), WORD(406)},
    {WORD(-49167), WORD(2092), WORD(-32915), WORD(4184)},
    {WORD(-20878222), WORD(586303210), WORD(-530581301), WORD(293151605)},
    {WORD(-1435637), WORD(3483), WORD(-4787), WORD(20898)},
    {WORD(-1129797), WORD(171620), WORD(238211), WORD(68648)},
    {WORD(-177095), WORD(2914), WORD(1132), WORD(1457)},
    {WORD(-343227551), WORD(1509), WORD(-3289), WORD(4527)},
    {WORD(57497376), WORD(1351), WORD(373), WORD(2702)},
    {WORD(3350543), WORD(5771893), WORD(-51196457), WORD(5771893)},
    {WORD(-44408), WORD(1670), WORD(367), WORD(1670)},
    {WORD(-4139), WORD(59959), WORD(-286689), WORD(119918)},
    {WORD(7397588), WORD(16695), WORD(-41627), WORD(20034)},
    {WORD(-78900791), WORD(10792), WORD(-30905), WORD(21584)},
    {WORD(-1204294), WORD(10134), WORD(-8945), WORD(30402)},
    {WORD(27649424), WORD(57014291), WORD(731583513), WORD(114028582)},
    {WORD(3275043), WORD(436410815), WORD(2018428417), WORD(174564326)},
#if FLINT64  /* skip on 32 bit only because of the literals */
    {WORD(61247), WORD(81381215), WORD(3622491319), WORD(32552486)},
    {WORD(-52118), WORD(125095621), WORD(-24931204413), WORD(125095621)},
    {WORD(201446493), WORD(951783261), WORD(2467429915), WORD(634522174)},
    {WORD(176112), WORD(72187934), WORD(2692844825), WORD(72187934)},
    {WORD(1272), WORD(8722219), WORD(9972821075), WORD(17444438)},
#endif
    {0, 0, 0, 0}
};

int main(void)
{
    fmpz_t hh, kk;
    fmpq_t s1, s2;
    slong i, h, k;

    FLINT_TEST_INIT(state);

    flint_printf("dedekind_sum....");
    fflush(stdout);
    
    fmpz_init(hh);
    fmpz_init(kk);
    fmpq_init(s1);
    fmpq_init(s2);

    for (k = -200; k < 200; k++)
    {
        for (h = -200; h < 200; h++)
        {
            fmpz_set_si(hh, h);
            fmpz_set_si(kk, k);

            arith_dedekind_sum(s1, hh, kk);
            arith_dedekind_sum_naive(s2, hh, kk);

            if (!fmpq_equal(s1, s2))
            {
                flint_printf("FAIL:\n");
                flint_printf("s(%wd,%wd)\n", h, k);
                flint_printf("s1: "); fmpq_print(s1); flint_printf("\n");
                flint_printf("s2: "); fmpq_print(s2); flint_printf("\n");
                abort();
            }
        }
    }

    /* Test large values, 10-30 bits */
    for (i = 0; testdata[i][0] != 0; i++)
    {
        h = testdata[i][0];
        k = testdata[i][1];

        fmpz_set_si(hh, h);
        fmpz_set_si(kk, k);

        arith_dedekind_sum(s1, hh, kk);

        fmpz_set_si(fmpq_numref(s2), testdata[i][2]);
        fmpz_set_si(fmpq_denref(s2), testdata[i][3]);

        if (!fmpq_equal(s1, s2))
        {
            flint_printf("FAIL:\n");
            flint_printf("s(%wd,%wd)\n", h, k);
            flint_printf("s1: "); fmpq_print(s1); flint_printf("\n");
            flint_printf("s2: "); fmpq_print(s2); flint_printf("\n");
            abort();
        }
    }

    /* Check a large value computed with Pari */
    fmpz_set_ui(hh, 1);
    fmpz_mul_2exp(hh, hh, 1000);
    fmpz_add_ui(hh, hh, 1);
    fmpz_set_ui(kk, 1);
    fmpz_mul_2exp(kk, kk, 1001);
    fmpz_add_ui(kk, kk, 1);

    arith_dedekind_sum(s1, hh, kk);
    if ((fmpz_fdiv_ui(fmpq_numref(s1), 1000000000) != 906445312) ||
        (fmpz_fdiv_ui(fmpq_denref(s1), 1000000000) != 8416259))
    {
        flint_printf("Wrong large value:\n");
        fmpq_print(s1);
        flint_printf("\n");
        abort();
    }


    /* Just check that nothing crashes with bignums */
    for (i = 0; i < 1000; i++)
    {
        fmpz_randtest(hh, state, 300);
        fmpz_randtest(kk, state, 300);

        arith_dedekind_sum(s1, hh, kk);
    }

    fmpz_clear(hh);
    fmpz_clear(kk);
    fmpq_clear(s1);
    fmpq_clear(s2);

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}
