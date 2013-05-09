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
    the naive implementation would take too long.
*/
static const len_t testdata[][4] =
{
    /* h, k,  p/q */
    {20816815L, 29229L, -10669L, 87687L},
    {-481962612L, 709105L, -910639L, 141821L},
    {-70965L, 3384L, 1785L, 752L},
    {1899905L, 6657L, -43795L, 5706L},
    {-1893L, 511167L, -3411568L, 170389L},
    {1417295L, 10180L, 3543L, 4072L},
    {-1149L, 9350L, 6971L, 9350L},
    {-15520L, 22977640L, -70331425L, 574441L},
    {3339L, 9873153L, 270746882L, 1097017L},
    {470645896L, 71754L, -21713L, 107631L},
    {1153L, 1332403L, 258755243L, 2664806L},
    {-501576L, 292801L, 269095L, 292801L},
    {1861L, 34440L, -723059L, 206640L},
    {-4278761L, 239321L, 791947L, 239321L},
    {9414763L, 30776409L, -93285463L, 92329227L},
    {4872687L, 2199L, 146L, 733L},
    {-22349505L, 60581653L, 27694241L, 60581653L},
    {85739724L, 9289L, 961L, 2654L},
    {-5616L, 124023L, -31447L, 41341L},
    {99382204L, 1378843L, -2537405L, 2757686L},
    {1903L, 15842L, 102L, 89L},
    {-907226L, 5818L, 5608L, 2909L},
    {-948920L, 4768L, -4815L, 1192L},
    {-352220914L, 15390287L, -171358081L, 30780574L},
    {-159206L, 3028284L, 12921745L, 4542426L},
    {61951448L, 1624L, -341L, 406L},
    {-49167L, 2092L, -32915L, 4184L},
    {-20878222L, 586303210L, -530581301L, 293151605L},
    {-1435637L, 3483L, -4787L, 20898L},
    {-1129797L, 171620L, 238211L, 68648L},
    {-177095L, 2914L, 1132L, 1457L},
    {-343227551L, 1509L, -3289L, 4527L},
    {57497376L, 1351L, 373L, 2702L},
    {3350543L, 5771893L, -51196457L, 5771893L},
    {-44408L, 1670L, 367L, 1670L},
    {-4139L, 59959L, -286689L, 119918L},
    {7397588L, 16695L, -41627L, 20034L},
    {-78900791L, 10792L, -30905L, 21584L},
    {-1204294L, 10134L, -8945L, 30402L},
    {27649424L, 57014291L, 731583513L, 114028582L},
    {3275043L, 436410815L, 2018428417L, 174564326L},
#if FLINT64  /* skip on 32 bit only because of the literals */
    {61247L, 81381215L, 3622491319L, 32552486L},
    {-52118L, 125095621L, -24931204413L, 125095621L},
    {201446493L, 951783261L, 2467429915L, 634522174L},
    {176112L, 72187934L, 2692844825L, 72187934L},
    {1272L, 8722219L, 9972821075L, 17444438L},
#endif
    {0, 0, 0, 0}
};

int main(void)
{
    flint_rand_t state;
    fmpz_t hh, kk;
    fmpq_t s1, s2;
    len_t i, h, k;

    printf("dedekind_sum....");
    fflush(stdout);

    flint_randinit(state);
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
                printf("FAIL:\n");
                printf("s(%ld,%ld)\n", h, k);
                printf("s1: "); fmpq_print(s1); printf("\n");
                printf("s2: "); fmpq_print(s2); printf("\n");
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
            printf("FAIL:\n");
            printf("s(%ld,%ld)\n", h, k);
            printf("s1: "); fmpq_print(s1); printf("\n");
            printf("s2: "); fmpq_print(s2); printf("\n");
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
        printf("Wrong large value:\n");
        fmpq_print(s1);
        printf("\n");
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

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
