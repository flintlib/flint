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
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "ulong_extras.h"

int
main(void)
{
    int i;
    flint_rand_t state;
    flint_randinit(state);

    printf("bsplit_sum_pq....");
    fflush(stdout);

    for (i = 0; i < 10000; i++)
    {
        fmpq *pq;
        fmpq_t s1, s2, pqp;
        fmpq_bsplit_t sum;
        len_t k, n;

        n = n_randint(state, 40);

        pq = _fmpq_vec_init(n);

        fmpq_init(s1);
        fmpq_init(s2);
        fmpq_init(pqp);

        for (k = 0; k < n; k++) fmpq_randtest(pq + k, state, 10);

        fmpq_bsplit_init(sum);
        fmpq_bsplit_sum_pq(sum, pq, 0, n);
        fmpq_bsplit_get_fmpq(s1, sum);

        fmpq_zero(s2);
        fmpq_one(pqp);

        for (k = 0; k < n; k++)
        {
            fmpq_mul(pqp, pqp, pq + k);
            fmpq_add(s2, s2, pqp);
        }

        if (!fmpq_is_canonical(s1) || !fmpq_equal(s1, s2))
        {
            printf("FAIL\n");
            printf("(p/q) = ");
            for (k = 0; k < n; k++) fmpq_print(pq+k), printf(" "); printf("\n");
            printf("s1: "); fmpq_print(s1); printf("\n"); 
            printf("s2: "); fmpq_print(s2); printf("\n"); 
            abort();
        }

        /* Check numerical evaluation */
        {
            mpfr_t f1, f2;
            mpfr_prec_t prec;

            prec = 5 + n_randint(state, 1000);

            mpfr_init2(f1, prec);
            mpfr_init2(f2, prec);

            fmpq_bsplit_get_mpfr(f1, sum);
            fmpq_get_mpfr(f2, s1, MPFR_RNDN);

            mpfr_sub(f1, f1, f2, MPFR_RNDN);
            if (!mpfr_zero_p(f1) &&
                !(mpfr_get_exp(f1) <= mpfr_get_exp(f2) - prec + 3))
            {
                printf("FAIL: numerical evaluation\n");
                printf("%ld, %ld, %ld\n", prec, mpfr_get_exp(f1),
                    mpfr_get_exp(f2) - prec + 3);
                abort();
            }


            mpfr_clear(f1);
            mpfr_clear(f2);
        }

        fmpq_bsplit_clear(sum);
        fmpq_clear(s1);
        fmpq_clear(s2);
        fmpq_clear(pqp);
        _fmpq_vec_clear(pq, n);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
