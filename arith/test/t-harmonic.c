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

    Copyright (C) 2010,2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "fmpq.h"
#include "ulong_extras.h"
#include "profiler.h"


void numerical_test(fmpq_t res, len_t n, double ans)
{
    const double tol = 1e-13;
    double err;

    mpq_t tmp;
    mpq_init(tmp);

    arith_harmonic_number(res, n);
    fmpq_get_mpq(tmp, res);
    err = mpq_get_d(tmp) - ans;
    err = FLINT_ABS(err);

    if (err > tol)
    {
        printf("FAIL: %ld %.16f %.16f\n", n, mpq_get_d(tmp), ans);
        abort();
    }

    mpq_clear(tmp);
}

void
mpq_harmonic_balanced(mpq_t res, len_t a, len_t b)
{
    len_t k;
    mpq_t t;

    mpq_init(t);

    if (b - a < 50)
    {
        mpq_set_ui(res, 0, 1UL);
        for (k = a; k <= b; k++)
        {
            mpq_set_ui(t, 1UL, k);
            mpq_add(res, res, t);
        }
    }
    else
    {
        mpq_harmonic_balanced(res, a, (a+b)/2);
        mpq_harmonic_balanced(t, (a+b)/2+1, b);
        mpq_add(res, res, t);
    }

    mpq_clear(t);
}


int main(void)
{
    len_t i;
    mpq_t x, y;
    fmpq_t t;

    printf("harmonic_number....");
    fflush(stdout);

    fmpq_init(t);
    mpq_init(x);
    mpq_init(y);

    for (i = -2; i < 1000; i++)
    {
        mpq_harmonic_balanced(x, 1, i);
        arith_harmonic_number(t, i);
        fmpq_get_mpq(y, t);

        if (!mpq_equal(x, y))
        {
            printf("FAIL: %ld\n", i);
            abort();
        }
    }

    numerical_test(t, 1000, 7.4854708605503449127);
    numerical_test(t, 1001, 7.4864698615493459117);
    numerical_test(t, 1002, 7.4874678655413618797);
    numerical_test(t, 1003, 7.4884648745144426375);

    numerical_test(t, 10000, 9.7876060360443822642);
    numerical_test(t, 10001, 9.7877060260453821642);
    numerical_test(t, 10002, 9.7878060060493813643);
    numerical_test(t, 10003, 9.7879059760583786652);
    numerical_test(t, 10004, 9.7880059360743722677);

    numerical_test(t, 20000, 10.480728217229327573);
    numerical_test(t, 30000, 10.886184992119899362);
    numerical_test(t, 40000, 11.173862897945522882);
    numerical_test(t, 50000, 11.397003949278482638);
    numerical_test(t, 60000, 11.579323839415955783);
    numerical_test(t, 70000, 11.733473328773164956);
    numerical_test(t, 80000, 11.867003828544530692);
    numerical_test(t, 90000, 11.984786169759202469);

    numerical_test(t, 100000, 12.090146129863427947);
    numerical_test(t, 100001, 12.090156129763428947);
    numerical_test(t, 100002, 12.090166129563432947);
    numerical_test(t, 100003, 12.090176129263441947);
    numerical_test(t, 100004, 12.090186128863457946);

    numerical_test(t, 300000, 13.188755085205611713);
    numerical_test(t, 500000, 13.699580042305528322);
    numerical_test(t, 700000, 14.036051993212618803);
    numerical_test(t, 900000, 14.287366262763433338);

    mpq_clear(x);
    mpq_clear(y);
    fmpq_clear(t);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
