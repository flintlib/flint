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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <mpir.h>
#include "flint.h"
#include "arith.h"
#include "ulong_extras.h"
#include "profiler.h"


void numerical_test(mpq_t res, long n, double ans)
{
    const double tol = 1e-13;
    double err;

    mpq_harmonic(res, n);
    err = mpq_get_d(res) - ans;
    err = FLINT_ABS(err);

    if (err > tol)
    {
        printf("FAIL: %ld %.16f %.16f\n", n, mpq_get_d(res), ans);
        abort();
    }
}

int main(void)
{
    long i;
    mpq_t x, y, z;

    printf("harmonic....");
    fflush(stdout);

    mpq_init(x);
    mpq_init(y);
    mpq_init(z);

    for (i = -2; i < 1000; i++)
    {
        _mpq_harmonic_odd_balanced(x, i);
        _mpq_harmonic_balanced(y, 1, i);
        mpq_harmonic(z, i);
        if (!mpq_equal(x, y) || !mpq_equal(x, z))
        {
            printf("FAIL: %ld\n", i);
            abort();
        }
    }

    numerical_test(x, 1000, 7.4854708605503449127);
    numerical_test(x, 1001, 7.4864698615493459117);
    numerical_test(x, 1002, 7.4874678655413618797);
    numerical_test(x, 1003, 7.4884648745144426375);

    numerical_test(x, 10000, 9.7876060360443822642);
    numerical_test(x, 10001, 9.7877060260453821642);
    numerical_test(x, 10002, 9.7878060060493813643);
    numerical_test(x, 10003, 9.7879059760583786652);
    numerical_test(x, 10004, 9.7880059360743722677);

    numerical_test(x, 20000, 10.480728217229327573);
    numerical_test(x, 30000, 10.886184992119899362);
    numerical_test(x, 40000, 11.173862897945522882);
    numerical_test(x, 50000, 11.397003949278482638);
    numerical_test(x, 60000, 11.579323839415955783);
    numerical_test(x, 70000, 11.733473328773164956);
    numerical_test(x, 80000, 11.867003828544530692);
    numerical_test(x, 90000, 11.984786169759202469);

    numerical_test(x, 100000, 12.090146129863427947);
    numerical_test(x, 100001, 12.090156129763428947);
    numerical_test(x, 100002, 12.090166129563432947);
    numerical_test(x, 100003, 12.090176129263441947);
    numerical_test(x, 100004, 12.090186128863457946);

    numerical_test(x, 300000, 13.188755085205611713);
    numerical_test(x, 500000, 13.699580042305528322);
    numerical_test(x, 700000, 14.036051993212618803);
    numerical_test(x, 900000, 14.287366262763433338);

    mpq_clear(x);
    mpq_clear(y);
    mpq_clear(z);

    printf("PASS\n");
    return 0;
}
