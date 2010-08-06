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

    Copyright (C) 2010 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include <mpfr.h>
#include <math.h>
#include "flint.h"
#include "mpfr_vec.h"
#include "mpfr_poly.h"
#include "ulong_extras.h"

int
main(void)
{
    int i, result;
    printf("FHT....");
    fflush(stdout);

    mpfr_poly_randinit();

    for (i = 0; i < 1000; i++)
    {
        mpfr_poly_t a, b;
        long j, n = n_randint(10);
        long length = (1L << n);
        mpfr_prec_t prec = n_randint(100) + 50 * MPFR_PREC_MIN;

        mpfr_poly_init2(a, length, prec);
        mpfr_poly_init2(b, length, prec);
        mpfr_poly_randtest(a, length);

        _mpfr_vec_copy(b->coeffs, a->coeffs, length);

        _mpfr_poly_FHT(a->coeffs, n, prec);
        _mpfr_poly_revbin(a->coeffs, n);

        _mpfr_poly_FHT(a->coeffs, n, prec);
        _mpfr_poly_revbin(a->coeffs, n);

        _mpfr_poly_scale(a->coeffs, n);

        for (j = 0; j < length; j++)
        {
            double d;
            mpfr_sub(a->coeffs + j, a->coeffs + j, b->coeffs + j, GMP_RNDN);
            d = mpfr_get_d(a->coeffs + j, GMP_RNDN);
            result = (fabs(d) < 0.1);
            if (!result)
            {
                printf("FAIL:\n");
                printf("d = %f\n", d);
                printf("length = %ld, j = %ld, prec = %ld\n", length, j,
                       (long) prec);
                abort();
            }
        }

        mpfr_poly_clear(a);
        mpfr_poly_clear(b);
    }

    mpfr_poly_randclear();
    printf("PASS\n");
    return 0;
}
