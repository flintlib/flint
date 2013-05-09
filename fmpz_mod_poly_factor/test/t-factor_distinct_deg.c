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

    Copyright (C) 2007 David Howden
    Copyright (C) 2007, 2008, 2009, 2010 William Hart
    Copyright (C) 2008 Richard Howell-Peak
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova

******************************************************************************/

#include <stdlib.h>
#include "fmpz_mod_poly_factor.h"
#include "ulong_extras.h"

int
main(void)
{
    int iter;
    flint_rand_t state;
    flint_randinit(state);

    printf("factor_distinct_deg....");
    fflush(stdout);

    for (iter = 0; iter < 200; iter++)
    {
        fmpz_mod_poly_t poly1, poly, q, r, product;
        fmpz_mod_poly_factor_t res;
        fmpz_t modulus;
        len_t i, length, num;
        len_t *degs;

        fmpz_init(modulus);
        fmpz_set_ui(modulus, n_randtest_prime(state, 0));

        fmpz_mod_poly_init(poly1, modulus);
        fmpz_mod_poly_init(poly, modulus);
        fmpz_mod_poly_init(q, modulus);
        fmpz_mod_poly_init(r, modulus);

        fmpz_mod_poly_zero(poly1);
        fmpz_mod_poly_set_coeff_ui(poly1, 0, 1);

        length = n_randint(state, 7) + 2;
        do
        {
            fmpz_mod_poly_randtest(poly, state, length);
            if (poly->length)
                fmpz_mod_poly_make_monic(poly, poly);
        }
        while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly)));

        fmpz_mod_poly_mul(poly1, poly1, poly);

        num = n_randint(state, 5) + 1;

        for (i = 1; i < num; i++)
        {
            do
            {
                length = n_randint(state, 7) + 2;
                fmpz_mod_poly_randtest(poly, state, length);
                if (poly->length)
                {
                    fmpz_mod_poly_make_monic(poly, poly);
                    fmpz_mod_poly_divrem(q, r, poly1, poly);
                }
            }
            while ((poly->length < 2) || (!fmpz_mod_poly_is_irreducible(poly))
                   || (r->length == 0));

            fmpz_mod_poly_mul(poly1, poly1, poly);
        }

        if (!(degs = flint_malloc((poly1->length - 1) * sizeof(len_t))))
        {
            printf("Fatal error: not enough memory.");
            abort();
        }
        fmpz_mod_poly_factor_init(res);
        fmpz_mod_poly_factor_distinct_deg(res, poly1, &degs);

        fmpz_mod_poly_init(product, &poly1->p);
        fmpz_mod_poly_set_coeff_ui(product, 0, 1);
        for (i = 0; i < res->num; i++)
            fmpz_mod_poly_mul(product, product, res->poly + i);

        fmpz_mod_poly_scalar_mul_fmpz(product, product,
                                      &(poly1->coeffs[poly1->length - 1]));

        if (!fmpz_mod_poly_equal(poly1, product))
        {
            printf
                ("Error: product of factors does not equal to the original polynomial\n");
            printf("poly:\n");
            fmpz_mod_poly_print(poly1);
            printf("\n");
            printf("product:\n");
            fmpz_mod_poly_print(product);
            printf("\n");
            abort();
        }

        flint_free(degs);
        fmpz_clear(modulus);
        fmpz_mod_poly_clear(product);
        fmpz_mod_poly_clear(q);
        fmpz_mod_poly_clear(r);
        fmpz_mod_poly_clear(poly1);
        fmpz_mod_poly_clear(poly);
        fmpz_mod_poly_factor_clear(res);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
