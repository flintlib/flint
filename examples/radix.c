/*
    Copyright (C) 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program to demonstrate some use of the 
    function fmpz_mod_poly_radix() for radix conversion 
    over $\mathbf{Z}/n \mathbf{Z}$.
*/

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <gmp.h>

#include "flint.h"
#include "fmpz_mod_poly.h"

int main(void)
{
    const slong n = 12376;
    const slong N = n / 26;

    clock_t c0, c1;
    double c;

    slong i;
    fmpz_t a, m;
    fmpz_mod_poly_t A, B, r, t;
    fmpz_mod_poly_radix_t S;
    fmpz_mod_poly_struct **b;

    FLINT_TEST_INIT(state);    

    fmpz_init(a);
    fmpz_init(m);

    fmpz_set_ui(m, 17);
    fmpz_pow_ui(m, m, 26);

    fmpz_mod_poly_init(A, m);
    fmpz_mod_poly_init(B, m);
    fmpz_mod_poly_init(r, m);
    fmpz_mod_poly_init(t, m);

    fmpz_mod_poly_set_coeff_ui(A, 3, 5);
    fmpz_mod_poly_set_coeff_ui(A, 4, 4);

    fmpz_mod_poly_set_coeff_ui(B, 0, 1);
    fmpz_mod_poly_set_coeff_ui(B, 2, 1);
    fmpz_mod_poly_set_coeff_ui(B, 3, 5);
    fmpz_mod_poly_set_coeff_ui(B, 4, 1);
    fmpz_mod_poly_set_coeff_ui(B, 5, 5);
    fmpz_mod_poly_set_coeff_ui(B, 8, 8);
    fmpz_mod_poly_set_coeff_ui(B, 9, 8);
    fmpz_mod_poly_set_coeff_ui(B, 10, 5);
    fmpz_mod_poly_set_coeff_ui(B, 12, 6);
    fmpz_mod_poly_set_coeff_ui(B, 13, 1);

    fmpz_mod_poly_pow(r, A, 3);
    fmpz_set_ui(a, 4);
    fmpz_mod_poly_scalar_mul_fmpz(r, r, a);

    fmpz_mod_poly_pow(t, B, 2);
    fmpz_set_ui(a, 27);
    fmpz_mod_poly_scalar_mul_fmpz(t, t, a);

    fmpz_mod_poly_add(r, r, t);

    b = flint_malloc((N + 1) * sizeof(fmpz_mod_poly_struct *));
    for (i = 0; i <= N; i++)
    {
        b[i] = flint_malloc(sizeof(fmpz_mod_poly_struct));
        fmpz_mod_poly_init(b[i], m);
    }

    fmpz_mod_poly_randtest(t, state, n + 1);

    flint_printf("Radix conversion\n");
    flint_printf("----------------\n");
    flint_printf("  Degree of the radix:     %wd\n", fmpz_mod_poly_degree(r));
    flint_printf("  Bit size of the modulus: %wd\n", (slong) fmpz_bits(fmpz_mod_poly_modulus(r)));
    flint_printf("  Degree of the input:     %wd\n", fmpz_mod_poly_degree(t));

    c0 = clock();
    fmpz_mod_poly_radix_init(S, r, n + 1);
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

    flint_printf("  Precomputation:          %fs\n", c);

    c0 = clock();
    fmpz_mod_poly_radix(b, t, S);
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

    flint_printf("  Conversion:              %fs\n", c);

    fmpz_clear(a);
    fmpz_clear(m);
    fmpz_mod_poly_clear(A);
    fmpz_mod_poly_clear(B);
    fmpz_mod_poly_clear(r);
    fmpz_mod_poly_clear(t);
    fmpz_mod_poly_radix_clear(S);

    for (i = 0; i <= N; i++)
    {
        fmpz_mod_poly_clear(b[i]);
        flint_free(b[i]);
    }
    flint_free(b);

    flint_randclear(state);

    return EXIT_SUCCESS;
}

