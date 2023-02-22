/*
    Copyright (C) 2012 Sebastian Pancratz
    Copyright (C) 2013 Tom Bachmann (C++ adaptation)

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

#include <cstdio>
#include <ctime>
#include "fmpz_mod_polyxx.h"

using namespace std;
using namespace flint;

int main(void)
{
    const slong n = 12376;
    const slong N = n / 26;

    frandxx state;

    fmpzxx m(17);
    m = m.pow(26u);

    fmpz_mod_polyxx A(m), B(m);

    A.set_coeff(3, 5);
    A.set_coeff(4, 4);

    B.set_coeff(0, 1);
    B.set_coeff(2, 1);
    B.set_coeff(3, 5);
    B.set_coeff(4, 1);
    B.set_coeff(5, 5);
    B.set_coeff(8, 8);
    B.set_coeff(9, 8);
    B.set_coeff(10, 5);
    B.set_coeff(12, 6);
    B.set_coeff(13, 1);

    fmpz_mod_polyxx r(A.pow(3u) * fmpzxx(4) + B.pow(2u) * fmpzxx(27));

    fmpz_mod_poly_vecxx b(N + 1, m);
    
    fmpz_mod_polyxx t = fmpz_mod_polyxx::randtest(m, state, n + 1);

    flint_printf("Radix conversion\n");
    flint_printf("----------------\n");
    flint_printf("  Degree of the radix:     %wd\n", r.degree());
    flint_printf("  Bit size of the modulus: %wd\n", (slong) bits(r.modulus()));
    flint_printf("  Degree of the input:     %wd\n", t.degree());

    clock_t c0 = clock();
    fmpz_mod_poly_radixxx S(r, n + 1);
    clock_t c1 = clock();
    double c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

    flint_printf("  Precomputation:          %fs\n", c);

    c0 = clock();
    b = t.radix(S);
    c1 = clock();
    c  = (double) (c1 - c0) / CLOCKS_PER_SEC;

    flint_printf("  Conversion:              %fs\n", c);

    return 0;
}

