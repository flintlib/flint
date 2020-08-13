/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Demo FLINT program for incremental multimodular reduction and
    reconstruction using the Chinese Remainder Theorem.
*/

#include <stdlib.h>
#include <stdio.h>
#include "flint.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main(int argc, char* argv[])
{
    slong i, bit_bound;
    mp_limb_t prime, res;
    fmpz_t x, y, prod;

    if (argc != 2)
    {
        flint_printf("Syntax: crt <integer>\n");
        return EXIT_FAILURE;
    }

    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(prod);

    fmpz_set_str(x, argv[1], 10);
    bit_bound = fmpz_bits(x) + 2;

    fmpz_zero(y);
    fmpz_one(prod);

    prime = 0;
    for (i = 0; fmpz_bits(prod) < bit_bound; i++)
    {
        prime = n_nextprime(prime, 0);
        res = fmpz_fdiv_ui(x, prime);
        fmpz_CRT_ui(y, y, prod, res, prime, 1);

        flint_printf("residue mod %wu = %wu; reconstruction = ", prime, res);
        fmpz_print(y);
        flint_printf("\n");

        fmpz_mul_ui(prod, prod, prime);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(prod);

    return EXIT_SUCCESS;
}
