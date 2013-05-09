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
    len_t i, bit_bound;
    mp_limb_t prime, res;
    fmpz_t x, y, prod;

    if (argc != 2)
    {
        printf("Syntax: crt <integer>\n");
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

        printf("residue mod %lu = %lu; reconstruction = ", prime, res);
        fmpz_print(y);
        printf("\n");

        fmpz_mul_ui(prod, prod, prime);
    }

    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_clear(prod);

    return EXIT_SUCCESS;
}
