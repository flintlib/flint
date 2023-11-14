/*
    Copyright (C) 2012 William Hart

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

int main(int argc, char *argv[])
{
    int result;
    slong exp, val, num, B1, B2, i;

    if (argc != 6)
    {
       printf("usage: %s exp val num B1 B2\n", argv[0]);
       printf("where D = -4*(10^exp + i) for i in [val..val + num)\n");
       printf("with prime bound B1 and large prime bound B2\n");
       return 1;
    }

    exp = atol(argv[1]);
    val = atol(argv[2]);
    num = atol(argv[3]);
    B1 = atol(argv[4]);
    B2 = atol(argv[5]);

    for (i = 0; i < num; i++)
    {
        fmpz_t D, exponent;
        slong e;

        fmpz_init(D);
        fmpz_init(exponent);

        fmpz_set_ui(D, 10);
        fmpz_pow_ui(D, D, exp);
        fmpz_add_ui(D, D, val + i);
        fmpz_mul_2exp(D, D, 2);
        fmpz_neg(D, D);

        if (qfb_exponent_grh(exponent, D, B1, B2))
        {
           printf("Discriminant: "); fmpz_print(D); printf("\n");
           printf("Exponent: "); fmpz_print(exponent); printf("\n\n");
        }

        fmpz_clear(D);
        fmpz_clear(exponent);
    }

    _fmpz_cleanup();
    return 0;
}
