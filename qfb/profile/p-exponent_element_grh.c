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

    Copyright (C) 2012 William Hart

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

int main(int argc, char *argv[])
{
    int result;
    long exp, val, num, iters, i;

    if (argc != 5)
    {
       printf("usage: %s exp val num iters\n", argv[0]);
       printf("where D = -4*(10^exp + i) for i in [val..val + num)\n");
       printf("and where O(iters) iterations are performed\n");
       return 1;
    }

    exp = atol(argv[1]);
    val = atol(argv[2]);
    num = atol(argv[3]);
    iters = atol(argv[4]);

    for (i = 0; i < num; i++) 
    {
        fmpz_t D, exponent;
        long e;
        
        fmpz_init(D);
        fmpz_init(exponent);
        
        fmpz_set_ui(D, 10);
        fmpz_pow_ui(D, D, exp);
        fmpz_add_ui(D, D, i);
        fmpz_mul_2exp(D, D, 2);
        fmpz_neg(D, D);

        if (qfb_exponent_grh(exponent, D, iters))
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
