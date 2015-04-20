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

    Copyright (C) 2015 Vladimir Glazachev

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "aprcl.h"
#include "fmpz_poly.h"

int main(void)
{
    FLINT_TEST_INIT(state);
   
    flint_printf("unity....");
    fflush(stdout);

    unity_root u1, u2, u3;
    ulong x;
    slong y;

    unity_init(u1, 7);
    unity_init(u2, 7);
    unity_init(u3, 7);

    unity_nth_root(u1, 0);
    unity_nth_root(u2, 3);
    unity_nth_root(u3, 6);

    flint_printf("\n");
    unity_print(u1);
    unity_print(u2);
    unity_print(u3);

    unity_roots_add(u1, u2, u3);

    unity_print(u1);

    unity_roots_mul(u3, u2, u1);

    unity_print(u3);

    unity_roots_mul_sub(u3, u2, u1);

    unity_print(u3);

    unity_clear(u1);
    unity_clear(u2);
    unity_clear(u3);

    unity_init(u3, 9);
    unity_nth_root(u3, 6);
    unity_nth_root(u3, 3);
    unity_nth_root(u3, 1);
    unity_print(u3);
    unity_roots_add(u3, u3, u3);
    unity_print(u3);
    unity_nth_root_inc(u3, 3);
    unity_print(u3);
    unity_roots_reduce_cyclotomic(u3, 3);
    unity_print(u3);

    FLINT_TEST_CLEANUP(state);
   
    flint_printf("PASS\n");
    return 0;
}

