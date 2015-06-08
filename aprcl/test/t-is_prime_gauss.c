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

int main(void)
{
    int i, j;
    FLINT_TEST_INIT(state);
   
    flint_printf("is_prime_gauss....");
    fflush(stdout);

    unity_zpq gauss;
    unity_zpq gausspower;
    unity_zpq gausssigma;

    ulong p = 2;
    ulong q = 11;

    fmpz_t n;
    fmpz_init_set_ui(n, 1026018272686390600);

    flint_printf("RESULT = %i\n", is_prime_gauss(n));

    fmpz_clear(n);

    FLINT_TEST_CLEANUP(state);

    flint_printf("NO TEST\n");
    return 0;
}

