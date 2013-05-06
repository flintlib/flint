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

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "arith.h"
#include "fmpz.h"
#include "ulong_extras.h"

int main()
{
    fmpz_t s, t;
    long n;

    printf("bernoulli_number_denom....");
    fflush(stdout);

    fmpz_init(s);
    fmpz_init(t);

    for (n = 0; n < 1000; n++)
    {
        arith_bernoulli_number_denom(t, n);
        fmpz_addmul_ui(s, t, n_nth_prime(n+1));
    }

    fmpz_set_str(t, "34549631155954474103407159", 10);

    if (!fmpz_equal(s, t))
    {
        printf("FAIL: Hash disagrees with known value\n");
        abort();
    }

    fmpz_clear(s);
    fmpz_clear(t);

    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
