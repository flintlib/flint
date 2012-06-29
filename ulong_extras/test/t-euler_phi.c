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

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"

int main(void)
{
    int n, k, t;

    printf("euler_phi....");
    fflush(stdout);

    for (n = 0; n < 20 * flint_test_multiplier(); n++)
    {
        t = 0;
        for (k = 1; k <= n; k++)
            t += (n_gcd(n, k) == 1);
        if (t != n_euler_phi(n))
        {
            printf("FAIL:\n");
            printf("phi(%d) = %d, got %lu\n", n, t, n_euler_phi(n)); 
            abort();
        }
    }

    printf("PASS\n");
    return 0;
}
