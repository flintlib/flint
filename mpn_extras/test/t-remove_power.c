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
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"
#include "ulong_extras.h"

void test_exact(int d)
{
    int i, j;
    ulong exp, exact;
    mpz_t a, a2, b, c;
    mpz_init(a);
    mpz_init(a2);
    mpz_init(b);
    mpz_init(c);
    for (i=0; i<100; i++)
    {
        for (j=1; j<100; j++)
        {
            exact = i / j;
            mpz_set_ui(a, d);
            mpz_pow_ui(a, a, i);
            mpz_set(a2, a);
            mpz_set_ui(b, d);
            mpz_pow_ui(b, b, j);
            a->_mp_size = flint_mpn_remove_power_ascending(a->_mp_d, a->_mp_size,
                b->_mp_d, b->_mp_size, &exp);
            mpz_pow_ui(b, b, exact);
            mpz_tdiv_q(c, a2, b);
            if (exp != i/j || mpz_cmp(a, c))
            {
                gmp_printf("%d^%d / %d^%d\n", d, i, d, j);
                abort();
            }
        }
    }

    mpz_clear(a);
    mpz_clear(a2);
    mpz_clear(b);
    mpz_clear(c);
}


int main(void)
{
    printf("remove_power....");
    fflush(stdout);

    test_exact(3);
    test_exact(10);
    test_exact(7429);

    printf("PASS\n");
    return 0;
}
