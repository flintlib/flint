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

int main(void)
{
    int zero, nonzero;
    mp_bitcnt_t check;
    mpz_t a;
    mpz_t b;

    printf("remove_2exp....");
    fflush(stdout);

    mpz_init(a);
    mpz_init(b);

    for (zero=0; zero<300; zero++)
    {
        for (nonzero=0; nonzero<300; nonzero++)
        {
            mpz_set_ui(a, 1);
            mpz_setbit(a, nonzero);
            mpz_set(b, a);
            mpz_mul_2exp(a, a, zero);
            a->_mp_size = flint_mpn_remove_2exp(a->_mp_d, a->_mp_size, &check);
            if (check != zero || mpz_cmp(a,b))
            {
                gmp_printf("%d %d \n", zero, nonzero);
                abort();
            }
        }
    }

    mpz_clear(a);
    mpz_clear(b);
    printf("PASS\n");
    return 0;
}
