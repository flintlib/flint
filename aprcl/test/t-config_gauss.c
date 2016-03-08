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
    int i;
    FLINT_TEST_INIT(state);
   
    flint_printf("config_gauss....");
    fflush(stdout);

    for (i = 0; i < 100; i++)
    {
        fmpz_t n, s2;
        aprcl_config conf;

        fmpz_init(s2);
        fmpz_init(n);
        fmpz_randtest_unsigned(n, state, 200);

        config_gauss_init(conf, n);

        fmpz_mul(s2, conf->s, conf->s);
        if (fmpz_cmp(s2, n) <= 0)
        {
            flint_printf("FAIL:\n");
            flint_printf("s^2 = ");
            fmpz_print(s2);
            flint_printf(" <= ");
            fmpz_print(n);
            abort();
        }

        config_gauss_clear(conf);
        fmpz_clear(n);
        fmpz_clear(s2);
    }

    FLINT_TEST_CLEANUP(state);
    
    flint_printf("PASS\n");
    return 0;
}

