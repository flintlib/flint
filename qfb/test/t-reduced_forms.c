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

    Copyright (C) 2011 Sebastian Pancratz

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <mpir.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "qfb.h"

int main(void)
{
    int result;
    flint_rand_t state;
    qfb * forms;
    long i, j, num;

    printf("reduced_forms....");
    fflush(stdout);

    flint_randinit(state);

    for (i = 1; i < 1000; i++) 
    {
        num = qfb_reduced_forms(&forms, -i);
        
        if (num)
        {
            printf("d = %ld: \n", -i);
            for (j = 0; j < num; j++)
                printf("  (%ld, %ld, %ld)\n", forms[j].a, forms[j].b, forms[j].c);
        }
        
        result = 1;
        if (!result)
        {
            abort();
        }

        flint_free(forms);
    }

    flint_randclear(state);
    _fmpz_cleanup();
    printf("PASS\n");
    return 0;
}
