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

    Copyright (C) 2008, 2009, William Hart 
    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"


void fmpz_comb_temp_clear(fmpz_comb_temp_t temp)
{
    long n, i, j;

    n = temp->n;
    j = (1L << (n - 1));

	for (i = 0; i < n; i++)
    {
        _fmpz_vec_clear(temp->comb_temp[i], j);
        j /= 2;
    }

	flint_free(temp->comb_temp);

    fmpz_clear(temp->temp);
    fmpz_clear(temp->temp2);
}

void
fmpz_comb_clear(fmpz_comb_t comb)
{
    long i, j, n;

    n = comb->n;

    /* Size of top level */
    j = (1L << (n - 1));

    /* Clear arrays at each level  */
    for (i = 0; i < n; i++)
    {
        _fmpz_vec_clear(comb->comb[i], j);
        _fmpz_vec_clear(comb->res[i], j);
        j /= 2;
    }
	
	if (n)
	{
        flint_free(comb->comb);
        flint_free(comb->res);
	}

    flint_free(comb->mod);
}
