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

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_vec.h"
#include "fmpz_vec.h"


void fmpz_comb_temp_init(fmpz_comb_temp_t temp, const fmpz_comb_t comb)
{
    long n, i, j;

    /* Allocate space for comb_temp */
    temp->n = n = comb->n;
    temp->comb_temp = (fmpz **) flint_malloc(n * sizeof(fmpz *));

    j = (1L << (n - 1));
    for (i = 0; i < n; i++)
    {
        temp->comb_temp[i] = _fmpz_vec_init(j);
        j /= 2;
    }

    fmpz_init(temp->temp);
    fmpz_init(temp->temp2);
}

void
fmpz_comb_init(fmpz_comb_t comb, mp_limb_t * primes, long num_primes)
{
    long i, j;
    long n, num;
    ulong log_comb, log_res;
    fmpz_t temp, temp2;

    comb->primes = primes;
    comb->num_primes = num_primes;

    n = FLINT_BIT_COUNT(num_primes);
    comb->n = n;

    /* Create nmod_poly modulus information */
	comb->mod = (nmod_t *) flint_malloc(sizeof(nmod_t) * num_primes);
    for (i = 0; i < num_primes; i++)
        nmod_init(&comb->mod[i], primes[i]);

    /* Nothing to do */
	if (n == 0)
        return;

	/* Allocate space for comb and res */
    comb->comb = (fmpz **) flint_malloc(n * sizeof(fmpz *));
    comb->res = (fmpz **) flint_malloc(n * sizeof(fmpz *));

    /* Size of top level */
    j = (1L << (n - 1));

    /* Initialise arrays at each level */
	for (i = 0; i < n; i++)
    {
        comb->comb[i] = _fmpz_vec_init(j);
        comb->res[i] = _fmpz_vec_init(j);
        j /= 2;
    }

	/* Compute products of pairs of primes and place in comb */
    for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
    {
        fmpz_set_ui(comb->comb[0] + j, primes[i]);
        fmpz_mul_ui(comb->comb[0] + j, comb->comb[0] + j, primes[i+1]);
    }

    /* In case number of primes is odd */
    if (i < num_primes)
    {
        fmpz_set_ui(comb->comb[0] + j, primes[i]);
	    i += 2;
	    j++;
	}

    /* Set the rest of the entries on that row of the comb to 1 */
    num = (1L << n);
	for (; i < num; i += 2, j++)
    {
        fmpz_one(comb->comb[0] + j);
    }

    /* Compute rest of comb by multiplying in pairs */
    log_comb = 1;
    num /= 2;
    while (num >= 2)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            fmpz_mul(comb->comb[log_comb] + j, comb->comb[log_comb-1] + i,
                comb->comb[log_comb-1] + i + 1);
        }
        log_comb++;
        num /= 2;
    }

    /* Compute inverses from pairs of primes */
    fmpz_init(temp);
    fmpz_init(temp2);

    for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
    {
        fmpz_set_ui(temp, primes[i]);
        fmpz_set_ui(temp2, primes[i+1]);
        fmpz_invmod(comb->res[0] + j, temp, temp2);
    }

    fmpz_clear(temp);
    fmpz_clear(temp2);

    /* Compute remaining inverses, each level
       combining pairs from the level below */
	log_res = 1;
    num = (1L << (n - 1));

    while (log_res < n)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            fmpz_invmod(comb->res[log_res] + j, comb->comb[log_res-1] + i,
                comb->comb[log_res-1] + i + 1);
        }
        log_res++;
        num /= 2;
    }
}
