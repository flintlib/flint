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

void
fmpz_multi_mod_ui_basecase(mp_limb_t * out, fmpz_t in, mp_limb_t * primes,
                           len_t num_primes)
{
    len_t i;
    for (i = 0; i < num_primes; i++)
    {
        out[i] = fmpz_fdiv_ui(in, primes[i]);
    }
}

void
fmpz_multi_mod_ui(mp_limb_t * out, const fmpz_t in, const fmpz_comb_t comb,
    fmpz_comb_temp_t temp)
{
    len_t i, j;
    len_t n = comb->n;
    len_t log_comb;
    len_t stride;
    len_t num;
    len_t num_primes = comb->num_primes;
    fmpz ** comb_temp = temp->comb_temp;

    if (num_primes == 1)
    {
        out[0] = fmpz_fdiv_ui(in, comb->primes[0]);
        return;
    }

    log_comb = n - 1;
   
    /* Find level in comb with entries bigger than the input integer */
    log_comb = 0;
    if (fmpz_sgn(in) < 0)
    {
        while ((fmpz_bits(in) >= fmpz_bits(comb->comb[log_comb]) - 1)
            && (log_comb < comb->n - 1)) log_comb++;
    }
    else
    {
        while (fmpz_cmpabs(in, comb->comb[log_comb]) >= 0 &&
            (log_comb < comb->n - 1))
            log_comb++;
    }

    num = (1L << (n - log_comb - 1));

    /* Set each entry of this level of temp to the input integer */
    for (i = 0; i < num; i++)
    {
        fmpz_set(comb_temp[log_comb] + i, in);
    }

    log_comb--;
    num *= 2;

    /* Fill in other entries of temp by taking entries of temp
        at higher level mod pairs from comb */

    /* keep going until we reach the basecase */
    while (log_comb > FLINT_FMPZ_LOG_MULTI_MOD_CUTOFF)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            fmpz_mod(comb_temp[log_comb] + i, comb_temp[log_comb + 1] + j,
                comb->comb[log_comb] + i);
            fmpz_mod(comb_temp[log_comb] + i + 1, comb_temp[log_comb + 1] + j,
                comb->comb[log_comb] + i + 1);
        }
        num *= 2;
        log_comb--;
    }

    /* Do basecase */
    num /= 2;
    log_comb++;

    stride = (1L << (log_comb + 1));
    for (i = 0, j = 0; j < num_primes; i++, j += stride)
    {
        fmpz_multi_mod_ui_basecase(out + j, comb_temp[log_comb] + i,
            comb->primes + j, FLINT_MIN(stride, num_primes - j));
    }
}
