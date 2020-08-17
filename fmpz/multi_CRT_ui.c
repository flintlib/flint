/*
    Copyright (C) 2008, 2009, William Hart 
    Copyright (C) 2010 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <gmp.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "nmod_vec.h"


void
__fmpz_multi_CRT_ui_sign(fmpz_t output, const fmpz_t input,
    const fmpz_comb_t comb, fmpz_t temp)
{
    slong n = comb->n;
    slong p;

    if (n == WORD(0))
    {
        if (fmpz_is_zero(input)) 
        {
            fmpz_zero(output);
            return;
        }

        /* XXX: overflow possible? */
        p = comb->primes[0];
        if ((p - (*input)) < (*input))
            fmpz_set_si(output, (slong) ((*input) - p));
        else
            fmpz_set_ui(output, (*input));
        return;
    }

    fmpz_sub(temp, input, comb->comb[comb->n - 1]);

    if (fmpz_cmpabs(temp, input) <= 0)
        fmpz_set(output, temp);
    else
        fmpz_set(output, input);

    return;
}

void fmpz_multi_CRT_ui(fmpz_t output, mp_srcptr residues,
    const fmpz_comb_t comb, fmpz_comb_temp_t ctemp, int sign)
{
    slong i, j;
    slong n = comb->n;
    slong num;
    slong log_res;
    slong num_primes = comb->num_primes;

    fmpz ** comb_temp = ctemp->comb_temp;
    fmpz * temp = ctemp->temp;
    fmpz * temp2 = ctemp->temp2;

    /* The output is less than a single prime, so just output the result */
    if (num_primes == 1)
    {
        if (sign)
        {
            mp_limb_t p = comb->primes[0];

            if ((p - residues[0]) < residues[0])
                fmpz_set_si(output, residues[0] - p);
            else
                fmpz_set_ui(output, residues[0]);
        }
        else
        {
            fmpz_set_ui(output, residues[0]);
        }
        return;
    }

    /* First layer of reconstruction */
    num = (WORD(1) << n);

    for (i = 0, j = 0; i + 2 <= num_primes; i += 2, j++)
    {
        fmpz_set_ui(temp, residues[i]);
        fmpz_mod_ui(temp2, temp, comb->primes[i+1]);
        fmpz_sub_ui(temp2, temp2, residues[i + 1]);
        fmpz_neg(temp2, temp2);
        fmpz_mul(temp, temp2, comb->res[0] + j);
        fmpz_mod_ui(temp2, temp, comb->primes[i+1]);
        fmpz_mul_ui(temp, temp2, comb->primes[i]); 
        fmpz_add_ui(comb_temp[0] + j, temp, residues[i]);
    }

    if (i < num_primes)
        fmpz_set_ui(comb_temp[0] + j, residues[i]);

    /* Compute other layers of reconstruction */
    num /= 2;
    log_res = 1;

    while (log_res < n)
    {
        for (i = 0, j = 0; i < num; i += 2, j++)
        {
            if (fmpz_is_one(comb->comb[log_res-1] + i + 1))
            {
                if (!fmpz_is_one(comb->comb[log_res-1] + i))
                    fmpz_set(comb_temp[log_res] + j, comb_temp[log_res-1] + i);
            }
            else
            {
                fmpz_mod(temp2, comb_temp[log_res-1] + i,
                    comb->comb[log_res-1] + i + 1);
                fmpz_sub(temp, comb_temp[log_res-1] + i + 1, temp2);
                fmpz_mul(temp2, temp, comb->res[log_res] + j);
                fmpz_mod(temp, temp2, comb->comb[log_res-1] + i + 1);
                fmpz_mul(temp2, temp, comb->comb[log_res-1] + i);
                fmpz_add(comb_temp[log_res] + j, temp2,
                    comb_temp[log_res-1] + i);
            }
        }
        log_res++;
        num /= 2; 
    }

    /* Write out the output */
    if (sign)
        __fmpz_multi_CRT_ui_sign(output, comb_temp[log_res - 1], comb, temp);
    else
        fmpz_set(output, comb_temp[log_res - 1]);
}
