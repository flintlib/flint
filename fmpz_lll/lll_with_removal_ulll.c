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

    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "fmpz_lll.h"

int
fmpz_lll_with_removal_ulll(fmpz_mat_t FM, fmpz_mat_t UM, slong new_size,
                           const fmpz_t gs_B, const fmpz_lll_t fl)
{
    int newd;
    if (fl->rt == Z_BASIS)
    {
        slong r, c, mbits, prev_mbits;
        int full_prec = 1, done = 0, is_U_I;
        fmpz_mat_t U, trunc_data;

        r = FM->r;
        c = FM->c;
        mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
        prev_mbits = mbits;

        fmpz_mat_init(U, r, r);
        fmpz_mat_init(trunc_data, r, c);

        if (mbits > new_size)
        {
            full_prec = 0;

            /* do some truncating */
            fmpz_mat_scalar_tdiv_q_2exp(trunc_data, FM,
                                        (ulong) (mbits - new_size));

            /* set U to the identity */
            fmpz_mat_one(U);
        }
        else
        {
            full_prec = 1;
        }

        while (done == 0)
        {
            if (full_prec == 0)
            {
                if (new_size > FLINT_BITS)
                {
                    fmpz_lll_with_removal_ulll(trunc_data, U, new_size / 2,
                                               gs_B, fl);
                }
                else
                {
                    fmpz_lll_wrapper_with_removal_knapsack(trunc_data, U, gs_B,
                                                           fl);
                }
            }
            else
            {
                newd =
                    fmpz_lll_wrapper_with_removal_knapsack(FM, UM, gs_B, fl);
            }

            if (full_prec == 1)
                done = 1;
            else
            {
                is_U_I = fmpz_mat_is_one(U);

                if (is_U_I == 0)
                {
                    fmpz_mat_mul(FM, U, FM);
                    if (UM != NULL)
                    {
                        fmpz_mat_mul(UM, U, UM);
                    }
                }

                mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
                /* make this condition better? */
                if ((mbits - new_size > 0)
                    && (mbits <= prev_mbits - (slong) (new_size / 4))
                    && is_U_I == 0)
                {
                    /* do some truncating */
                    fmpz_mat_scalar_tdiv_q_2exp(trunc_data, FM,
                                                (ulong) (mbits - new_size));
                    /* set U to the identity */
                    fmpz_mat_one(U);
                }
                else
                {
                    /* can switch to FM, no need for a new identity */
                    full_prec = 1;
                }

                prev_mbits = mbits;
            }
        }

        fmpz_mat_clear(trunc_data);
        fmpz_mat_clear(U);
    }
    else
    {
        newd = fmpz_lll_wrapper_with_removal_knapsack(FM, UM, gs_B, fl);
    }
    return newd;
}
