/*
    Copyright (C) 2009, 2010 William Hart
    Copyright (C) 2009, 2010 Andy Novocin
    Copyright (C) 2014 Abhinav Baid

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_lll.h"

void
fmpz_lll_storjohann_ulll(fmpz_mat_t FM, slong new_size, const fmpz_lll_t fl)
{
    if (fl->rt == Z_BASIS)
    {
        slong r, c, mbits, prev_mbits, i, j;
        int full_prec = 1, done = 0, is_U_I;
        fmpz_mat_t U, big_td, trunc_data;
        mpq_t deltax, etax;
        fmpq_t delta, eta;

        r = FM->r;
        c = FM->c;
        mbits = FLINT_ABS(fmpz_mat_max_bits(FM));
        prev_mbits = mbits;

        fmpz_mat_init(big_td, r, c + r);
        fmpz_mat_init(trunc_data, r, c);

        mpq_init(deltax);
        mpq_init(etax);

        fmpq_init(delta);
        fmpq_init(eta);

        mpq_set_d(deltax, fl->delta);
        mpq_set_d(etax, fl->eta);
        fmpq_set_mpq(delta, deltax);
        fmpq_set_mpq(eta, etax);
        mpq_clears(deltax, etax, NULL);

        if (mbits > new_size)
        {
            full_prec = 0;

            /* do some truncating */
            fmpz_mat_scalar_tdiv_q_2exp(trunc_data, FM,
                                        (ulong) (mbits - new_size));

            /* make a large lattice which has identity in one corner and trunc_data in the other */
            for (i = 0; i < r; i++)
            {
                fmpz_one(fmpz_mat_entry(big_td, i, i));

                for (j = r; j < r + c; j++)
                    fmpz_set(fmpz_mat_entry(big_td, i, j),
                             fmpz_mat_entry(trunc_data, i, j - r));
            }
        }
        else
        {
            full_prec = 1;
        }

        while (done == 0)
        {
            if (full_prec == 0)
            {
                fmpz_mat_lll_storjohann(big_td, delta, eta);
            }
            else
            {
                fmpz_mat_lll_storjohann(FM, delta, eta);
            }

            if (full_prec == 1)
                done = 1;
            else
            {
                /* get U and compare it to the identity */
                fmpz_mat_window_init(U, big_td, 0, 0, r, r);
                is_U_I = fmpz_mat_is_one(U);

                if (is_U_I == 0)
                {
                    fmpz_mat_mul(FM, U, FM);
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

                    /* keep with the big_td concept */
                    for (i = 0; i < r; i++)
                    {
                        for (j = 0; j < i; j++)
                            fmpz_zero(fmpz_mat_entry(big_td, i, j));
                        fmpz_one(fmpz_mat_entry(big_td, i, i));

                        for (j = i + 1; j < r; j++)
                            fmpz_zero(fmpz_mat_entry(big_td, i, j));

                        for (j = r; j < r + c; j++)
                            fmpz_set(fmpz_mat_entry
                                     (big_td, i, j),
                                     fmpz_mat_entry(trunc_data, i, j - r));
                    }
                }
                else
                {
                    /* can switch to FM, no need for a new identity */
                    full_prec = 1;
                }

                prev_mbits = mbits;
                fmpz_mat_window_clear(U);
            }
        }

        fmpz_mat_clear(trunc_data);
        fmpz_mat_clear(big_td);
        fmpq_clear(delta);
        fmpq_clear(eta);
    }
    else
    {
        fmpz_lll_with_removal_ulll(FM, NULL, new_size, NULL, fl);
    }
}
