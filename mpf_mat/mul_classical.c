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

    Copyright (C) 2010,2011 Fredrik Johansson
    Copyright (C) 2014 Abhinav Baid

******************************************************************************/

#include "mpf_mat.h"

void
mpf_mat_mul_classical(mpf_mat_t C, const mpf_mat_t A, const mpf_mat_t B,
                      mpf_rnd_t rnd)
{
    slong ar, bc, br;
    slong i, j, k, exp;
    mpf_t tmp, rtmp;

    ar = A->r;
    br = B->r;
    bc = B->c;

    if (C == A || C == B)
    {
        mpf_mat_t t;
        mpf_mat_init(t, ar, bc, C->prec);
        mpf_mat_mul_classical(t, A, B, rnd);
        mpf_mat_swap(C, t);
        mpf_mat_clear(t);
        return;
    }

    if (C->r != ar || C->c != bc)
    {
        flint_printf
            ("Exception (mpf_mat_mul_classical). Incompatible dimensions.\n");
        abort();
    }

    if (br == 0)
    {
        mpf_mat_zero(C);
        return;
    }

    mpf_init2(tmp, C->prec);
    mpf_init2(rtmp, C->prec);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mpf_mul(mpf_mat_entry(C, i, j), mpf_mat_entry(A, i, 0),
                    mpf_mat_entry(B, 0, j));
            if (rnd != MPF_RNDZ)
            {
                mpf_get_d_2exp(&exp, mpf_mat_entry(C, i, j));
                mpf_set_ui(rtmp, 1);
                mpf_mul_2exp(rtmp, rtmp,
                             exp - mpf_get_prec(mpf_mat_entry(C, i, j)));
                if (rnd == MPF_RNDU)
                {
                    mpf_add(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j),
                            rtmp);
                }
                else
                {
                    mpf_sub(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j),
                            rtmp);
                }
            }

            for (k = 1; k < br; k++)
            {
                mpf_mul(tmp, mpf_mat_entry(A, i, k), mpf_mat_entry(B, k, j));
                if (rnd != MPF_RNDZ)
                {
                    mpf_get_d_2exp(&exp, tmp);
                    mpf_set_ui(rtmp, 1);
                    mpf_mul_2exp(rtmp, rtmp, exp - mpf_get_prec(tmp));
                    if (rnd == MPF_RNDU)
                    {
                        mpf_add(tmp, tmp, rtmp);
                    }
                    else
                    {
                        mpf_sub(tmp, tmp, rtmp);
                    }
                }
                mpf_add(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j), tmp);
                if (rnd != MPF_RNDZ)
                {
                    mpf_get_d_2exp(&exp, mpf_mat_entry(C, i, j));
                    mpf_set_ui(rtmp, 1);
                    mpf_mul_2exp(rtmp, rtmp,
                                 exp - mpf_get_prec(mpf_mat_entry(C, i, j)));
                    if (rnd == MPF_RNDU)
                    {
                        mpf_add(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j),
                                rtmp);
                    }
                    else
                    {
                        mpf_sub(mpf_mat_entry(C, i, j), mpf_mat_entry(C, i, j),
                                rtmp);
                    }
                }
            }
        }
    }

    mpf_clears(tmp, rtmp, '\0');
}
