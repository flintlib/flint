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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_mat.h"
#include "longlong.h"

void
fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
                                                const fmpz_mat_t B)
{
    len_t ar, bc, br;
    len_t i, j, k;

    fmpz a, b;
    mpz_t t;

    mp_limb_t au, bu;
    mp_limb_t pos[3];
    mp_limb_t neg[3];

    ar = A->r;
    br = B->r;
    bc = B->c;

    mpz_init(t);

    for (i = 0; i < ar; i++)
    {
        for (j = 0; j < bc; j++)
        {
            mpz_set_ui(t, 0UL);

            pos[2] = pos[1] = pos[0] = neg[2] = neg[1] = neg[0] = 0UL;

            for (k = 0; k < br; k++)
            {
                a = A->rows[i][k];
                b = B->rows[k][j];

                if (a == 0 || b == 0)
                    continue;

                if (!COEFF_IS_MPZ(a))   /* a is small */
                {
                    if (!COEFF_IS_MPZ(b))  /* both are small */
                    {
                        au = FLINT_ABS(a);
                        bu = FLINT_ABS(b);

                        umul_ppmm(au, bu, au, bu);

                        if ((a ^ b) >= 0L)
                            add_sssaaaaaa(pos[2], pos[1], pos[0],
                                          pos[2], pos[1], pos[0], 0, au, bu);
                        else
                            add_sssaaaaaa(neg[2], neg[1], neg[0],
                                          neg[2], neg[1], neg[0], 0, au, bu);
                    }
                    else
                    {
                        if (a >= 0)
                            mpz_addmul_ui(t, COEFF_TO_PTR(b), a);
                        else
                            mpz_submul_ui(t, COEFF_TO_PTR(b), -a);
                    }
                }
                else if (!COEFF_IS_MPZ(b))  /* b is small */
                {
                    if (b >= 0)
                        mpz_addmul_ui(t, COEFF_TO_PTR(a), b);
                    else
                        mpz_submul_ui(t, COEFF_TO_PTR(a), -b);
                }
                else
                {
                    mpz_addmul(t, COEFF_TO_PTR(a), COEFF_TO_PTR(b));
                }
            }

            if (mpz_sgn(t) != 0 || pos[2] || neg[2] || pos[1] || neg[1])
            {
                __mpz_struct r;

                r._mp_size = pos[2] ? 3 : (pos[1] ? 2 : pos[0] != 0);
                r._mp_alloc = r._mp_size;
                r._mp_d = pos;

                mpz_add(t, t, &r);

                r._mp_size = neg[2] ? 3 : (neg[1] ? 2 : neg[0] != 0);
                r._mp_alloc = r._mp_size;
                r._mp_d = neg;

                mpz_sub(t, t, &r);

                fmpz_set_mpz(fmpz_mat_entry(C, i, j), t);
            }
            else
            {
                if (neg[0] > pos[0])
                {
                    fmpz_set_ui(fmpz_mat_entry(C, i, j), neg[0] - pos[0]);
                    fmpz_neg(fmpz_mat_entry(C, i, j), fmpz_mat_entry(C, i, j));
                }
                else
                {
                    fmpz_set_ui(fmpz_mat_entry(C, i, j), pos[0] - neg[0]);
                }
            }
        }
    }

    mpz_clear(t);
}
