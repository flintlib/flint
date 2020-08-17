/*
    Copyright (C) 2011 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mat.h"
#include "longlong.h"

void
fmpz_mat_mul_classical_inline(fmpz_mat_t C, const fmpz_mat_t A,
                                                const fmpz_mat_t B)
{
    slong ar, bc, br;
    slong i, j, k;

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
            flint_mpz_set_ui(t, UWORD(0));

            pos[2] = pos[1] = pos[0] = neg[2] = neg[1] = neg[0] = UWORD(0);

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

                        if ((a ^ b) >= WORD(0))
                            add_sssaaaaaa(pos[2], pos[1], pos[0],
                                          pos[2], pos[1], pos[0], 0, au, bu);
                        else
                            add_sssaaaaaa(neg[2], neg[1], neg[0],
                                          neg[2], neg[1], neg[0], 0, au, bu);
                    }
                    else
                    {
                        if (a >= 0)
                            flint_mpz_addmul_ui(t, COEFF_TO_PTR(b), a);
                        else
                            flint_mpz_submul_ui(t, COEFF_TO_PTR(b), -a);
                    }
                }
                else if (!COEFF_IS_MPZ(b))  /* b is small */
                {
                    if (b >= 0)
                        flint_mpz_addmul_ui(t, COEFF_TO_PTR(a), b);
                    else
                        flint_mpz_submul_ui(t, COEFF_TO_PTR(a), -b);
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
                    fmpz_neg_ui(fmpz_mat_entry(C, i, j), neg[0] - pos[0]);
                else
                    fmpz_set_ui(fmpz_mat_entry(C, i, j), pos[0] - neg[0]);
            }
        }
    }

    mpz_clear(t);
}
