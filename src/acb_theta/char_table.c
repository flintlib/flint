/*
    Copyright (C) 2024 Jean Kieffer

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_mat.h"
#include "acb_theta.h"

void
acb_theta_char_table(ulong * ch, slong * e, const fmpz_mat_t mat, slong ab)
{
    slong g = sp2gz_dim(mat);
    slong n2 = 1 << (2 * g);
    fmpz_mat_t a, b, c, d, cdt, abt;
    fmpz_mat_t mat_tp;
    fmpz_mat_t diags, alphabeta, alpha, beta;
    fmpz_mat_t Cvec_1, Cvec_2, Lvec;
    fmpz_mat_t coef;
    fmpz_t eps, x;
    slong i;
    slong start = (ab < 0 ? 0 : ab);
    slong end = (ab < 0 ? n2 : ab + 1);

    fmpz_mat_window_init(a, mat, 0, 0, g, g);
    fmpz_mat_window_init(b, mat, 0, g, g, 2 * g);
    fmpz_mat_window_init(c, mat, g, 0, 2 * g, g);
    fmpz_mat_window_init(d, mat, g, g, 2 * g, 2 * g);
    fmpz_mat_init(cdt, g, g);
    fmpz_mat_init(abt, g, g);
    fmpz_mat_init(mat_tp, 2 * g, 2 * g);
    fmpz_mat_init(diags, 2 * g, 1);
    fmpz_mat_init(alphabeta, 2 * g, 1);
    fmpz_mat_init(Cvec_1, g, 1);
    fmpz_mat_init(Cvec_2, g, 1);
    fmpz_mat_init(Lvec, 1, g);
    fmpz_mat_init(coef, 1, 1);
    fmpz_init(eps);
    fmpz_init(x);

    /* Set diags */
    fmpz_mat_transpose(mat_tp, mat);
    fmpz_mat_transpose(cdt, d);
    fmpz_mat_mul(cdt, c, cdt);
    fmpz_mat_transpose(abt, b);
    fmpz_mat_mul(abt, a, abt);
    for (i = 0; i < g; i++)
    {
        fmpz_neg(fmpz_mat_entry(diags, i, 0), fmpz_mat_entry(cdt, i, i));
        fmpz_neg(fmpz_mat_entry(diags, g + i, 0), fmpz_mat_entry(abt, i, i));
    }

    for (ab = start; ab < end; ab++)
    {
        ch[ab - start] = 0;

        /* Turn ab into a 2g x 1 fmpz matrix, set alphabeta = diags + ab */
        for (i = 0; i < 2 * g; i++)
        {
            fmpz_add_si(fmpz_mat_entry(alphabeta, 2 * g - 1 - i, 0),
                fmpz_mat_entry(diags, 2 * g - 1 - i, 0), (ab >> i) % 2);
        }

        /* Perform matrix-vector multiplication */
        fmpz_mat_mul(alphabeta, mat_tp, alphabeta);

        /* Compute eps */
        fmpz_mat_window_init(alpha, alphabeta, 0, 0, g, 1);
        fmpz_mat_window_init(beta, alphabeta, g, 0, 2 * g, 1);

        fmpz_zero(eps);

        fmpz_mat_mul(Cvec_1, c, beta);
        fmpz_mat_mul(Cvec_2, b, alpha);
        fmpz_mat_transpose(Lvec, Cvec_2);
        fmpz_mat_mul(coef, Lvec, Cvec_1);
        fmpz_addmul_ui(eps, fmpz_mat_entry(coef, 0, 0), 2);

        fmpz_mat_mul(Cvec_1, b, alpha);
        fmpz_mat_mul(Cvec_2, d, alpha);
        fmpz_mat_transpose(Lvec, Cvec_2);
        fmpz_mat_mul(coef, Lvec, Cvec_1);
        fmpz_sub(eps, eps, fmpz_mat_entry(coef, 0, 0));

        fmpz_mat_mul(Cvec_1, a, beta);
        fmpz_mat_mul(Cvec_2, c, beta);
        fmpz_mat_transpose(Lvec, Cvec_2);
        fmpz_mat_mul(coef, Lvec, Cvec_1);
        fmpz_sub(eps, eps, fmpz_mat_entry(coef, 0, 0));

        for (i = 0; i < g; i++)
        {
            fmpz_set(fmpz_mat_entry(Lvec, 0, i), fmpz_mat_entry(abt, i, i));
        }
        fmpz_mat_mul(Cvec_1, d, alpha);
        fmpz_mat_mul(Cvec_2, c, beta);
        fmpz_mat_sub(Cvec_1, Cvec_1, Cvec_2);
        fmpz_mat_mul(coef, Lvec, Cvec_1);
        fmpz_addmul_ui(eps, fmpz_mat_entry(coef, 0, 0), 2);

        /* Convert alphabeta mod 2 to ulong */
        for (i = 0; i < 2 * g; i++)
        {
            ch[ab - start] = ch[ab - start] << 1;
            ch[ab - start] += fmpz_tstbit(fmpz_mat_entry(alphabeta, i, 0), 0);
        }
        /* Adjust sign of eps and reduce mod 8 */
        for (i = 0; i < g; i++)
        {
            if (fmpz_mod_ui(x, fmpz_mat_entry(alphabeta, i, 0), 2) == 1
                && fmpz_mod_ui(x, fmpz_mat_entry(alphabeta, i + g, 0), 4) > 1)
            {
                fmpz_add_ui(eps, eps, 4);
            }
        }
        e[ab - start] = fmpz_mod_ui(eps, eps, 8);
    }

    fmpz_mat_window_clear(a);
    fmpz_mat_window_clear(b);
    fmpz_mat_window_clear(c);
    fmpz_mat_window_clear(d);
    fmpz_mat_clear(cdt);
    fmpz_mat_clear(abt);
    fmpz_mat_clear(mat_tp);
    fmpz_mat_clear(diags);
    fmpz_mat_clear(alphabeta);
    fmpz_mat_window_clear(alpha);
    fmpz_mat_window_clear(beta);
    fmpz_mat_clear(Cvec_1);
    fmpz_mat_clear(Cvec_2);
    fmpz_mat_clear(Lvec);
    fmpz_mat_clear(coef);
    fmpz_clear(eps);
    fmpz_clear(x);
}
