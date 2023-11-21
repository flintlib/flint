/*
    Copyright (C) 2012 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "arb.h"

/* With parameter n, the error is bounded by 3/(3+sqrt(8))^n */
#define ERROR_A 1.5849625007211561815 /* log2(3) */
#define ERROR_B 2.5431066063272239453 /* log2(3+sqrt(8)) */

void mag_borwein_error(mag_t err, slong n);

void
arb_zeta_ui_vec_borwein(arb_ptr z, ulong start, slong num, ulong step, slong prec)
{
    slong j, k, s, n, wp;
    fmpz_t c, d, t, u;
    fmpz * zeta;
    mag_t err;

    if (num < 1)
        return;

    wp = prec + FLINT_BIT_COUNT(prec);
    n = wp / 2.5431066063272239453 + 1;

    fmpz_init(c);
    fmpz_init(d);
    fmpz_init(t);
    fmpz_init(u);
    zeta = _fmpz_vec_init(num);

    fmpz_set_ui(c, 1);
    fmpz_mul_2exp(c, c, 2 * n - 1);
    fmpz_set(d, c);

    for (k = n; k > 0; k--)
    {
        /* divide by first k^s */
        fmpz_ui_pow_ui(u, k, start);
        fmpz_tdiv_q(t, d, u);
        if (k % 2 == 0)
            fmpz_neg(t, t);
        fmpz_add(zeta, zeta, t);

        /* remaining k^s */
        fmpz_ui_pow_ui(u, k, step);
        for (j = 1; j < num; j++)
        {
            fmpz_tdiv_q(t, t, u);
            fmpz_add(zeta + j, zeta + j, t);
        }

        /* hypergeometric recurrence */
        fmpz_mul2_uiui(c, c, k, 2 * k - 1);
        fmpz_divexact2_uiui(c, c, 2 * (n - k + 1), n + k - 1);
        fmpz_add(d, d, c);
    }

    mag_init(err);
    mag_borwein_error(err, n);

    for (k = 0; k < num; k++)
    {
        arb_ptr x = z + k;
        s = start + step * k;

        arb_set_fmpz(x, zeta + k);
        /* the error in each term in the main loop is < 2 */
        mag_set_ui(arb_radref(x), 2 * n);
        arb_div_fmpz(x, x, d, wp);

        /* mathematical error for eta(s), bounded by 3/(3+sqrt(8))^n */
        mag_add(arb_radref(x), arb_radref(x), err);

        /* convert from eta(s) to zeta(s) */
        arb_div_2expm1_ui(x, x, s - 1, wp);
        arb_mul_2exp_si(x, x, s - 1);
    }

    mag_clear(err);

    fmpz_clear(c);
    fmpz_clear(d);
    fmpz_clear(t);
    fmpz_clear(u);
    _fmpz_vec_clear(zeta, num);
}
