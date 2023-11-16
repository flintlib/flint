/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq_poly.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

/* Algorithm based on Cohen section 4.5.1 */

/* C = A o P  mod B */
/* todo: rem P by B initially if larger */
/* todo: brent-kung or something more clever */
void
fmpq_poly_compose_fmpz_poly_mod_fmpz_poly(fmpq_poly_t C, const fmpz_poly_t A, const fmpq_poly_t P, const fmpz_poly_t B)
{
    slong i, m;
    fmpq_poly_t B2;
    fmpq_poly_init(B2);

    fmpq_poly_set_fmpz_poly(B2, B);

    m = fmpz_poly_degree(A);

    fmpq_poly_set_fmpz(C, A->coeffs + m);

    for (i = m - 1; i >= 0; i--)
    {
        fmpq_poly_mul(C, C, P);
        fmpq_poly_add_fmpz(C, C, A->coeffs + i);
        fmpq_poly_rem(C, C, B2);
    }

    fmpq_poly_clear(B2);
}

int
qqbar_equal_fmpq_poly_val(const qqbar_t x, const fmpq_poly_t f, const qqbar_t y)
{
    slong prec;
    acb_t z1, z2, z3;
    fmpq_poly_t C;
    int res;

    /* todo: should this be used for f->length == 2 too? */
    if (f->length <= 1 || qqbar_degree(y) == 1)
    {
        qqbar_t v;
        qqbar_init(v);
        qqbar_evaluate_fmpq_poly(v, f, y);
        res = qqbar_equal(v, x);
        qqbar_clear(v);
        return res;
    }

    acb_init(z1);
    acb_init(z2);
    acb_init(z3);
    fmpq_poly_init(C);

    acb_set(z1, QQBAR_ENCLOSURE(x));
    acb_set(z2, QQBAR_ENCLOSURE(y));

    res = 0;
    for (prec = QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        _qqbar_enclosure_raw(z1, QQBAR_POLY(x), z1, prec);
        _qqbar_enclosure_raw(z2, QQBAR_POLY(y), z2, prec);

        _arb_fmpz_poly_evaluate_acb(z3, f->coeffs, f->length, z2, 2 * prec);
        acb_div_fmpz(z3, z3, f->den, 2 * prec);

        if (!acb_overlaps(z1, z3))
        {
            res = 0;
            break;
        }

        if (prec == QQBAR_DEFAULT_PREC / 2)
        {
            fmpq_poly_compose_fmpz_poly_mod_fmpz_poly(C, QQBAR_POLY(x), f, QQBAR_POLY(y));
            if (!fmpq_poly_is_zero(C))
            {
                res = 0;
                break;
            }
        }

        acb_union(z3, z1, z3, prec);

        if (_qqbar_validate_uniqueness(z3, QQBAR_POLY(x), z3, 2 * prec))
        {
            res = 1;
            break;
        }
    }

    acb_clear(z1);
    acb_clear(z2);
    acb_clear(z3);
    fmpq_poly_clear(C);

    return res;
}

