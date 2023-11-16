/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly_factor.h"
#include "fmpz_lll.h"
#include "arb_fmpz_poly.h"
#include "qqbar.h"

#define QQBAR_GUESS_TRY_SMALLER 1
#define QQBAR_GUESS_CHECK 2

int
qqbar_guess(qqbar_t res, const acb_t z, slong max_deg, slong max_bits, int flags, slong prec)
{
    acb_ptr zpow;
    int found;
    slong i, j, fac_bits, prec2;
    fmpz_poly_t poly;
    fmpz_poly_factor_t fac;
    acb_t z2;
    mag_t rad;

    if (!acb_is_finite(z))
        return 0;

    if (max_deg > 8 && (flags & QQBAR_GUESS_TRY_SMALLER))
    {
        if (qqbar_guess(res, z, max_deg / 4, max_bits, flags, prec))
            return 1;
    }

    found = 0;

    fmpz_poly_init2(poly, max_deg + 1);
    fmpz_poly_factor_init(fac);
    acb_init(z2);
    mag_init(rad);
    zpow = _acb_vec_init(max_deg + 1);

    _acb_vec_set_powers(zpow, z, max_deg + 1, prec);

    _fmpz_poly_set_length(poly, max_deg + 1);
    if (_qqbar_acb_lindep(poly->coeffs, zpow, max_deg + 1, 1, prec))
    {
        _fmpz_poly_normalise(poly);

        fmpz_poly_factor(fac, poly);

        for (i = 0; i < fac->num; i++)
        {
            fac_bits = fmpz_poly_max_bits(fac->p + i);
            fac_bits = FLINT_ABS(fac_bits);

            if (fac_bits <= max_bits)
            {
                slong deg;
                qqbar_ptr roots;

                /* Rule out p(x) != 0 */
                arb_fmpz_poly_evaluate_acb(z2, fac->p + i, z, prec);

                if (acb_contains_zero(z2))
                {
                    /* Try to use the original interval, to avoid computing the polynomial roots... */
                    if (acb_rel_accuracy_bits(z) >= QQBAR_DEFAULT_PREC - 3)
                    {
                        for (prec2 = QQBAR_DEFAULT_PREC / 2; prec2 < 2 * prec; prec2 *= 2)
                        {
                            acb_set(z2, z);
                            acb_get_mag(rad, z);
                            mag_mul_2exp_si(rad, rad, -prec2);
                            acb_add_error_mag(z2, rad);

                            if (_qqbar_validate_existence_uniqueness(z2, fac->p + i, z2, 2 * prec2))
                            {
                                fmpz_poly_set(QQBAR_POLY(res), fac->p + i);
                                acb_set(QQBAR_ENCLOSURE(res), z2);
                                found = 1;
                                break;
                            }
                        }
                    }

                    deg = fmpz_poly_degree(fac->p + i);
                    roots = _qqbar_vec_init(deg);

                    qqbar_roots_fmpz_poly(roots, fac->p + i, QQBAR_ROOTS_IRREDUCIBLE);

                    for (j = 0; j < deg; j++)
                    {
                        qqbar_get_acb(z2, roots + j, prec);
                        if (acb_overlaps(z, z2))
                        {
                            qqbar_swap(res, roots + j);
                            found = 1;
                            break;
                        }
                    }

                    _qqbar_vec_clear(roots, deg);

                    if (found)
                        break;
                }
            }
        }
    }

    fmpz_poly_clear(poly);
    fmpz_poly_factor_clear(fac);
    _acb_vec_clear(zpow, max_deg + 1);
    acb_clear(z2);
    mag_clear(rad);

    return found;
}

