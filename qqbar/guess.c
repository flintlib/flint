/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "flint/fmpz_lll.h"
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
                    roots = qqbar_vec_init(deg);

                    qqbar_roots_fmpz_poly(roots, fac->p + i, QQBAR_ROOTS_IRREDUCIBLE);

                    for (j = 0; j < deg; j++)
                    {
                        qqbar_get_acb(z2, roots + j, prec);
                        if (acb_overlaps(z, z2))
                        {
                            found = 1;
                            break;
                        }
                    }

                    qqbar_vec_clear(roots, deg);

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


/*

    @staticmethod
    def guess(z, deg=None, verbose=False, try_smaller=True, check=True):
        """
        Use LLL to try to find an algebraic number matching the real or
        complex enclosure *z* (given as a *flint.arb* or *flint.acb*).
        Returns *None* if unsuccessful.

        To determine the number of bits to use in the LLL matrix, the
        algorithm accounts for the accuracy of the input as well as the
        current Arb working precision (*flint.ctx.prec*).

        The maximum degree is determined by the *deg* parameter.
        By default the degree is set to roughly the square root of the working
        precision.
        If *try_smaller* is set, then the guessing is attempted recursively
        with lower degree bounds to speed up detection of lower-degree
        numbers.

        If *check* is True (default), the algorithm verifies that the
        result is contained in the complex interval *z*. (At high precision,
        this reduces the likelihood of spurious results.)

        """
        z = acb(z)
        if not z.is_finite():
            return None
        if deg is None:
            deg = max(1, int(ctx.prec ** 0.5))
        if try_smaller and deg > 8:
            g = alg.guess(z, deg // 4, verbose=verbose, try_smaller=True, check=check)
            if g is not None:
                return g
        # todo: early detect exact (dyadic) real and imaginary parts?
        # (avoid reliance on enclosure(pretty=True)
        prec = ctx.prec
        z = acb(z)
        zpow = [z**i for i in range(deg+1)]
        magn = max(abs(t) for t in zpow)
        if magn >= arb(2)**prec:
            return None   # too large
        magnrad = max(abs(t.rad()) for t in zpow)
        if magnrad >= 0.125:
            return None   # too imprecise
        magnrad = max(2**(max(1,prec//20)) * magnrad, arb(2)**(-2*prec))
        scale = 1 / magnrad
        nonreal = z.imag != 0
        A = fmpz_mat(deg + 1, deg + 2 + nonreal)
        for i in range(deg+1):
            A[i,i] = 1
        for i in range(deg+1):
            v = zpow[i] * scale
            # fixme: don't use unique_fmpz
            A[i,deg + 1] = v.real.mid().floor().unique_fmpz()
            if nonreal:
                A[i,deg + 2] = v.imag.mid().floor().unique_fmpz()
        A = A.lll()
        poly = fmpz_poly([A[0,i] for i in range(deg+1)])
        if verbose:
            print("LLL-reduced matrix:")
            print(A)
            print("poly(z) =", poly(z))
        if not check:
            candidates = alg.polynomial_roots(poly) or [(alg(0), 1)]
            return min([cand for (cand, mult) in candidates], key=lambda x: abs(x.enclosure()-z).mid())
        if not poly(z).contains(0):
            return None
        orig = ctx.prec
        try:
            candidates = alg.polynomial_roots(poly)
            while ctx.prec <= 16 * orig:
                for cand, mult in candidates:
                    r = cand.enclosure(pretty=True)
                    if z.contains(r):
                        return cand
                ctx.prec *= 2
        finally:
            ctx.prec = orig
        return None

*/
