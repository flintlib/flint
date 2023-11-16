/*
    Copyright (C) 2015 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "arb_poly.h"
#include "acb_poly.h"
#include "acb_hypgeom.h"

void
_acb_poly_reciprocal_majorant(arb_ptr res, acb_srcptr vec, slong len, slong prec)
{
    slong i;

    for (i = 0; i < len; i++)
    {
        if (i == 0)
        {
            acb_get_abs_lbound_arf(arb_midref(res + i), vec + i, prec);
            mag_zero(arb_radref(res + i));
        }
        else
        {
            acb_get_abs_ubound_arf(arb_midref(res + i), vec + i, prec);
            arf_neg(arb_midref(res + i), arb_midref(res + i));
            mag_zero(arb_radref(res + i));
        }
    }
}

void
acb_poly_reciprocal_majorant(arb_poly_t res, const acb_poly_t poly, slong prec)
{
    arb_poly_fit_length(res, poly->length);
    _acb_poly_reciprocal_majorant(res->coeffs, poly->coeffs, poly->length, prec);
    _arb_poly_set_length(res, poly->length);
}


/* F = 1 + U + U^2 + ... = 1/(1-U) assuming that U[0] is positive;
   indeterminate if not convergent */
static void
arb_poly_geometric_sum(arb_poly_t F, const arb_poly_t U, slong len, slong prec)
{
    if (U->length == 0)
    {
        arb_poly_one(F);
        return;
    }

    arb_poly_add_si(F, U, -1, prec);
    arb_poly_neg(F, F);

    if (F->length > 0 && arb_is_positive(F->coeffs))
    {
        arb_poly_inv_series(F, F, len, prec);
    }
    else
    {
        arb_poly_fit_length(F, len);
        _arb_vec_indeterminate(F->coeffs, len);
        _arb_poly_set_length(F,  len);
    }
}

/* F = 1 + U + U^2 + U^3 + ... = 1/(1-U)

   U = product of (1 + |A-B|/(|B[0] - |B[1:]|)
       product of (1 / (|B[0] - |B[1:]|))
       * |Z|
*/
void
acb_hypgeom_pfq_series_bound_factor(arb_poly_t F,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z,
    slong n, slong len, slong prec)
{
    slong i;

    arb_poly_t T, U, V;
    acb_poly_t BN, AB;

    /* not convergent */
    if (p > q)
    {
        arb_poly_fit_length(F, len);
        _arb_vec_indeterminate(F->coeffs, len);
        _arb_poly_set_length(F, len);
        return;
    }

    arb_poly_init(T);
    arb_poly_init(U);
    arb_poly_init(V);

    acb_poly_init(BN);
    acb_poly_init(AB);

    acb_poly_majorant(U, z, prec);

    for (i = 0; i < q; i++)
    {
        acb_poly_add_si(BN, b + i, n, prec);

        if (acb_poly_length(BN) != 0 &&
                arb_is_positive(acb_realref(BN->coeffs)))
        {
            if (i < p)
            {
                /* 1 + |a-b|/reciprocal_majorant(b + n) */
                acb_poly_sub(AB, a + i, b + i, prec);
                acb_poly_majorant(T, AB, prec);
                acb_poly_reciprocal_majorant(V, BN, prec);
                arb_poly_div_series(T, T, V, len, prec);
                arb_poly_add_si(T, T, 1, prec);
                arb_poly_mullow(U, U, T, len, prec);
            }
            else
            {
                acb_poly_reciprocal_majorant(T, BN, prec);
                arb_poly_div_series(U, U, T, len, prec);
            }
        }
        else
        {
            arb_poly_fit_length(U, len);
            _arb_vec_indeterminate(U->coeffs, len);
            _arb_poly_set_length(U,  len);
            break;
        }
    }

    /* F = 1/(1-U) */
    arb_poly_geometric_sum(F, U, len, prec);

    arb_poly_clear(T);
    arb_poly_clear(U);
    arb_poly_clear(V);

    acb_poly_clear(BN);
    acb_poly_clear(AB);
}

void
acb_hypgeom_pfq_series_direct(acb_poly_t res,
    const acb_poly_struct * a, slong p,
    const acb_poly_struct * b, slong q,
    const acb_poly_t z, int regularized,
    slong n, slong len, slong prec)
{
    acb_poly_t s, t, err;
    arb_poly_t C, T;
    slong i;
    int is_real;
    int terminating;

    /* default algorithm to choose number of terms */
    if (n < 0)
    {
        n = acb_hypgeom_pfq_series_choose_n(a, p, b, q, z, len, prec);
    }

    terminating = 0;

    /* check if it terminates due to a root of the numerator */
    for (i = 0; i < p; i++)
    {
        if (acb_poly_length(a + i) == 0 && n > 0)
        {
            terminating = 1;
        }
        else if (acb_poly_length(a + i) == 1)
        {
            acb_srcptr c = acb_poly_get_coeff_ptr(a + i, 0);

            if (acb_is_int(c) && arb_is_negative(acb_realref(c)) &&
                arf_cmpabs_ui(arb_midref(acb_realref(c)), n) < 0)
            {
                terminating = 1;
            }
        }
    }

    /* check if it terminates (to order n) due to z */
    /* the following tests could be made stronger... */
    if (z->length == 0 && n >= 1)
    {
        terminating = 1;
    }
    else if (!terminating && z->length > 0 && acb_is_zero(z->coeffs) && n >= len)
    {
        if (regularized)
        {
            terminating = 1;
        }
        else
        {
            terminating = 1;

            for (i = 0; i < q; i++)
            {
                acb_srcptr c = acb_poly_get_coeff_ptr(b + i, 0);

                if (!arb_is_positive(acb_realref(c)) && acb_contains_int(c))
                    terminating = 0;
            }
        }
    }

    acb_poly_init(s);
    acb_poly_init(t);
    acb_poly_init(err);
    arb_poly_init(C);
    arb_poly_init(T);

    acb_hypgeom_pfq_series_sum(s, t, a, p, b, q, z, regularized, n, len, prec);

    if (!terminating)
    {
        is_real = acb_poly_is_real(z);
        for (i = 0; i < p; i++)
            is_real = is_real && acb_poly_is_real(a + i);
        for (i = 0; i < q; i++)
            is_real = is_real && acb_poly_is_real(b + i);

        acb_poly_majorant(T, t, MAG_BITS);
        acb_hypgeom_pfq_series_bound_factor(C, a, p, b, q, z, n, len, MAG_BITS);

        if (!_arb_vec_is_finite(T->coeffs, T->length) ||
            !_arb_vec_is_finite(C->coeffs, C->length))
        {
            arb_poly_fit_length(T, len);
            _arb_vec_indeterminate(T->coeffs, len);
            _arb_poly_set_length(T, len);
        }
        else
        {
            arb_poly_mullow(T, T, C, len, MAG_BITS);
        }

        /* create polynomial of errors */
        acb_poly_fit_length(err, len);

        for (i = 0; i < FLINT_MIN(len, T->length); i++)
        {
            arb_add_error(acb_realref(err->coeffs + i), T->coeffs + i);
            if (!is_real)
                arb_add_error(acb_imagref(err->coeffs + i), T->coeffs + i);
        }

        _acb_poly_set_length(err, len);
        _acb_poly_normalise(err);

        acb_poly_add(s, s, err, prec);
    }

    acb_poly_set(res, s);

    acb_poly_clear(s);
    acb_poly_clear(t);
    acb_poly_clear(err);
    arb_poly_clear(C);
    arb_poly_clear(T);
}

