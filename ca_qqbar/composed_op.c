/*
    Copyright (C) 2020 Fredrik Johansson

    This file is part of Calcium.

    Calcium is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "arb_fmpz_poly.h"
#include "ca_qqbar.h"

#define OP_ADD 0
#define OP_SUB 1
#define OP_MUL 2
#define OP_DIV 3

static void
fmpq_poly_hadamard_product(fmpq_poly_t res, const fmpq_poly_t poly1, const fmpq_poly_t poly2)
{
    slong i, len;

    len = FLINT_MIN(fmpq_poly_length(poly1), fmpq_poly_length(poly2));

    fmpq_poly_fit_length(res, len);

    for (i = 0; i < len; i++)
        fmpz_mul(res->coeffs + i, poly1->coeffs + i, poly2->coeffs + i);

    fmpz_mul(res->den, poly1->den, poly2->den);
    _fmpq_poly_set_length(res, len);
    _fmpq_poly_canonicalise(res->coeffs, res->den, len);
}

static void
fmpq_poly_borel_transform(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong i, len;

    len = fmpq_poly_length(poly);

    if (len <= 2)
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpz_t c;
        fmpz_init(c);

        fmpz_one(c);
        fmpq_poly_fit_length(res, len);

        for (i = len - 1; i >= 0; i--)
        {
            fmpz_mul(res->coeffs + i, poly->coeffs + i, c);
            if (i > 1)
                fmpz_mul_ui(c, c, i);
        }

        fmpz_mul(fmpq_poly_denref(res), fmpq_poly_denref(poly), c);

        _fmpq_poly_set_length(res, len);
        _fmpq_poly_canonicalise(res->coeffs, res->den, len);

        fmpz_clear(c);
    }
}

static void
fmpq_poly_inv_borel_transform(fmpq_poly_t res, const fmpq_poly_t poly)
{
    slong i, len;

    len = fmpq_poly_length(poly);

    if (len <= 2)
    {
        fmpq_poly_set(res, poly);
    }
    else
    {
        fmpz_t c;
        fmpz_init(c);

        fmpz_one(c);

        fmpq_poly_fit_length(res, len);
        fmpz_set(fmpq_poly_denref(res), fmpq_poly_denref(poly));
        fmpz_set(res->coeffs, poly->coeffs);
        fmpz_set(res->coeffs + 1, poly->coeffs + 1);

        for (i = 2; i < len; i++)
        {
            fmpz_mul_ui(c, c, i);
            fmpz_mul(res->coeffs + i, poly->coeffs + i, c);
        }

        _fmpq_poly_set_length(res, len);
        _fmpq_poly_canonicalise(res->coeffs, res->den, len);

        fmpz_clear(c);
    }
}

void
ca_qqbar_fmpz_poly_composed_op(fmpz_poly_t res, const fmpz_poly_t A, const fmpz_poly_t B, int op)
{
    slong d1, d2, n, i;
    fmpq_poly_t P1, P2, P1rev, P1drev, P2rev, P2drev;

    d1 = fmpz_poly_degree(A);
    d2 = fmpz_poly_degree(B);

    if (d1 <= 0 || d2 <= 0)
    {
        flint_printf("composed_op: inputs must not be constants\n");
        flint_abort();
    }

    n = d1 * d2 + 1;

    fmpq_poly_init(P1);
    fmpq_poly_init(P2);
    fmpq_poly_init(P1rev);
    fmpq_poly_init(P1drev);
    fmpq_poly_init(P2rev);
    fmpq_poly_init(P2drev);

    fmpq_poly_set_fmpz_poly(P1, A);
    fmpq_poly_set_fmpz_poly(P2, B);

    if (op == OP_DIV)
    {
        if (fmpz_is_zero(P2->coeffs))
        {
            flint_printf("composed_op: division by zero\n");
            flint_abort();
        }

        fmpq_poly_reverse(P2, P2, d2 + 1);
    }

    if (op == OP_SUB)
        for (i = 1; i <= d2; i += 2)
            fmpz_neg(P2->coeffs + i, P2->coeffs + i);

    fmpq_poly_reverse(P1rev, P1, d1 + 1);
    fmpq_poly_derivative(P1drev, P1);
    fmpq_poly_reverse(P1drev, P1drev, d1);

    fmpq_poly_reverse(P2rev, P2, d2 + 1);
    fmpq_poly_derivative(P2drev, P2);
    fmpq_poly_reverse(P2drev, P2drev, d2);

    fmpq_poly_div_series(P1, P1drev, P1rev, n);
    fmpq_poly_div_series(P2, P2drev, P2rev, n);

    if (op == OP_MUL || op == OP_DIV)
    {
        fmpq_poly_hadamard_product(P1, P1, P2);
        fmpq_poly_shift_right(P1, P1, 1);
        fmpq_poly_neg(P1, P1);
        fmpq_poly_integral(P1, P1);
    }
    else
    {
        fmpq_poly_borel_transform(P1, P1);
        fmpq_poly_borel_transform(P2, P2);
        fmpq_poly_mullow(P1, P1, P2, n);
        fmpq_poly_shift_right(P1, P1, 1);
        fmpq_poly_inv_borel_transform(P1, P1);
        fmpq_poly_neg(P1, P1);
        fmpq_poly_shift_left(P1, P1, 1);
    }

    fmpq_poly_exp_series(P1, P1, n);
    fmpq_poly_reverse(P1, P1, n);

    fmpq_poly_get_numerator(res, P1);

    fmpq_poly_clear(P1);
    fmpq_poly_clear(P2);
    fmpq_poly_clear(P1rev);
    fmpq_poly_clear(P1drev);
    fmpq_poly_clear(P2rev);
    fmpq_poly_clear(P2drev);
}

void
ca_qqbar_binary_op(ca_qqbar_t res, const ca_qqbar_t x, const ca_qqbar_t y, int op)
{
    slong i, prec, found;
    fmpz_poly_t H;
    fmpz_poly_factor_t fac;
    acb_t z1, z2, w, t;

    fmpz_poly_init(H);
    fmpz_poly_factor_init(fac);
    acb_init(z1);
    acb_init(z2);
    acb_init(w);
    acb_init(t);

    /* flint_printf("BEGIN COMPOSED OP %wd %wd %wd %wd\n",
        fmpz_poly_degree(CA_QQBAR_POLY(x)),
        fmpz_poly_max_bits(CA_QQBAR_POLY(x)),
        fmpz_poly_degree(CA_QQBAR_POLY(y)),
        fmpz_poly_max_bits(CA_QQBAR_POLY(y))); */
    ca_qqbar_fmpz_poly_composed_op(H, CA_QQBAR_POLY(x), CA_QQBAR_POLY(y), op);
    fmpz_poly_factor(fac, H);
    acb_set(z1, CA_QQBAR_ENCLOSURE(x));
    acb_set(z2, CA_QQBAR_ENCLOSURE(y));

/*
    ca_qqbar_print(x); printf("\n");
    ca_qqbar_print(y); printf("\n");
*/

    for (prec = CA_QQBAR_DEFAULT_PREC / 2; ; prec *= 2)
    {
        /* printf("binop %ld\n", prec); */

        _ca_qqbar_enclosure_raw(z1, CA_QQBAR_POLY(x), z1, prec);
        _ca_qqbar_enclosure_raw(z2, CA_QQBAR_POLY(y), z2, prec);

        /*
        acb_printd(z1, 30); printf("\n");
        acb_printd(z2, 30); printf("\n");
        */

        if (op == 0)
            acb_add(w, z1, z2, prec);
        else if (op == 1)
            acb_sub(w, z1, z2, prec);
        else if (op == 2)
            acb_mul(w, z1, z2, prec);
        else
            acb_div(w, z1, z2, prec);

        /* Look for potential roots -- we want exactly one */
        found = -1;
        for (i = 0; i < fac->num && found != -2; i++)
        {
            arb_fmpz_poly_evaluate_acb(t, fac->p + i, w, prec);
            if (acb_contains_zero(t))
            {
                if (found == -1)
                    found = i;
                else
                    found = -2;
            }
        }

        /* printf("found: %ld\n", found); */

        /* Check if the enclosure is good enough */
        if (found >= 0)
        {
            if (_ca_qqbar_validate_uniqueness(t, fac->p + found, w, 2 * prec))
            {
                fmpz_poly_set(CA_QQBAR_POLY(res), fac->p + found);
                acb_set(CA_QQBAR_ENCLOSURE(res), t);
                break;
            }
        }
    }

    fmpz_poly_clear(H);
    fmpz_poly_factor_clear(fac);
    acb_clear(z1);
    acb_clear(z2);
    acb_clear(w);
    acb_clear(t);
}

