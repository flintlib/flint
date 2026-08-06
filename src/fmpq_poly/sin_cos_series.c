/*
    Copyright (C) 2016, 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

void
_fmpq_poly_sin_cos_series_basecase(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    fmpz_t v, Aden_copy;
    fmpz * Aprime;
    fmpz * S_alloc = NULL, * C_alloc = NULL;
    fmpz_t Sden_alloc, Cden_alloc;
    slong k, l;

    Alen = FLINT_MIN(Alen, n);

    if (Alen == 1 || n == 1)
    {
        if (S != NULL)
        {
            fmpz_zero(S);
            _fmpz_vec_zero(S, n);
            fmpz_one(Sden);
        }

        if (C != NULL)
        {
            fmpz_one(C);
            _fmpz_vec_zero(C + 1, n - 1);
            fmpz_one(Cden);
        }
        return;
    }

    if (S == NULL)
    {
        fmpz_init(Sden_alloc);
        S_alloc = _fmpz_vec_init(n);
        S = S_alloc;
        Sden = Sden_alloc;
    }
    if (C == NULL)
    {
        fmpz_init(Cden_alloc);
        C_alloc = _fmpz_vec_init(n);
        C = C_alloc;
        Cden = Cden_alloc;
    }

    Aprime = _fmpz_vec_init(Alen - 1);
    _fmpz_poly_derivative(Aprime, A, Alen);

    fmpz_init(v);
    fmpz_init(Aden_copy);

    /* Save Aden in case of aliasing with Sden/Cden. */
    fmpz_set(Aden_copy, Aden);

    fmpz_fac_ui(Sden, n - 1);
    fmpz_pow_ui(v, Aden_copy, n - 1);
    fmpz_mul(Sden, Sden, v);
    fmpz_set(C, Sden);
    fmpz_set(Cden, Sden);
    fmpz_zero(S);

    {
        fmpz_mul(C + 1, Aprime, S);
        fmpz_mul(S + 1, Aprime, C);
        fmpz_divexact(C + 1, C + 1, Aden_copy);
        fmpz_divexact(S + 1, S + 1, Aden_copy);
    }

    for (k = 2; k < n; k++)
    {
        l = FLINT_MIN(Alen - 1, k);
        _fmpz_vec_dot_general(C + k, NULL, 1, Aprime, S + k - l, 1, l);
        _fmpz_vec_dot_general(S + k, NULL, 0, Aprime, C + k - l, 1, l);

        fmpz_mul_ui(v, Aden_copy, k);
        fmpz_divexact(C + k, C + k, v);
        fmpz_divexact(S + k, S + k, v);
    }

    _fmpz_vec_clear(Aprime, Alen - 1);

    fmpz_clear(v);
    fmpz_clear(Aden_copy);

    if (S_alloc == NULL)
    {
        _fmpq_poly_canonicalise(S, Sden, n);
    }
    else
    {
        _fmpz_vec_clear(S_alloc, n);
        fmpz_clear(Sden_alloc);
    }

    if (C_alloc == NULL)
    {
        _fmpq_poly_canonicalise(C, Cden, n);
    }
    else
    {
        _fmpz_vec_clear(C_alloc, n);
        fmpz_clear(Cden_alloc);
    }
}

static void
_fmpq_poly_neg(fmpz * r, fmpz_t rden, const fmpz * P, const fmpz_t Pden, slong n)
{
    _fmpz_vec_neg(r, P, n);
    fmpz_set(rden, Pden);
}

static void
_fmpq_poly_set(fmpz * r, fmpz_t rden, const fmpz * P, const fmpz_t Pden, slong n)
{
    _fmpz_vec_set(r, P, n);
    fmpz_set(rden, Pden);
}

static void
MULLOW(fmpz * z, fmpz_t zden, const fmpz * x, const fmpz_t xden, slong xn, const fmpz * y, const fmpz_t yden, slong yn, slong n)
{
    FLINT_ASSERT(xn + yn - 1 >= n);

    if (xn >= yn)
        _fmpz_poly_mullow(z, x, xn, y, yn, n);
    else
        _fmpz_poly_mullow(z, y, yn, x, xn, n);

    fmpz_mul(zden, xden, yden);
    _fmpq_poly_canonicalise(z, zden, n);
}

static void
MULMID(fmpz * z, fmpz_t zden, const fmpz * x, const fmpz_t xden, slong xn, const fmpz * y, const fmpz_t yden, slong yn, slong nlo, slong nhi)
{
    FLINT_ASSERT(xn + yn - 1 >= nhi);

    _fmpz_poly_mulmid(z, x, xn, y, yn, nlo, nhi);
    fmpz_mul(zden, xden, yden);
    _fmpq_poly_canonicalise(z, zden, nhi - nlo);
}

/* Assuming that the low m coefficients of poly have denominator
   den in canonical form and that the high n - m coefficients have
   denominator high_den in canonical form, combine the high and lower
   parts and put the polynomial in canonical form {poly, den, n}. */
static void
CONCATENATE(fmpz * poly, fmpz_t den, const fmpz_t high_den, slong m, slong n)
{
    fmpz_t gcd, d1, d2;

    fmpz_init(gcd);
    fmpz_init(d1);
    fmpz_init(d2);

    fmpz_gcd(gcd, den, high_den);
    fmpz_divexact(d1, high_den, gcd);
    fmpz_divexact(d2, den, gcd);

    _fmpz_vec_scalar_mul_fmpz(poly, poly, m, d1);
    _fmpz_vec_scalar_mul_fmpz(poly + m, poly + m, n - m, d2);

    fmpz_mul(den, d2, high_den);

    fmpz_clear(gcd);
    fmpz_clear(d1);
    fmpz_clear(d2);
}

/* Adapted from _gr_poly_sin_cos_series_newton */
void
_fmpq_poly_sin_cos_series_newton(
    fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden,
    const fmpz * h, const fmpz_t hden, slong hlen, slong cutoff, slong n)
{
    slong a[FLINT_BITS];
    slong original_n, i, m, l, r, nm, Alen;
    fmpz * hprime, * P, * Q, * t, * u;
    fmpz_t hprimeden, hdenin, Pden, Qden, tden, uden, high_den;
    fmpz * S_scratch = NULL, * C_scratch = NULL;
    fmpz_t Sden_scratch, Cden_scratch;
    int want_S, want_C;

    want_S = (S != NULL);
    want_C = (C != NULL);

    /* Allocate scratch vectors for the unwanted output(s). */
    if (!want_S)
    {
        fmpz_init(Sden_scratch);
        S_scratch = _fmpz_vec_init(n);
        S = S_scratch;
        Sden = Sden_scratch;
    }
    if (!want_C)
    {
        fmpz_init(Cden_scratch);
        C_scratch = _fmpz_vec_init(n);
        C = C_scratch;
        Cden = Cden_scratch;
    }

    original_n = n;
    hlen = FLINT_MIN(hlen, n);

    cutoff = FLINT_MAX(cutoff, 2);

    for (i = 1; (WORD(1) << i) < n; i++);
    a[i = 0] = n;
    while (n >= cutoff || i == 0)
        a[++i] = (n = (n + 1) / 2);

    hprime = _fmpz_vec_init(hlen - 1);
    P      = _fmpz_vec_init(original_n);
    Q      = _fmpz_vec_init(original_n);
    t      = _fmpz_vec_init(original_n);
    u      = _fmpz_vec_init(original_n);

    fmpz_init(hprimeden);
    fmpz_init(hdenin);
    fmpz_init(Pden);
    fmpz_init(Qden);
    fmpz_init(tden);
    fmpz_init(uden);
    fmpz_init(high_den);

    fmpz_set(hprimeden, hden);
    fmpz_set(hdenin, hden);
    _fmpz_poly_derivative(hprime, h, hlen);
    _fmpq_poly_canonicalise(hprime, hprimeden, FLINT_MIN(n, hlen) - 1);

    _fmpq_poly_sin_cos_series_basecase(S, Sden, C, Cden, h, hden, hlen, n);

    for (i--; i >= 0; i--)
    {
        m = n;
        n = a[i];
        nm = n - m;

        /* Karatsuba seems to be slower only in the basecase range. */
        int use_karatsuba = 1;

        l    = FLINT_MIN(hlen, n) - 1;
        r    = FLINT_MIN(l + m - 1, n - 1);
        Alen = r - m + 1;

        /* Extend hprime if more coefficients become relevant this step. */
        if (l > m - 1)
        {
            fmpz_set(uden, hdenin);
            _fmpq_poly_canonicalise(hprime + m - 1, uden, l - m + 1);
            CONCATENATE(hprime, hprimeden, uden, m - 1, l);
        }

        /* Remark: some of the canonicalisations between operations could
           be avoided, particularly in the Karatsuba blocks, but the
           speedup is small/inconsistent. */

        /* Compute B = [h'S]_{m-1..r-1} -> t/tden, A = [h'C]_{m-1..r-1} -> u/uden. */
        MULMID(t, tden, hprime, hprimeden, l, S, Sden, m, m - 1, r);
        MULMID(u, uden, hprime, hprimeden, l, C, Cden, m, m - 1, r);

        /* Compute P = A*S_low - B*C_low  and  Q = A*C_low + B*S_low */
        if (!use_karatsuba)
        {
            MULLOW(P, Pden, u, uden, Alen, S, Sden, nm, nm);
            MULLOW(Q, Qden, u, uden, Alen, C, Cden, nm, nm);
            MULLOW(u, uden, t, tden, Alen, C, Cden, nm, nm);
            _fmpq_poly_sub(P, Pden, P, Pden, nm, u, uden, nm);
            MULLOW(u, uden, t, tden, Alen, S, Sden, nm, nm);
            _fmpq_poly_add(Q, Qden, Q, Qden, nm, u, uden, nm);
        }
        else
        {
            MULLOW(P, Pden, t, tden, Alen, C, Cden, nm, nm);
            _fmpq_poly_neg(P, Pden, P, Pden, nm);
            MULLOW(Q, Qden, u, uden, Alen, S, Sden, nm, nm);
            _fmpq_poly_neg(Q, Qden, Q, Qden, nm);
            _fmpq_poly_sub(u, uden, u, uden, Alen, t, tden, Alen);
            _fmpq_poly_sub(t, tden, C, Cden, nm, S, Sden, nm);
            MULLOW(C + m, high_den, u, uden, Alen, t, tden, nm, nm);
            _fmpq_poly_sub(C + m, high_den, C + m, high_den, nm, P, Pden, nm);
            _fmpq_poly_sub(C + m, high_den, C + m, high_den, nm, Q, Qden, nm);
            _fmpq_poly_sub(P, Pden, P, Pden, nm, Q, Qden, nm);
            _fmpq_poly_set(Q, Qden, C + m, high_den, nm);
        }

        /* Integrate P and Q */
        _fmpq_poly_integral_offset(P, Pden, P, Pden, nm, m);
        _fmpq_poly_integral_offset(Q, Qden, Q, Qden, nm, m);

        /* Update C and S */
        if (use_karatsuba && (i != 0 || (want_S && want_C)))
        {
            MULLOW(C + m, high_den, C, Cden, nm, P, Pden, nm, nm);
            MULLOW(S + m, uden, S, Sden, nm, Q, Qden, nm, nm);
            _fmpq_poly_add(P, Pden, P, Pden, nm, Q, Qden, nm);
            _fmpq_poly_add(Q, Qden, C, Cden, nm, S, Sden, nm);
            MULLOW(t, tden, Q, Qden, nm, P, Pden, nm, nm);
            _fmpq_poly_sub(t, tden, t, tden, nm, C + m, high_den, nm);
            _fmpq_poly_sub(t, tden, t, tden, nm, S + m, uden,     nm);
            _fmpq_poly_sub(C + m, high_den,
                           C + m, high_den, nm, S + m, uden, nm);
            _fmpz_vec_set(S + m, t, nm);
            fmpz_set(uden, tden);
            CONCATENATE(C, Cden, high_den, m, n);
            CONCATENATE(S, Sden, uden,     m, n);
        }
        else
        {
            if (want_C || i != 0)
            {
                MULLOW(t, tden, C, Cden, nm, P, Pden, nm, nm);
                MULLOW(u, uden, S, Sden, nm, Q, Qden, nm, nm);
                _fmpq_poly_sub(C + m, high_den, t, tden, nm, u, uden, nm);
                CONCATENATE(C, Cden, high_den, m, n);
            }

            if (want_S || i != 0)
            {
                MULLOW(t, tden, S, Sden, nm, P, Pden, nm, nm);
                MULLOW(u, uden, C, Cden, nm, Q, Qden, nm, nm);
                _fmpq_poly_add(S + m, high_den, t, tden, nm, u, uden, nm);
                CONCATENATE(S, Sden, high_den, m, n);
            }
        }
    }

    _fmpz_vec_clear(hprime, hlen - 1);
    _fmpz_vec_clear(P, original_n);
    _fmpz_vec_clear(Q, original_n);
    _fmpz_vec_clear(t, original_n);
    _fmpz_vec_clear(u, original_n);

    fmpz_clear(hprimeden);
    fmpz_clear(hdenin);
    fmpz_clear(Pden);
    fmpz_clear(Qden);
    fmpz_clear(tden);
    fmpz_clear(uden);
    fmpz_clear(high_den);

    if (S_scratch != NULL)
    {
        _fmpz_vec_clear(S_scratch, original_n);
        fmpz_clear(Sden_scratch);
    }
    if (C_scratch != NULL)
    {
        _fmpz_vec_clear(C_scratch, original_n);
        fmpz_clear(Cden_scratch);
    }
}

void
_fmpq_poly_sin_cos_series(fmpz * S, fmpz_t Sden,
    fmpz * C, fmpz_t Cden, const fmpz * A, const fmpz_t Aden,
    slong Alen, slong n)
{
    slong cutoff = 20;

    if (Alen < cutoff || n < cutoff)
        _fmpq_poly_sin_cos_series_basecase(S, Sden, C, Cden, A, Aden, Alen, n);
    else
        _fmpq_poly_sin_cos_series_newton(S, Sden, C, Cden, A, Aden, Alen, cutoff, n);
}

void
fmpq_poly_sin_cos_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_zero(res2);
        return;
    }

    if (fmpq_poly_is_zero(poly) || n == 1)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_one(res2);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sin_cos_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res1, n);
    fmpq_poly_fit_length(res2, n);
    _fmpq_poly_sin_cos_series(res1->coeffs, res1->den,
        res2->coeffs, res2->den, poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res1, n);
    _fmpq_poly_normalise(res1);
    _fmpq_poly_set_length(res2, n);
    _fmpq_poly_normalise(res2);
}

void
_fmpq_poly_sin_series(fmpz * g, fmpz_t gden,
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    _fmpq_poly_sin_cos_series(g, gden, NULL, NULL, h, hden, hlen, n);
}

void fmpq_poly_sin_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (poly->length && !fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_sin_series). Constant term != 0.\n");
    }

    if (poly->length == 0 || n < 2)
    {
        fmpq_poly_zero(res);
        return;
    }

    fmpq_poly_fit_length(res, n);
    _fmpq_poly_sin_series(res->coeffs, res->den,
        poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}

void
_fmpq_poly_cos_series(fmpz * g, fmpz_t gden,
                       const fmpz * h, const fmpz_t hden, slong hlen, slong n)
{
    _fmpq_poly_sin_cos_series(NULL, NULL, g, gden, h, hden, hlen, n);
}

void fmpq_poly_cos_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res);
        return;
    }

    if (poly->length == 0 || n == 1)
    {
        fmpq_poly_one(res);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_cos_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res, n);
    _fmpq_poly_cos_series(res->coeffs, res->den,
        poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}

