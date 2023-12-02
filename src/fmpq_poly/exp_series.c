/*
    Copyright (C) 2010 Sebastian Pancratz
    Copyright (C) 2011, 2014, 2021 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"

static ulong _fmpz_gcd_big_small(const fmpz_t g, ulong h)
{
    __mpz_struct * z = COEFF_TO_PTR(*g);

    return n_gcd(mpn_mod_1(z->_mp_d, FLINT_ABS(z->_mp_size), h), h);
}

static ulong _fmpz_gcd_small(const fmpz_t g, ulong h)
{
    if (!COEFF_IS_MPZ(*g))
        return n_gcd(FLINT_ABS(*g), h);
    else
        return _fmpz_gcd_big_small(g, h);
}


/* Basecase algorithm, given a precomputed derivative of
   of the input series (Alen still refers to the length
   of the original series). */
void
_fmpq_poly_exp_series_basecase_deriv(fmpz * B, fmpz_t Bden,
    const fmpz * Aprime, const fmpz_t Aden, slong Alen, slong n)
{
    fmpz_t t, u;
    slong j, k;

    Alen = FLINT_MIN(Alen, n);

    fmpz_init(t);
    fmpz_init(u);

    fmpz_fac_ui(t, n - 1);
    fmpz_pow_ui(u, Aden, n - 1);

    fmpz_mul(Bden, t, u);
    fmpz_set(B, Bden);

    for (k = 1; k < n; k++)
    {
        fmpz_mul(t, Aprime, B + k - 1);

        for (j = 2; j < FLINT_MIN(Alen, k + 1); j++)
            fmpz_addmul(t, Aprime + j - 1, B + k - j);

        fmpz_mul_ui(u, Aden, k);
        fmpz_divexact(B + k, t, u);
    }

    _fmpq_poly_canonicalise(B, Bden, n);

    fmpz_clear(t);
    fmpz_clear(u);
}

/* Basecase algorithm; supports aliasing and guarantees canonical output. */
void
_fmpq_poly_exp_series_basecase(fmpz * B, fmpz_t Bden,
    const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    fmpz * Aprime;
    fmpz_t Aden2;

    Alen = FLINT_MIN(Alen, n);

    Aprime = _fmpz_vec_init(Alen - 1);
    fmpz_init(Aden2);

    /* There is probably not much content, so avoid the overhead of canonicalising. */
    if (Alen <= 6)
    {
        _fmpz_poly_derivative(Aprime, A, Alen);
        fmpz_set(Aden2, Aden);
    }
    else
    {
        _fmpq_poly_derivative(Aprime, Aden2, A, Aden, Alen);
    }

    _fmpq_poly_exp_series_basecase_deriv(B, Bden, Aprime, Aden2, Alen, n);

    _fmpz_vec_clear(Aprime, Alen - 1);
    fmpz_clear(Aden2);
}

/* c_k x^k -> c_k x^k / (m+k) */
void _fmpq_poly_integral_offset(fmpz * rpoly, fmpz_t rden,
                           const fmpz * poly, const fmpz_t den, slong len, slong m)
{
    slong k;
    ulong v, c, d;
    mp_ptr divisors;
    fmpz_t t, u;
    TMP_INIT;

    TMP_START;
    divisors = TMP_ALLOC(sizeof(ulong) * len);

    fmpz_init(t);
    fmpz_one(t);

    for (k = len - 1; k >= 0; k--)
    {
        if (fmpz_is_zero(poly + k))
        {
            fmpz_zero(rpoly + k);
        }
        else
        {
            c = _fmpz_gcd_small(poly + k, k + m);

            if (c == k + m)
            {
                fmpz_divexact_ui(rpoly + k, poly + k, k + m);
                divisors[k] = 1;
            }
            else
            {
                if (c == 1)
                {
                    fmpz_set(rpoly + k, poly + k);
                    divisors[k] = k + m;
                }
                else
                {
                    fmpz_divexact_ui(rpoly + k, poly + k, c);
                    divisors[k] = (k + m) / c;
                }

                c = divisors[k];
                d = _fmpz_gcd_small(t, c);
                if (d != c)
                    fmpz_mul_ui(t, t, c / d);
            }
        }
    }

    fmpz_mul(rden, den, t);

    if (!fmpz_is_one(t))
    {
        if (!COEFF_IS_MPZ(*t))
        {
            v = *t;
            for (k = len - 1; k >= 0; k--)
            {
                if (!fmpz_is_zero(rpoly + k) && v != divisors[k])
                    fmpz_mul_ui(rpoly + k, rpoly + k, divisors[k] == 1 ? v : v / divisors[k]);
            }
        }
        else
        {
            fmpz_init(u);

            for (k = len - 1; k >= 0; k--)
            {
                if (!fmpz_is_zero(rpoly + k))
                {
                    if (divisors[k] == 1)
                    {
                        fmpz_mul(rpoly + k, rpoly + k, t);
                    }
                    else
                    {
                        fmpz_divexact_ui(u, t, divisors[k]);
                        fmpz_mul(rpoly + k, rpoly + k, u);
                    }
                }
            }

            fmpz_clear(u);
        }
    }

    fmpz_clear(t);
    TMP_END;
}

static void
MULLOW(fmpz * z, fmpz_t zden, const fmpz * x, const fmpz_t xden, slong xn, const fmpz * y, const fmpz_t yden, slong yn, slong n)
{
    if (xn + yn - 1 < n)
        flint_throw(FLINT_ERROR, "(%s)\n", __func__);

    if (xn >= yn)
        _fmpz_poly_mullow(z, x, xn, y, yn, n);
    else
        _fmpz_poly_mullow(z, y, yn, x, xn, n);

    fmpz_mul(zden, xden, yden);
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

/* Newton iteration. If g == NULL, computes {f, fden, n} = exp({h, hden, hlen}).
   If g != NULL, simultaneously computes {g, gden, n} = exp(-{h, hden, hlen}).
   Allows aliasing between (f, fden) and (h, hden) but not with (g, gden). */
void
_fmpq_poly_exp_series_newton(fmpz * f, fmpz_t fden,
    fmpz * g, fmpz * gden,
    const fmpz * h, const fmpz_t hden,
    slong hlen, slong n)
{
    slong a[FLINT_BITS];
    slong original_n, i, m, l, r, cutoff;
    fmpz * t, * hprime;
    fmpz_t tden, hprimeden, uden, d, hdenin;
    int inverse;

    /* If g is provided, we compute g = exp(-h), and we can use g as
       scratch space. Otherwise, we still need to compute exp(-h) to length
       (n+1)/2 for intermediate use, and we still need n coefficients of
       scratch space. */
    original_n = n;
    inverse = (g != NULL);
    if (!inverse)
    {
        g = _fmpz_vec_init(n + 1);
        gden = g + n;
    }

    hlen = FLINT_MIN(hlen, n);

    t = _fmpz_vec_init(n);
    hprime = _fmpz_vec_init(hlen - 1);
    fmpz_init(tden);
    fmpz_init(hprimeden);
    fmpz_init(uden);
    fmpz_init(d);
    fmpz_init(hdenin);

    /* Precompute h' which is needed throughout. We will not
       canonicalise immediately; we want the truncated series
       to have minimal content in each step of the Newton iteration,
       so we canonicalise gradually. */
    fmpz_set(hdenin, hden);
    fmpz_set(hprimeden, hden);
    _fmpz_poly_derivative(hprime, h, hlen);

    cutoff = 20 + 1000 / n_sqrt(fmpz_bits(hden));

    for (i = 1; (WORD(1) << i) < n; i++);
    a[i = 0] = n;
    while (n >= cutoff || i == 0)
        a[++i] = (n = (n + 1) / 2);

    /* Canonicalise h' for first step. */
    _fmpq_poly_canonicalise(hprime, hprimeden, FLINT_MIN(n, hlen) - 1);

    /* Initial approximation f := exp(h) + O(x^n) using basecase algorithm. */
    _fmpq_poly_exp_series_basecase_deriv(f, fden, hprime, hprimeden, hlen, n);

    /* Initial approximation of inverse g := exp(-h) + O(x^n) */
    _fmpq_poly_inv_series(g, gden, f, fden, n, n);

    for (i--; i >= 0; i--)
    {
        m = n;             /* previous length */
        n = a[i];          /* new length */

        l = FLINT_MIN(hlen, n) - 1;
        r = FLINT_MIN(l + m - 1, n - 1);

        /* Extend h' */
        if (l > m - 1)
        {
            fmpz_set(uden, hdenin);
            _fmpq_poly_canonicalise(hprime + m - 1, uden, l - m + 1);
            CONCATENATE(hprime, hprimeden, uden, m - 1, l);
        }

        MULLOW(t, tden, hprime, hprimeden, l, f, fden, m, r);
        _fmpq_poly_canonicalise(t + m - 1, tden, r + 1 - m);
        MULLOW(g + m, uden, g, gden, n - m, t + m - 1, tden, r + 1 - m, n - m);
        _fmpq_poly_canonicalise(g + m, uden, n - m);
        _fmpq_poly_integral_offset(g + m, uden, g + m, uden, n - m, m);
        MULLOW(f + m, uden, f, fden, n - m, g + m, uden, n - m, n - m);
        /* Assuming that the low part is canonicalised on input,
           we just need to canonicalise the high part. */
        _fmpq_poly_canonicalise(f + m, uden, n - m);
        CONCATENATE(f, fden, uden, m, n);

        /* g := exp(-h) + O(x^n); not needed if we only want exp(x) */
        if (i != 0 || inverse)
        {
            MULLOW(t, tden, f, fden, n, g, gden, m, n);
            _fmpq_poly_canonicalise(t + m, tden, n - m);
            MULLOW(g + m, uden, g, gden, m, t + m, tden, n - m, n - m);
            /* Assuming that the low part is canonicalised on input,
               we just need to canonicalise the high part. */
            _fmpq_poly_canonicalise(g + m, uden, n - m);
            CONCATENATE(g, gden, uden, m, n);
            _fmpz_vec_neg(g + m, g + m, n - m);
        }
    }

    _fmpz_vec_clear(hprime, hlen - 1);
    _fmpz_vec_clear(t, original_n);
    fmpz_clear(tden);
    fmpz_clear(hprimeden);
    fmpz_clear(uden);
    fmpz_clear(d);
    fmpz_clear(hdenin);
    if (!inverse)
        _fmpz_vec_clear(g, original_n + 1);
}

void
_fmpq_poly_exp_series(fmpz * B, fmpz_t Bden,
    const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    Alen = FLINT_MIN(Alen, n);

    if (Alen == 1)
    {
        fmpz_one(B);
        fmpz_one(Bden);
        _fmpz_vec_zero(B + 1, n - 1);
        return;
    }

    /* A is a monomial (p/q) * x^d */
    if (_fmpz_vec_is_zero(A + 1, Alen - 2))
    {
        fmpz * R;
        ulong v;
        slong i, d, m;

        d = Alen - 1;       /* Degree of input monomial. */
        m = (n - 1) / d;    /* m*d is highest degree in output. */
        R = _fmpz_vec_init(m + 1);

        /* A[d]/Aden could be non-canonical due to truncation */
        fmpz_gcd(R, A + d, Aden);
        fmpz_divexact(B + d, A + d, R);
        fmpz_divexact(R, Aden, R);  /* store q in R[0] */

        fmpz_set(R + 1, R);
        fmpz_set(Bden, R);

        for (i = 2; i <= m; i++)
        {
            /* Computing (p/q)^i / i! from the previous term, we only
               need to remove the gcd between the numerator and i. */
            fmpz_mul(B + i * d, B + (i - 1) * d, B + d);
            fmpz_mul(Bden, Bden, R);

            v = _fmpz_gcd_small(B + i * d, i);
            fmpz_divexact_ui(B + i * d, B + i * d, v);
            fmpz_mul_ui(Bden, Bden, i / v);
            fmpz_mul_ui(R + i, R, i / v);
        }

        /* Put all terms on the same denominator as the last term. */
        for (i = m - 1; i > 0; i--)
        {
            fmpz_mul(B + i * d, B + i * d, R + m);
            fmpz_mul(R + m, R + m, R + i);
        }

        /* Constant term = 1. */
        fmpz_set(B, Bden);

        for (i = 0; d != 1 && i < n; i++)
            if (i % d != 0)
                fmpz_zero(B + i);

        _fmpz_vec_clear(R, m + 1);
        return;
    }

    if (Alen <= 12 || n <= 10 + 1000 / n_sqrt(fmpz_bits(Aden)))
    {
        _fmpq_poly_exp_series_basecase(B, Bden, A, Aden, Alen, n);
    }
    else
    {
        _fmpq_poly_exp_series_newton(B, Bden, NULL, NULL, A, Aden, Alen, n);
    }
}

void
_fmpq_poly_exp_expinv_series(fmpz * B, fmpz_t Bden, fmpz * C, fmpz_t Cden,
    const fmpz * A, const fmpz_t Aden, slong Alen, slong n)
{
    Alen = FLINT_MIN(Alen, n);

    if (Alen == 1)
    {
        fmpz_one(B);
        fmpz_one(Bden);
        fmpz_one(C);
        fmpz_one(Cden);
        _fmpz_vec_zero(B + 1, n - 1);
        _fmpz_vec_zero(C + 1, n - 1);
        return;
    }

    if (_fmpz_vec_is_zero(A + 1, Alen - 2))
    {
        slong i;
        _fmpq_poly_exp_series(B, Bden, A, Aden, Alen, n);
        _fmpz_vec_set(C, B, n);
        for (i = Alen - 1; i < n; i += 2 * (Alen - 1))
            fmpz_neg(C + i, C + i);
        fmpz_set(Cden, Bden);
        return;
    }

    /* todo: tweak tuning for this function */
    if (Alen <= 12 || n <= 10 + 1000 / n_sqrt(fmpz_bits(Aden)))
    {
        _fmpq_poly_exp_series_basecase(B, Bden, A, Aden, Alen, n);
        _fmpq_poly_inv_series(C, Cden, B, Bden, n, n);
    }
    else
    {
        fmpz * tmp;
        if (A == C || Aden == Cden)
        {
            tmp = _fmpz_vec_init(n + 1);
            _fmpq_poly_exp_series_newton(B, Bden, tmp, tmp + n, A, Aden, Alen, n);
            _fmpz_vec_swap(C, tmp, n);
            fmpz_swap(Cden, tmp + n);
            _fmpz_vec_clear(tmp, n);
        }
        else
        {
            _fmpq_poly_exp_series_newton(B, Bden, C, Cden, A, Aden, Alen, n);
        }
    }
}

void fmpq_poly_exp_series(fmpq_poly_t res, const fmpq_poly_t poly, slong n)
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
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_exp_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res, n);
    _fmpq_poly_exp_series(res->coeffs, res->den,
        poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res, n);
    _fmpq_poly_normalise(res);
}

void fmpq_poly_exp_expinv_series(fmpq_poly_t res1, fmpq_poly_t res2, const fmpq_poly_t poly, slong n)
{
    if (n == 0)
    {
        fmpq_poly_zero(res1);
        fmpq_poly_zero(res2);
        return;
    }

    if (poly->length == 0 || n == 1)
    {
        fmpq_poly_one(res1);
        fmpq_poly_one(res2);
        return;
    }

    if (!fmpz_is_zero(poly->coeffs))
    {
        flint_throw(FLINT_ERROR, "Exception (fmpq_poly_exp_expinv_series). Constant term != 0.\n");
    }

    fmpq_poly_fit_length(res1, n);
    fmpq_poly_fit_length(res2, n);
    _fmpq_poly_exp_expinv_series(res1->coeffs, res1->den,
                          res2->coeffs, res2->den,
        poly->coeffs, poly->den, poly->length, n);
    _fmpq_poly_set_length(res1, n);
    _fmpq_poly_set_length(res2, n);
    _fmpq_poly_normalise(res1);
    _fmpq_poly_normalise(res2);
}
