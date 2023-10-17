/*
    Copyright (C) 2011, 2012 Sebastian Pancratz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "qadic.h"

// Function to compute the power of a q-adic number
void _qadic_pow(fmpz *result, const fmpz *base, slong length, const fmpz_t exponent,
               const fmpz *a, const slong *j, slong a_length, const fmpz_t p) {
    const slong degree = j[a_length - 1];

    if (fmpz_is_zero(exponent)) {
        fmpz_one(result);
        _fmpz_vec_zero(result + 1, 2 * degree - 1 - 1);
    } else if (fmpz_is_one(exponent)) {
        _fmpz_vec_set(result, base, length);
        _fmpz_vec_zero(result + length, 2 * degree - 1 - length);
    } else {
        ulong msb = fmpz_bits(exponent) - 1;
        ulong bit;
        fmpz *temp = _fmpz_vec_init(2 * degree - 1);
        fmpz *R, *S, *T;

        _fmpz_vec_zero(result, 2 * degree - 1);

        // Set bits to the bitmask with a 1 one place lower than the msb of exponent
        bit = msb - 1;

        // Trial run without any polynomial arithmetic to determine the parity
        // of the number of swaps; then set R and S accordingly
        {
            unsigned int swaps = 0U;
            ulong bit2 = bit;
            if (fmpz_tstbit(exponent, bit2))
                swaps = ~swaps;
            while (bit2--)
                if (!fmpz_tstbit(exponent, bit2))
                    swaps = ~swaps;

            if (swaps == 0U) {
                R = result;
                S = temp;
            } else {
                R = temp;
                S = result;
            }
        }

        // Unroll the first step of the loop
        _fmpz_poly_sqr(R, base, length);
        _fmpz_mod_poly_reduce(R, 2 * length - 1, a, j, a_length, p);

        if (fmpz_tstbit(exponent, bit)) {
            _fmpz_poly_mul(S, R, degree, base, length);
            _fmpz_mod_poly_reduce(S, degree + length - 1, a, j, a_length, p);
            T = R;
            R = S;
            S = T;
        }

        while (bit--) {
            if (fmpz_tstbit(exponent, bit)) {
                _fmpz_poly_sqr(S, R, degree);
                _fmpz_mod_poly_reduce(S, 2 * degree - 1, a, j, a_length, p);
                _fmpz_poly_mul(R, S, degree, base, length);
                _fmpz_mod_poly_reduce(R, degree + length - 1, a, j, a_length, p);
            } else {
                _fmpz_poly_sqr(S, R, degree);
                _fmpz_mod_poly_reduce(S, 2 * degree - 1, a, j, a_length, p);
                T = R;
                R = S;
                S = T;
            }
        }

        _fmpz_vec_clear(temp, 2 * degree - 1);
    }
}

// Function to compute the power of a q-adic number
void qadic_pow(qadic_t result, const qadic_t base, const fmpz_t exponent, const qadic_ctx_t ctx) {
    const slong N = qadic_prec(result);

    if (fmpz_sgn(exponent) < 0) {
        flint_printf("Exception (qadic_pow). e < 0.\n");
        flint_abort();
    }

    if (fmpz_is_zero(exponent)) {
        qadic_one(result);
    } else if (qadic_is_zero(base)) {
        qadic_zero(result);
    } else {
        fmpz_t val;  /* N - exponent * val(base) */

        fmpz_init_set(val, exponent);
        fmpz_mul_si(val, val, base->val);

        if (fmpz_cmp_si(val, N) >= 0) {
            qadic_zero(result);
        } else if (fmpz_is_one(exponent)) {
            qadic_set(result, base, ctx);
        } else {
            const slong degree = qadic_ctx_degree(ctx);
            fmpz *t;
            fmpz_t pow;
            int alloc;

            alloc = _padic_ctx_pow_ui(pow, N - fmpz_get_si(val), &ctx->pctx);

            if (result == base) {
                t = _fmpz_vec_init(2 * degree - 1);
            } else {
                padic_poly_fit_length(result, 2 * degree - 1);
                t = result->coeffs;
            }

            _qadic_pow(t, base->coeffs, base->length, exponent, ctx->a, ctx->j, ctx->len, pow);
            result->val = fmpz_get_si(val);

            if (result == base) {
                _fmpz_vec_clear(result->coeffs, result->alloc);
                result->coeffs = t;
                result->alloc  = 2 * degree - 1;
                result->length = degree;
            } else {
                _padic_poly_set_length(result, degree);
            }
            _padic_poly_normalise(result);

            if (alloc)
                fmpz_clear(pow);
        }
        fmpz_clear(val);
    }
}

