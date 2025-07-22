/*
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "longlong.h"
#include "flint.h"
#include "fmpz.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "fmpq_poly.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "ulong_extras.h"
#include "gr.h"
#include "gr_poly.h"

flint_bitcnt_t
_fmpq_vec_max_height_bits(const fmpq * x, slong len)
{
    flint_bitcnt_t bits, max_bits = 0;

    for (slong i = 0; i < len; i++)
    {
        bits = fmpq_height_bits(x + i);
        if (bits > max_bits)
            max_bits = bits;
    }

    return max_bits;
}

int _is_lucky_prime_ui(ulong p, const fmpq * xs, slong n)
{
    int ok = 1;
    slong i;
    for (i = 0; i < n && ok; i++)
        ok = !fmpz_divisible_ui(fmpq_denref(xs + i), p);
    return ok;
}


static void
_fmpz_ui_vec_prod(fmpz_t res, nn_srcptr x, slong n)
{
    if (n < 16)
    {
        slong i;
        fmpz_set_ui(res, x[0]);
        for (i = 1; i < n; i++)
            fmpz_mul_ui(res, res, x[i]);
    }
    else
    {
        fmpz_t t;
        fmpz_init(t);
        _fmpz_ui_vec_prod(res, x, n / 2);
        _fmpz_ui_vec_prod(t, x + n / 2, n - n / 2);
        fmpz_mul(res, res, t);
        fmpz_clear(t);
    }
}

static int
_fmpq_poly_check_interpolant(const fmpz * poly, const fmpz_t den,
                                const fmpq * xs, const fmpq * ys, slong n)
{
    fmpz_t ynum, yden;
    fmpq_t y1;
    slong i;
    int ok = 1;

    fmpq_init(y1);
    fmpz_init(ynum);
    fmpz_init(yden);

    for (i = 0; i < n && ok; i++)
    {
        _fmpz_poly_evaluate_fmpq(ynum, yden, poly, n,
                                            fmpq_numref(xs + i), fmpq_denref(xs + i));
        fmpq_mul_fmpz(y1, ys + i, den);
        fmpq_mul_fmpz(y1, y1, yden);
        ok = fmpq_equal_fmpz(y1, ynum);
    }
    fmpq_clear(y1);
    fmpz_clear(ynum);
    fmpz_clear(yden);

    return ok;
}

static int
_fmpq_vec_has_unique_entries(const fmpq * x, slong n)
{
    fmpq * t;
    slong i;
    int ok = 1;

    t = _fmpq_vec_init(n);
    for (i = 0; i < n; i++)
        fmpq_set(t + i, x + i);
    _fmpq_vec_sort(t, n);
    for (i = 1; i < n; i++)
        if (fmpq_equal(t + i - 1, t + i))
            ok = 0;

    _fmpq_vec_clear(t, n);

    return ok;
}

/* _nmod_poly_interpolate_nmod_vec currently lacks error handling */
int
_checked_nmod_poly_interpolate(nn_ptr r, nn_srcptr x, nn_srcptr y, slong n, nmod_t mod);

void
_fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2,
                   const fmpz_t m2, const fmpz_t m1m2, fmpz_t c, int sign);

int
_fmpq_poly_interpolate_multi_mod(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    ulong p, mden;
    slong i, j, k;
    nmod_t mod;
    nn_ptr xm, ym;
    fmpz_t M, t, u, c, M2, M1M2;
    slong total_primes, num_primes, count_good;
    flint_bitcnt_t xbits, ybits, bound_bits;
    int ok = 1, rat_rec;
    int checked_unique = 0;
    nn_ptr primes = NULL, xmod = NULL, ymod = NULL, residues = NULL, residuesden = NULL;
    int * good = NULL;
    fmpq * coeffs;

    if (n == 0)
        return 1;

    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(c);
    fmpz_init(M2);
    fmpz_init(M1M2);

    coeffs = _fmpq_vec_init(n);

    /* Lagrange interpolation gives:

        w_i = prod_{j != i} (x - x_i) / (x_j - x_i)
        hb = heightbits

        hb(w_i) <= max_i hb(x_i) * 2(n-1)

        hb(f) <= sum_i hb(y_i) + hb(w_i) * 2(n -1)
              <= n * (hb(y_i) + hb(x_i) * 4(n -1))

    */

    xbits = _fmpq_vec_max_height_bits(xs, n);
    ybits = _fmpq_vec_max_height_bits(ys, n);
    bound_bits = n * (xbits * 4 * (n - 1) + ybits);
    /* Rational reconstruction needs modulus M >= 1/2*height(f)^2 + 1 */
    bound_bits *= 2;

    xm = _nmod_vec_init(n);
    ym = _nmod_vec_init(n);

    total_primes = 0;
    /* Important: when changing this, change the list of adversarial primes
       in the test code to match. */
    p = UWORD(1) << (FLINT_BITS - 1);
    for (;;)
    {
        if (total_primes < 16)
        {
            do /* p must not divide any denominator */
            {
                p = n_nextprime(p, 1);
            } while (!(_is_lucky_prime_ui(p, xs, n) && _is_lucky_prime_ui(p, ys, n)));

            nmod_init(&mod, p);
            for (i = 0; i < n; i++) {
                xm[i] = fmpz_get_nmod(fmpq_numref(xs + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(xs + i), mod);
                xm[i] = nmod_div(xm[i], mden, mod);

                ym[i] = fmpz_get_nmod(fmpq_numref(ys + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(ys + i), mod);
                ym[i] = nmod_div(ym[i], mden, mod);
            }
            if (_checked_nmod_poly_interpolate(ym, xm, ym, n, mod))
            {
                num_primes = 1;
                if (total_primes == 0)
                {
                    _fmpz_vec_set_nmod_vec(poly, ym, n, mod);
                    fmpz_set_ui(M, p);
                }
                else
                {
                    _fmpz_poly_CRT_ui(poly, poly, n, M, ym, n, mod.n, mod.ninv, 0);
                    fmpz_mul_ui(M, M, p);
                }
            }
            else
                num_primes = 0;
        }
        else  /* Do a batch of primes at once with fast reduction/CRT */
        {
            fmpz_comb_t comb;
            fmpz_comb_temp_t temp;

            /* Todo: when we have a good bound, don't wildly overshoot. */
            num_primes = FLINT_MAX(1, total_primes / 2);

            primes = flint_realloc(primes, sizeof(ulong) * num_primes);
            xmod = flint_realloc(xmod, sizeof(ulong) * n * num_primes);
            ymod = flint_realloc(ymod, sizeof(ulong) * n * num_primes);
            residues = flint_realloc(residues, sizeof(ulong) * num_primes);
            residuesden = flint_realloc(residuesden, sizeof(ulong) * num_primes);
            good = flint_realloc(good, sizeof(int) * num_primes);

            for (k = 0; k < num_primes;)
            {
                p = n_nextprime(p, 1);
                if (_is_lucky_prime_ui(p, xs, n) && _is_lucky_prime_ui(p, ys, n))
                {
                    primes[k] = p;
                    k++;
                }
            }

            _fmpz_ui_vec_prod(M2, primes, num_primes);

            fmpz_comb_init(comb, primes, num_primes);
            fmpz_comb_temp_init(temp, comb);

            // Reduce each j-th coeff mod all primes
            for (j = 0; j < n; j++)
                {
                fmpz_multi_mod_ui(residues, fmpq_numref(xs + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(xs + j), comb, temp);
                for (k = 0; k < num_primes; k++) {
                    nmod_init(&mod, primes[k]);
                    xmod[k * n + j] = nmod_div(residues[k], residuesden[k], mod);
                }

                fmpz_multi_mod_ui(residues, fmpq_numref(ys + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(ys + j), comb, temp);
                for (k = 0; k < num_primes; k++) {
                    nmod_init(&mod, primes[k]);
                    ymod[k * n + j] = nmod_div(residues[k], residuesden[k], mod);
                }
            }

            count_good = 0;
            for (k = 0; k < num_primes; k++)
            {
                nmod_init(&mod, primes[k]);
                good[k] = _checked_nmod_poly_interpolate(ymod + k * n, xmod + k * n, ymod + k * n, n, mod);
                count_good += (good[k] != 0);
            }

            /* Some primes failed; keep the ones that worked */
            if (count_good < num_primes)
            {
                count_good = 0;
                for (k = 0; k < num_primes; k++)
                {
                    if (good[k])
                    {
                        primes[count_good] = primes[k];
                        if (count_good != k)
                            _nmod_vec_set(ymod + count_good * n, ymod + k * n, n);
                        count_good++;
                    }
                }

                /* Reinitialize CRT for the product of correct primes */
                num_primes = count_good;

                if (num_primes != 0)
                {
                    _fmpz_ui_vec_prod(M2, primes, num_primes);
                    fmpz_comb_temp_clear(temp);
                    fmpz_comb_clear(comb);
                    fmpz_comb_init(comb, primes, num_primes);
                    fmpz_comb_temp_init(temp, comb);
                }
            }

            if (num_primes != 0)
            {
                fmpz_mul(M1M2, M, M2);
                fmpz_mod(c, M, M2);
                fmpz_invmod(c, c, M2);

                for (j = 0; j < n; j++)
                {
                    for (k = 0; k < num_primes; k++)
                        residues[k] = ymod[k * n + j];
                    // reconstruct modulo M2 = prod(primes)
                    fmpz_multi_CRT_ui(t, residues, comb, temp, 0);
                    // reconstruct modulo M1M2 = M*M2
                    _fmpz_CRT(u, poly + j, M, t, M2, M1M2, c, 0);
                    fmpz_swap(poly + j, u);
                }
                fmpz_swap(M, M1M2);
            }

            fmpz_comb_temp_clear(temp);
            fmpz_comb_clear(comb);
        }

        total_primes += num_primes;

        /* At that step we have in poly, the polynomial modulo M,
           which a product of total_primes primes */

        /* No primes succeeded -- verify that evaluation points are actually OK */
        if (num_primes == 0 && !checked_unique)
        {
            if (!_fmpq_vec_has_unique_entries(xs, n))
            {
                ok = 0;
                break;
            }
            checked_unique = 1;
        }
        else
        {
            /* Rational reconstruction */
            for (i = 0, rat_rec = 1; rat_rec && i < n; i++)
                rat_rec = fmpq_reconstruct_fmpz(coeffs + i, poly + i, M);

            if (rat_rec)
            {/* Check for termination */
                _fmpq_vec_get_fmpz_vec_fmpz(poly, den, coeffs, n);
                if (_fmpq_poly_check_interpolant(poly, den, xs, ys, n))
                    break;
            }
        }

        /* No solution */
        if (fmpz_bits(M) > bound_bits + 1)
        {
            ok = 0;
            break;
        }
    }

    fmpz_clear(M);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(c);
    fmpz_clear(M2);
    fmpz_clear(M1M2);

    _nmod_vec_clear(xm);
    _nmod_vec_clear(ym);

    _fmpq_vec_clear(coeffs, n);

    flint_free(primes);
    flint_free(xmod);
    flint_free(ymod);
    flint_free(residues);
    flint_free(residuesden);
    flint_free(good);

    return ok;
}

int
fmpq_poly_interpolate_multi_mod(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    int ok;
    fmpq_poly_fit_length(poly, n);
    ok = _fmpq_poly_interpolate_multi_mod(poly->coeffs, poly->den, xs, ys, n);
    _fmpq_poly_set_length(poly, n);
    _fmpq_poly_normalise(poly);
    return ok;
}
