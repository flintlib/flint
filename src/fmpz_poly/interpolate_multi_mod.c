/*
    Copyright (C) 2011, 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz.h"
#include "fmpz/impl.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"
#include "nmod_vec.h"
#include "gr_poly.h"

static int
_fmpz_poly_check_interpolant(const fmpz * poly, const fmpz * xs, const fmpz * ys, slong n)
{
    fmpz_t y;
    slong i;
    int ok = 1;

    fmpz_init(y);

    for (i = 0; i < n && ok; i++)
    {
        _fmpz_poly_evaluate_fmpz(y, poly, n, xs + i);
        ok = fmpz_equal(y, ys + i);
    }

    fmpz_clear(y);

    return ok;
}

static int
_fmpz_vec_has_unique_entries(const fmpz * x, slong n)
{
    fmpz * t;
    slong i;
    int ok = 1;

    t = _fmpz_vec_init(n);
    _fmpz_vec_set(t, x, n);
    _fmpz_vec_sort(t, n);
    for (i = 1; i < n; i++)
        if (fmpz_equal(t + i - 1, t + i))
            ok = 0;

    _fmpz_vec_clear(t, n);

    return ok;
}

/* _nmod_poly_interpolate_nmod_vec currently lacks error handling */
int
_checked_nmod_poly_interpolate(nn_ptr r, nn_srcptr x, nn_srcptr y, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, mod.n);
    return (_gr_poly_interpolate_fast(r, x, y, n, ctx) == GR_SUCCESS);
}

int
_fmpz_poly_interpolate_multi_mod(fmpz * poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    ulong p;
    slong j, k;
    nmod_t mod;
    nn_ptr xm, ym;
    fmpz_t M, t, u, c, M2, M1M2;
    slong total_primes, num_primes, count_good;
    slong xbits, ybits;
    int ok = 1;
    int checked_unique = 0;
    slong bound_bits;
    nn_ptr primes = NULL, xmod = NULL, ymod = NULL, residues = NULL;
    int * good = NULL;

    if (n == 0)
        return 1;

    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(c);
    fmpz_init(M2);
    fmpz_init(M1M2);

    /* Lagrange interpolation gives:

        w_i = prod_{j != i} (x - x_i) / (x_j - x_i)

        height(w_i) <= height((x - max |x_i|))^(n-1))
                    <= (2 max |x_i|)^(n-1)

        height(f) <= sum_i |y_i| height(w_i)
    */

    xbits = _fmpz_vec_max_bits(xs, n);
    xbits = FLINT_ABS(xbits);
    ybits = _fmpz_vec_max_bits(ys, n);
    ybits = FLINT_ABS(ybits);
    bound_bits = (xbits + 1) * (n - 1) + FLINT_BIT_COUNT(n) * ybits;

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
            p = n_nextprime(p, 1);
            nmod_init(&mod, p);

            _fmpz_vec_get_nmod_vec(xm, xs, n, mod);
            _fmpz_vec_get_nmod_vec(ym, ys, n, mod);

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
                    _fmpz_poly_CRT_ui(poly, poly, n, M, ym, n, mod.n, mod.ninv, 1);
                    fmpz_mul_ui(M, M, p);
                }
            }
            else
            {
                num_primes = 0;
            }
        }
        else        /* Do a batch of primes at once with fast reduction/CRT */
        {
            fmpz_comb_t comb;
            fmpz_comb_temp_t temp;

            /* Todo: when we have a good bound, don't wildly overshoot. */
            num_primes = FLINT_MAX(1, total_primes / 2);

            primes = flint_realloc(primes, sizeof(ulong) * num_primes);
            xmod = flint_realloc(xmod, sizeof(ulong) * n * num_primes);
            ymod = flint_realloc(ymod, sizeof(ulong) * n * num_primes);
            residues = flint_realloc(residues, sizeof(ulong) * num_primes);
            good = flint_realloc(good, sizeof(int) * num_primes);

            for (k = 0; k < num_primes; k++)
            {
                p = n_nextprime(p, 1);
                primes[k] = p;
            }

            _fmpz_ui_vec_prod(M2, primes, num_primes);

            fmpz_comb_init(comb, primes, num_primes);
            fmpz_comb_temp_init(temp, comb);

            for (j = 0; j < n; j++)
            {
                /* todo: do iterative evaluation may perform better if small values */
                fmpz_multi_mod_ui(residues, xs + j, comb, temp);
                for (k = 0; k < num_primes; k++)
                    xmod[k * n + j] = residues[k];

                fmpz_multi_mod_ui(residues, ys + j, comb, temp);
                for (k = 0; k < num_primes; k++)
                    ymod[k * n + j] = residues[k];
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

                    fmpz_multi_CRT_ui(t, residues, comb, temp, 0);
                    _fmpz_CRT(u, poly + j, M, t, M2, M1M2, c, 1);
                    fmpz_set(poly + j, u);
                }

                fmpz_swap(M, M1M2);
            }

            fmpz_comb_temp_clear(temp);
            fmpz_comb_clear(comb);
        }

        total_primes += num_primes;

        /* No primes succeeded -- verify that evaluation points are actually OK */
        if (num_primes == 0 && !checked_unique)
        {
            if (!_fmpz_vec_has_unique_entries(xs, n))
            {
                ok = 0;
                break;
            }
            checked_unique = 1;
        }
        else
        {
            /* Check for termination */
            if (_fmpz_poly_check_interpolant(poly, xs, ys, n))
                break;
    
        }

        /* No solution */
        if ((slong) fmpz_bits(M) > bound_bits + 1)
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

    flint_free(primes);
    flint_free(xmod);
    flint_free(ymod);
    flint_free(residues);
    flint_free(good);

    return ok;
}

int
fmpz_poly_interpolate_multi_mod(fmpz_poly_t poly,
                                    const fmpz * xs, const fmpz * ys, slong n)
{
    int ok;
    fmpz_poly_fit_length(poly, n);
    ok = _fmpz_poly_interpolate_multi_mod(poly->coeffs, xs, ys, n);
    _fmpz_poly_set_length(poly, n);
    _fmpz_poly_normalise(poly);
    return ok;
}

