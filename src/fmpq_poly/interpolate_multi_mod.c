/*
    Copyright (C) 2011, 2025 Fredrik Johansson
    Copyright (C) 2025 Rémi Prébet

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz/impl.h"
#include "fmpq.h"
#include "fmpz_vec.h"
#include "fmpq_vec.h"
#include "fmpz_poly.h"
#include "fmpz_poly/impl.h"
#include "fmpq_poly.h"
#include "nmod_vec.h"

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

void
_fmpq_poly_interpolate_multi_mod(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    ulong p, mden;
    slong i, j, k;
    nmod_t mod;
    nmod_t * Lmod = NULL;
    nn_ptr xm, ym;
    fmpz_t M, t, u, c, M2, M1M2;
    slong total_primes, num_primes, count_good;
    int rat_rec;
    nn_ptr primes = NULL, xmod = NULL, ymod = NULL, residues = NULL, residuesden = NULL;
    int * good = NULL;
    fmpq * coeffs;
    fmpz * rescoeffs;

    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(c);
    fmpz_init(M2);
    fmpz_init(M1M2);

    coeffs = _fmpq_vec_init(n);
    rescoeffs = _fmpz_vec_init(n);

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
            mden = 0;
            nmod_init(&mod, p);
            /* Modular reduction */
            for (i = 0; i < n; i++) {
                xm[i] = fmpz_get_nmod(fmpq_numref(xs + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(xs + i), mod);
                if (!mden) break;
                xm[i] = nmod_div(xm[i], mden, mod);
                ym[i] = fmpz_get_nmod(fmpq_numref(ys + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(ys + i), mod);
                if (!mden) break;
                ym[i] = nmod_div(ym[i], mden, mod);
            }
            /* Modular interpolation and CRT */
            if (mden && _checked_nmod_poly_interpolate(ym, xm, ym, n, mod))
            {
                num_primes = 1;
                if (total_primes == 0)
                {
                    for (i = 0; i < n; i++)
                        fmpz_set_ui(rescoeffs + i, ym[i]);
                    fmpz_set_ui(M, p);
                }
                else
                {
                    _fmpz_poly_CRT_ui(rescoeffs, rescoeffs, n, M, ym, n, mod.n, mod.ninv, 0);
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
            Lmod = flint_realloc(Lmod, sizeof(nmod_t) * num_primes);
            xmod = flint_realloc(xmod, sizeof(ulong) * n * num_primes);
            ymod = flint_realloc(ymod, sizeof(ulong) * n * num_primes);
            residues = flint_realloc(residues, sizeof(ulong) * num_primes);
            residuesden = flint_realloc(residuesden, sizeof(ulong) * num_primes);
            good = flint_realloc(good, sizeof(int) * num_primes);

            for (k = 0; k < num_primes; k++)
            {
                p = n_nextprime(p, 1);
                nmod_init(Lmod + k, p);
                primes[k] = p;
                good[k] = 1;
            }

            _fmpz_ui_vec_prod(M2, primes, num_primes);

            fmpz_comb_init(comb, primes, num_primes);
            fmpz_comb_temp_init(temp, comb);

            // Reduce each j-th coeff mod all primes (report zero denominators)
            for (j = 0; j < n; j++)
            {
                fmpz_multi_mod_ui(residues, fmpq_numref(xs + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(xs + j), comb, temp);
                for (k = 0; k < num_primes; k++)
                {
                    if (residuesden[k] && good[k])
                        xmod[k * n + j] = nmod_div(residues[k], residuesden[k], Lmod[k]);
                    else
                        good[k] = 0;
                }

                fmpz_multi_mod_ui(residues, fmpq_numref(ys + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(ys + j), comb, temp);
                for (k = 0; k < num_primes; k++)
                {
                    if (residuesden[k] && good[k])
                        ymod[k * n + j] = nmod_div(residues[k], residuesden[k], Lmod[k]);
                    else
                        good[k] = 0;
                }
            }

            count_good = 0;
            for (k = 0; k < num_primes; k++)
            {
                if (!good[k]) continue;
                good[k] = _checked_nmod_poly_interpolate(ymod + k * n, xmod + k * n, ymod + k * n, n, Lmod[k]);
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
                    _fmpz_CRT(u, rescoeffs + j, M, t, M2, M1M2, c, 0);
                    fmpz_swap(rescoeffs + j, u);
                }
                fmpz_swap(M, M1M2);
            }

            fmpz_comb_temp_clear(temp);
            fmpz_comb_clear(comb);
        }

        total_primes += num_primes;

        /* At that step we have in poly, the polynomial modulo M,
           which a product of total_primes primes */

        if (num_primes > 0)
        {
            /* Rational reconstruction */
            for (i = 0, rat_rec = 1; rat_rec && i < n; i++)
                rat_rec = fmpq_reconstruct_fmpz(coeffs + i, rescoeffs + i, M);

            if (rat_rec)
            {/* Check for termination */
                _fmpq_vec_get_fmpz_vec_fmpz(poly, den, coeffs, n);
                if (_fmpq_poly_check_interpolant(poly, den, xs, ys, n))
                    break;
            }
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

    _fmpz_vec_clear(rescoeffs,n);
    _fmpq_vec_clear(coeffs, n);

    flint_free(primes);
    flint_free(xmod);
    flint_free(ymod);
    flint_free(residues);
    flint_free(residuesden);
    flint_free(good);
    flint_free(Lmod);
}

void
fmpq_poly_interpolate_multi_mod(fmpq_poly_t poly,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    if (n == 0)
        fmpq_poly_zero(poly);
    else
    {
        fmpq_poly_fit_length(poly, n);
        _fmpq_poly_interpolate_multi_mod(poly->coeffs, poly->den, xs, ys, n);
        _fmpq_poly_set_length(poly, n);
        _fmpq_poly_normalise(poly);
    }
}
