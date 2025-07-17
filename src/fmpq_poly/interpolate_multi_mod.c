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
#include "profiler.h"

ulong _fmpq_vec_max_height(const fmpq * vec, slong n) {
    ulong h = 0;

    for (slong i = 0; i < n; i++)
        h = FLINT_MAX(h, fmpq_height_bits(vec + i));

    return h;
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
_fmpq_poly_check_interpolant(const fmpq * coeffs,
                                const fmpq * xs, const fmpq * ys, slong n)
{
    fmpq_t y;
    fmpz * poly;
    fmpz_t den;
    slong i;
    int ok = 1;

    fmpq_init(y);
    fmpz_init(den);
    poly = _fmpz_vec_init(n);

    _fmpq_vec_get_fmpz_vec_fmpz(poly, den, coeffs, n);

    for (i = 0; i < n && ok; i++)
    {
        _fmpq_poly_evaluate_fmpq(fmpq_numref(y), fmpq_denref(y), poly, den, n,
                                 fmpq_numref(xs + i), fmpq_denref(xs + i));
        ok = fmpq_equal(y, ys + i);
    }

    fmpq_clear(y);

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
_checked_nmod_poly_interpolate_bis(nn_ptr r, nn_srcptr x, nn_srcptr y, slong n, nmod_t mod)
{
    gr_ctx_t ctx;
    gr_ctx_init_nmod(ctx, mod.n);
    return (_gr_poly_interpolate_fast(r, x, y, n, ctx) == GR_SUCCESS);
}


void
_fmpz_CRT(fmpz_t out, const fmpz_t r1, const fmpz_t m1, const fmpz_t r2,
                   const fmpz_t m2, const fmpz_t m1m2, fmpz_t c, int sign);

slong
_fmpq_vec_max_height_bits(const fmpq * x, slong len)
{
    slong i, bits, max_bits = 0;

    for (i = 0; i < len; i++)
    {
        bits = fmpq_height_bits(x + i);
        if (bits > max_bits)
            max_bits = bits;
    }

    return max_bits;
}

int
_fmpq_poly_interpolate_multi_mod(fmpz * poly, fmpz_t den,
                                    const fmpq * xs, const fmpq * ys, slong n)
{
    ulong p, mden;
    slong i, j, k;
    nmod_t mod;
    nn_ptr xm, ym;
    fmpz_t M, t, u, c, M2, M1M2;
    nn_ptr tmp;
    slong total_primes, num_primes, total_skipped, count_good;
    slong xbits, ybits;
    int ok = 1, rat_rec;
    int checked_unique = 0;
    ulong bound_bits;
    nn_ptr primes = NULL, xmod = NULL, ymod = NULL, residues = NULL, residuesden = NULL;
    int * good = NULL;
    fmpq * coeffs;
    fmpz_t old_res, curr_res;
    fmpq_t old_coeffs, curr_coeff; // We track the last coeff for now
    nn_ptr Lprimes, Lres;

    fmpz *as, *bs, *cs, *ds;

    if (n == 0)
        return 1;

    fmpz_init(M);
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(c);
    fmpz_init(M2);
    fmpz_init(M1M2);

    fmpz_init(old_res);
    fmpz_init(curr_res);
    fmpq_init(old_coeffs);
    fmpq_init(curr_coeff);

    as = _fmpz_vec_init(n);
    bs = _fmpz_vec_init(n);
    cs = _fmpz_vec_init(n);
    ds = _fmpz_vec_init(n);

    timeit_t t0, t1;

    coeffs = _fmpq_vec_init(n);
    for (i = 0; i < n; i++) {
        fmpz_set(as + i, fmpq_numref(xs + i));
        fmpz_set(bs + i, fmpq_denref(xs + i));
        fmpz_set(cs + i, fmpq_numref(ys + i));
        fmpz_set(ds + i, fmpq_denref(ys + i));
    }
    //_fmpq_vec_print(xs, n); flint_printf("\n");
    //_fmpq_vec_print(ys, n); flint_printf("\n");

    /* Lagrange interpolation gives: (to correct)

        w_i = prod_{j != i} (x - x_i) / (x_j - x_i)

        height(w_i) <= height(max height(x_i^(n-1)) / max height(x_i^(2(n-1)))
                    <= (max height(x_i)) * 2(n-1)

        height(f) <= sum_i height(y_i) + height(w_i) * 2(n -1)
                  <= n * ( height(y_i) + height(w_i) * 2(n -1) )

    */

    xbits = _fmpq_vec_max_height(xs, n);
    ybits = _fmpq_vec_max_height(ys, n);
    bound_bits = n * (xbits * 4 * (n - 1) + ybits);
    /* Rational reconstruction needs modulus M >= 1/2*height(f)^2 + 1 */
    bound_bits *= 2;

    /* Max nb of primes and residues needed */
    slong boundprimes = bound_bits / (FLINT_BITS - 1) + 1;
    Lprimes = _nmod_vec_init(boundprimes);
    Lres = _nmod_vec_init(boundprimes * n);
    xm = _nmod_vec_init(n);
    ym = _nmod_vec_init(n);

    total_skipped = 0;
    total_primes = 0;
    /* Important: when changing this, change the list of adversarial primes
       in the test code to match. */
    p = UWORD(1) << (FLINT_BITS - 1);
    flint_printf("(%wd, %wd)\n", xbits, ybits);
    for (;;)
    {
        timeit_start(t1);
        if (total_primes < 16)
        {
            /*todo: divides no den? */
            p = n_nextprime(p, 1);
            //flint_printf("prime = %lu\n", p);
            nmod_init(&mod, p);
            //flint_printf("Reduce:\n");
            for (i = 0; i < n; i++) {
                xm[i] = fmpz_get_nmod(fmpq_numref(xs + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(xs + i), mod);
                xm[i] = nmod_div(xm[i], mden, mod);

                ym[i] = fmpz_get_nmod(fmpq_numref(ys + i), mod);
                mden = fmpz_get_nmod(fmpq_denref(ys + i), mod);
                ym[i] = nmod_div(ym[i], mden, mod);
            }
            //flint_printf("Interpolate\n");
            if (_checked_nmod_poly_interpolate_bis(ym, xm, ym, n, mod))
            {
                num_primes = 1;
                /* CRT */
                //flint_printf("CRT\n");
                if (total_primes == 0)
                {
                    flint_printf("plop\n");
                    // Save this turn
                    for (i = 0; i < n; i++)
                        Lres[i * boundprimes] = ym[i];
                    Lprimes[0] = p;
                    // Initialize tracking coeff and current modulus
                    fmpz_set_ui(curr_res, ym[n-1]);
                    fmpz_set_ui(M, p);
                    flint_printf("plop\n");
                }
                else
                {
                    // Save this turn
                    for (i = 0; i < n; i++)
                        Lres[total_primes + i * boundprimes] = ym[i];
                    Lprimes[total_primes] = p;
                    // Update tracking coeff and current modulus
                    fmpz_CRT_ui(curr_res, curr_res, M, ym[n-1], mod.n, 0);
                    fmpz_mul_ui(M, M, p);
                }
            }
            else
            {
                total_skipped++;
                num_primes = 0;
            }
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

            for (k = 0; k < num_primes; k++)
            {
                /* todo: make sure divides no dens? */
                p = n_nextprime(p, 1);
                primes[k] = p;
            }

            _fmpz_ui_vec_prod(M2, primes, num_primes);

            fmpz_comb_init(comb, primes, num_primes);
            fmpz_comb_temp_init(temp, comb);

            for (j = 0; j < n; j++)
            {
                /* todo: do iterative evaluation may perform better if small values */
               //flint_printf("modular reduction of xs\n");
                fmpz_multi_mod_ui(residues, fmpq_numref(xs + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(xs + j), comb, temp);
                for (k = 0; k < num_primes; k++) {
                    nmod_init(&mod, primes[k]);
                    xmod[k * n + j] = nmod_div(residues[k], residuesden[k], mod);
                }
                //flint_printf("modular reduction of ys\n");
                fmpz_multi_mod_ui(residues, fmpq_numref(ys + j), comb, temp);
                fmpz_multi_mod_ui(residuesden, fmpq_denref(ys + j), comb, temp);
                for (k = 0; k < num_primes; k++) {
                    nmod_init(&mod, primes[k]);
                    //flint_printf("\t%lu %lu %lu\n", residuesden[k], mod.n, primes[k]);
                    ymod[k * n + j] = nmod_div(residues[k], residuesden[k], mod);
                    //flint_printf("%ld %ld %lu\n", j, k, ymod[k*n+j]);
                }
            }

            count_good = 0;
            for (k = 0; k < num_primes; k++)
            {
                nmod_init(&mod, primes[k]);
                good[k] = _checked_nmod_poly_interpolate_bis(ymod + k * n, xmod + k * n, ymod + k * n, n, mod);
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

                total_skipped += (num_primes - count_good);

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
                //flint_printf("CRT\n");
                fmpz_mul(M1M2, M, M2);
                fmpz_mod(c, M, M2);
                fmpz_invmod(c, c, M2);

                for (j = 0; j < n; j++) {
                    i = total_primes + j * boundprimes;
                    for (k = 0; k < num_primes; k++) {
                        Lres[i + k] = ymod[k * n + j];
                        //flint_printf("(%ld, %ld)  %lu\n", i+k, k * n + j, ymod[k * n + j]);
                    }
                }
                for (k = 0; k < num_primes; k++)
                    Lprimes[total_primes + k] = primes[k];

                /* Update tracking coeff and current modulus */
                // Get the slice with the new residues
                tmp = Lres + (n - 1) * boundprimes + total_primes;
                // reconstruct modulo M2 = prod(primes)
                fmpz_multi_CRT_ui(curr_res, tmp, comb, temp, 0);
                // reconstruct modulo M1M2 = M*M2
                _fmpz_CRT(curr_res, curr_res, M, t, M2, M1M2, c, 0);
                fmpz_set(M, M1M2);
            }

            fmpz_comb_temp_clear(temp);
            fmpz_comb_clear(comb);
        }

        total_primes += num_primes;
        timeit_stop(t1);
        flint_printf("interp + CRT: cpu = %ld ms wall = %ld ms\n", t1->cpu, t1->wall);

        /* At that step we have in poly, the polynomial modulo M,
           which a product of total_primes primes */
        //flint_printf("Reconstruct\n");
        //_fmpz_vec_print(poly, n);flint_printf("\n modulo= ");
        //fmpz_print(M); flint_printf("\n");
        /* Rational reconstruction */
        /* todo: handle failure */
        rat_rec = 1;
        timeit_start(t0);
        //rat_rec = fmpz_equal(curr_res, old_res);
        //fmpz_set(old_res, curr_res); // Residue stabilization? Is that relevant?
        //flint_printf(rat_rec ? "RStab\n" : "Not Rstab\n");
        //fmpz_print(old_res); flint_printf("\n");
        if (rat_rec) {
            rat_rec = fmpq_reconstruct_fmpz(curr_coeff, curr_res, M);
            if (rat_rec && total_primes > 1) {
                rat_rec = fmpq_equal(curr_coeff, old_coeffs); // Rational stabilization
                flint_printf(rat_rec ? "Stab\n" : "Not stab\n");
                fmpq_set(old_coeffs, curr_coeff);
                fmpq_print(old_coeffs); flint_printf("\n");
            }
        }
        // todo: wait for stabilization in two steps
        // todo: if ratrec doesn't work use computed CRT (reinitialize Lres and Lprimes)
        // (maybe just do one addition prime)
        if (rat_rec)
        {   // Possible stabilization run complete CRT on each coeff
            // (except the tracker)
            fmpz_comb_t comb;
            fmpz_comb_temp_t temp;
            fmpz_comb_init(comb, Lprimes, total_primes);
            fmpz_comb_temp_init(temp, comb);
            for (i= 0; i < total_primes; i++)
                flint_printf("%lu %lu \n", Lprimes[i], Lres[i+boundprimes]);
            for (i = 0; i < n - 1; i++) {
                fmpz_multi_CRT_ui(poly + i, Lres + i*boundprimes, comb, temp, 0);
                // Then run rational reconstruction
                //fmpz_mul(poly + i, poly + i, fmpq_denref(curr_coeff));
                //fmpz_mod(poly + i, poly + i, M);
                rat_rec = fmpq_reconstruct_fmpz(coeffs + i, poly + i, M);
                if (!rat_rec) break;
                //fmpz_mul(fmpq_denref(coeffs + i), fmpq_denref(coeffs + i), fmpq_denref(curr_coeff));
            }
            timeit_stop(t0);
            flint_printf("Rat rec: stopped after %ld / %ld iteration; cpu = %ld ms wall = %ld ms\n", i + 1, n, t0->cpu, t0->wall);
            fmpz_comb_temp_clear(temp);
            fmpz_comb_clear(comb);
        }

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
        else if (rat_rec)
        {
            //_fmpq_vec_print(coeffs, n); flint_printf("\n");
            /* Check for termination */
            flint_printf("Check");
            timeit_start(t0);
            if (_fmpq_poly_check_interpolant(coeffs, xs, ys, n)) {
                _fmpq_vec_get_fmpz_vec_fmpz(poly, den, coeffs, n);
                timeit_stop(t0);
                flint_printf(" done cpu = %ld ms wall = %ld ms\n", t0->cpu, t0->wall);
                break;
            }
            flint_printf("\n");
            timeit_stop(t0);
        }

        /* No solution */
        if ((slong) fmpz_bits(M) > bound_bits + 1 || total_primes > 100)
        {
            ok = 0;
            break;
        }
        flint_printf("nprime = %wd  height = %wd\n\n", total_primes, fmpz_bits(M));
    }

    flint_printf("total %wd  skipped %wd  bits %wd bound %ld \n", total_primes, total_skipped, fmpz_bits(M), bound_bits);

    fmpz_clear(M);
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(c);
    fmpz_clear(M2);
    fmpz_clear(M1M2);

    _nmod_vec_clear(xm);
    _nmod_vec_clear(ym);

    _fmpz_vec_clear(as, n);
    _fmpz_vec_clear(bs, n);
    _fmpz_vec_clear(cs, n);
    _fmpz_vec_clear(ds, n);

    _fmpq_vec_clear(coeffs, n);

    flint_free(primes);
    flint_free(xmod);
    flint_free(ymod);
    flint_free(residues);
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

