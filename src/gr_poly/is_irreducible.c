/*
    Copyright (C) 2011 Fredrik Johansson
    Copyright (C) 2012 Lina Kulakova
    Copyright (C) 2013 Mike Hansen
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "gr_vec.h"
#include "gr_mat.h"
#include "gr_poly.h"
#include "factor_impl.h"

/*
    Ben-Or's irreducibility test: f of degree n is irreducible iff
    gcd(f, x^(q^i) - x) = 1 for 1 <= i <= n/2. The Frobenius iterates
    are computed via modular composition or exponentiation.
*/
truth_t
gr_poly_is_irreducible_ben_or(const gr_poly_t f, gr_ctx_t ctx)
{
    truth_t res;
    gr_poly_t v, x, xq, xqi, g;
    gr_poly_preinv_t P;
    gr_mat_t H;
    fmpz_t q;
    slong i, n;
    int status = GR_SUCCESS;
    int have_H = 0;

    n = f->length - 1;

    /* zero and constants are not irreducible */
    if (n < 1)
        return (n < 0 || gr_is_zero(f->coeffs, ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;

    if (gr_is_zero(gr_poly_coeff_srcptr(f, n, ctx), ctx) != T_FALSE)
        return T_UNKNOWN;

    if (n == 1)
        return T_TRUE;

    if (_gr_poly_factor_ff_info(NULL, NULL, NULL, ctx) != GR_SUCCESS)
        return T_UNKNOWN;

    res = gr_poly_is_squarefree(f, ctx);
    if (res != T_TRUE)
        return res;

    fmpz_init(q);
    gr_poly_init(v, ctx);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(x, ctx);
    gr_poly_init(xq, ctx);
    gr_poly_init(xqi, ctx);
    gr_poly_init(g, ctx);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_make_monic(v, f, ctx);
    status |= gr_poly_preinv_set(P, v, ctx);
    status |= gr_poly_gen(x, ctx);
    status |= gr_poly_preinv_powmod_x_fmpz(xq, q, P, ctx);
    status |= gr_poly_set(xqi, xq, ctx);

    if (status == GR_SUCCESS && fmpz_bits(q) > ((n_sqrt(n) + 1) * 3) / 4 && n / 2 > 1)
    {
        gr_mat_init(H, n_sqrt(n) + 1, n, ctx);
        have_H = 1;
        status |= gr_poly_preinv_precompute_matrix(H, xq, P, ctx);
    }

    res = T_TRUE;

    for (i = 1; i <= n / 2 && status == GR_SUCCESS; i++)
    {
        status |= gr_poly_sub(g, xqi, x, ctx);
        status |= gr_poly_gcd(g, v, g, ctx);

        if (status != GR_SUCCESS)
            break;

        if (g->length != 1)
        {
            res = (g->length > 1 && gr_is_zero(gr_poly_coeff_srcptr(g, g->length - 1, ctx), ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;
            break;
        }

        if (i == n / 2)
            break;

        if (have_H)
            status |= gr_poly_preinv_compose_mod_brent_kung_precomp(g, xqi, H, P, ctx);
        else
            status |= gr_poly_preinv_powmod_fmpz_sliding(g, xqi, q, 0, P, ctx);

        gr_poly_swap(xqi, g, ctx);
    }

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    fmpz_clear(q);
    gr_poly_clear(v, ctx);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(x, ctx);
    gr_poly_clear(xq, ctx);
    gr_poly_clear(xqi, ctx);
    gr_poly_clear(g, ctx);
    if (have_H)
        gr_mat_clear(H, ctx);

    return res;
}

/*
    Irreducibility test using the coarse phase of the baby-step giant-step
    distinct degree factorization: f is irreducible iff no interval
    polynomial has a nontrivial gcd with f.
*/
truth_t
gr_poly_is_irreducible_ddf(const gr_poly_t f, gr_ctx_t ctx)
{
    truth_t res;
    gr_poly_t v, frob, tmp;
    gr_poly_preinv_t P;
    gr_poly_struct * h, * H, * I;
    gr_mat_t HH;
    fmpz_t q;
    slong i, j, l, m, n, d;
    double beta;
    int status = GR_SUCCESS;
    int have_HH = 0;

    n = f->length - 1;

    /* zero and constants are not irreducible */
    if (n < 1)
        return (n < 0 || gr_is_zero(f->coeffs, ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;

    if (gr_is_zero(gr_poly_coeff_srcptr(f, n, ctx), ctx) != T_FALSE)
        return T_UNKNOWN;

    if (n == 1)
        return T_TRUE;

    if (_gr_poly_factor_ff_info(NULL, NULL, NULL, ctx) != GR_SUCCESS)
        return T_UNKNOWN;

    res = gr_poly_is_squarefree(f, ctx);
    if (res != T_TRUE)
        return res;

    beta = 0.5 * (1. - (log(2) / log(n)));
    l = ceil(pow(n, beta));
    m = ceil(0.5 * n / l);

    fmpz_init(q);
    gr_poly_init(v, ctx);
    gr_poly_preinv_init(P, ctx);
    gr_poly_init(frob, ctx);
    gr_poly_init(tmp, ctx);

    h = flint_malloc((2 * m + l + 1) * sizeof(gr_poly_struct));
    H = h + (l + 1);
    I = H + m;

    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_init(h + i, ctx);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_make_monic(v, f, ctx);
    status |= gr_poly_preinv_set(P, v, ctx);
    status |= gr_poly_preinv_powmod_x_fmpz(frob, q, P, ctx);

    /* baby steps: h[i] = x^{q^i} mod v */
    status |= _gr_poly_iterated_frobenius_preinv(h, l + 1, frob, P, q, ctx);

    if (status != GR_SUCCESS)
    {
        res = T_UNKNOWN;
        goto cleanup;
    }

    status |= gr_poly_set(H + 0, h + l, ctx);
    gr_mat_init(HH, n_sqrt(n) + 1, n, ctx);
    have_HH = 1;
    status |= gr_poly_preinv_precompute_matrix(HH, H + 0, P, ctx);

    res = T_TRUE;
    d = 1;
    for (j = 0; j < m && status == GR_SUCCESS; j++)
    {
        /* giant steps: H[j] = x^{q^(lj)} mod v */
        if (j > 0)
            status |= gr_poly_preinv_compose_mod_brent_kung_precomp(H + j, H + j - 1, HH, P, ctx);

        /* interval polynomial */
        status |= gr_poly_one(I + j, ctx);
        for (i = l - 1; i >= 0 && 2 * d <= n; i--, d++)
        {
            status |= gr_poly_sub(tmp, H + j, h + i, ctx);
            status |= gr_poly_preinv_mulmod(I + j, tmp, I + j, P, ctx);
        }

        status |= gr_poly_gcd(I + j, v, I + j, ctx);

        if (status != GR_SUCCESS)
            break;

        if (I[j].length != 1)
        {
            res = (I[j].length > 1 && gr_is_zero(gr_poly_coeff_srcptr(I + j, I[j].length - 1, ctx), ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;
            break;
        }

        if (n < 2 * d)
            break;
    }

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

cleanup:
    fmpz_clear(q);
    gr_poly_clear(v, ctx);
    gr_poly_preinv_clear(P, ctx);
    gr_poly_clear(frob, ctx);
    gr_poly_clear(tmp, ctx);

    if (have_HH)
        gr_mat_clear(HH, ctx);

    for (i = 0; i < 2 * m + l + 1; i++)
        gr_poly_clear(h + i, ctx);
    flint_free(h);

    return res;
}

/*
    Rabin's irreducibility test: f of degree n is irreducible iff
    x^(q^n) = x mod f and gcd(x^(q^(n/r)) - x, f) = 1 for every prime
    r dividing n. The Frobenius powers x^(q^e) are computed by modular
    composition from the ladder x^(q^(2^i)) (each step a composition of
    the previous with itself), so that the test costs
    O((1 + omega(n)) log n) compositions, independent of q apart from
    the initial x^q. This beats the distinct degree test for irreducible
    inputs of large degree, but does not exit early on inputs with small
    factors; see gr_poly_is_irreducible for the combination.
*/
truth_t
gr_poly_is_irreducible_rabin(const gr_poly_t f, gr_ctx_t ctx)
{
    truth_t res;
    gr_poly_t v, x, t, g;
    gr_poly_preinv_t P;
    gr_poly_struct * ladder;
    fmpz_t q;
    n_factor_t nfac;
    slong i, n, nbits, e;
    int status = GR_SUCCESS;

    n = f->length - 1;

    if (n < 1)
        return (n < 0 || gr_is_zero(f->coeffs, ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;

    if (gr_is_zero(gr_poly_coeff_srcptr(f, n, ctx), ctx) != T_FALSE)
        return T_UNKNOWN;

    if (n == 1)
        return T_TRUE;

    if (_gr_poly_factor_ff_info(NULL, NULL, NULL, ctx) != GR_SUCCESS)
        return T_UNKNOWN;

    res = gr_poly_is_squarefree(f, ctx);
    if (res != T_TRUE)
        return res;

    nbits = FLINT_BIT_COUNT(n);

    fmpz_init(q);
    gr_poly_init(v, ctx);
    gr_poly_init(x, ctx);
    gr_poly_init(t, ctx);
    gr_poly_init(g, ctx);
    gr_poly_preinv_init(P, ctx);
    ladder = flint_malloc(nbits * sizeof(gr_poly_struct));
    for (i = 0; i < nbits; i++)
        gr_poly_init(ladder + i, ctx);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_make_monic(v, f, ctx);
    status |= gr_poly_preinv_set(P, v, ctx);
    status |= gr_poly_gen(x, ctx);

    /* ladder[i] = x^(q^(2^i)) mod v */
    status |= gr_poly_preinv_powmod_x_fmpz(ladder + 0, q, P, ctx);
    for (i = 1; i < nbits && status == GR_SUCCESS; i++)
        status |= gr_poly_preinv_compose_mod(ladder + i, ladder + i - 1, ladder + i - 1, P, ctx);

    res = T_TRUE;

    /* t = x^(q^n) mod v */
#define FROB_POWER(dest, exp) \
    do { \
        slong _e = (exp), _i; \
        int first = 1; \
        for (_i = 0; _e != 0 && status == GR_SUCCESS; _e >>= 1, _i++) \
        { \
            if (_e & 1) \
            { \
                if (first) \
                    status |= gr_poly_set(dest, ladder + _i, ctx); \
                else \
                    status |= gr_poly_preinv_compose_mod(dest, dest, ladder + _i, P, ctx); \
                first = 0; \
            } \
        } \
    } while (0)

    FROB_POWER(t, n);

    if (status == GR_SUCCESS && gr_poly_equal(t, x, ctx) != T_TRUE)
    {
        res = (gr_poly_equal(t, x, ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;
    }
    else if (status == GR_SUCCESS)
    {
        n_factor_init(&nfac);
        n_factor(&nfac, n, 1);

        for (i = 0; i < nfac.num && status == GR_SUCCESS; i++)
        {
            e = n / nfac.p[i];
            FROB_POWER(t, e);
            status |= gr_poly_sub(t, t, x, ctx);
            status |= gr_poly_gcd(g, t, v, ctx);

            if (status == GR_SUCCESS && g->length != 1)
            {
                res = (g->length > 1 && gr_is_zero(gr_poly_coeff_srcptr(g, g->length - 1, ctx), ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;
                break;
            }
        }
    }

#undef FROB_POWER

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    fmpz_clear(q);
    gr_poly_clear(v, ctx);
    gr_poly_clear(x, ctx);
    gr_poly_clear(t, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_preinv_clear(P, ctx);
    for (i = 0; i < nbits; i++)
        gr_poly_clear(ladder + i, ctx);
    flint_free(ladder);

    return res;
}

/* Quick test for linear factors over a small prime field. */
static truth_t
_gr_poly_has_root_small_prime_field(const gr_poly_t f, gr_ctx_t ctx)
{
    fmpz_t p;
    slong deg;
    ulong x, pp;
    gr_ptr t, u;
    truth_t res = T_FALSE;
    int status = GR_SUCCESS;

    fmpz_init(p);

    if (gr_ctx_fq_degree(&deg, ctx) != GR_SUCCESS || deg != 1 ||
        gr_ctx_fq_prime(p, ctx) != GR_SUCCESS ||
        !fmpz_fits_si(p) || fmpz_get_ui(p) > FLINT_MAX(200, 2 * (ulong) f->length))
    {
        fmpz_clear(p);
        return T_UNKNOWN;
    }

    pp = fmpz_get_ui(p);
    fmpz_clear(p);

    GR_TMP_INIT2(t, u, ctx);

    for (x = 0; x < pp; x++)
    {
        status |= gr_set_ui(u, x, ctx);
        status |= gr_poly_evaluate(t, f, u, ctx);

        if (status != GR_SUCCESS)
        {
            res = T_UNKNOWN;
            break;
        }

        res = gr_is_zero(t, ctx);

        if (res != T_FALSE)
            break;
    }

    GR_TMP_CLEAR2(t, u, ctx);

    return res;
}

/* Quick test for irreducible factors of degree at most k: computes
   gcd(f, x^(q^i) - x) for i = 1, ..., k. */
static truth_t
_gr_poly_has_small_factor(const gr_poly_t f, slong k, gr_ctx_t ctx)
{
    gr_poly_t v, x, h, g;
    gr_poly_preinv_t P;
    fmpz_t q;
    slong i;
    truth_t res = T_FALSE;
    int status = GR_SUCCESS;

    fmpz_init(q);
    gr_poly_init(v, ctx);
    gr_poly_init(x, ctx);
    gr_poly_init(h, ctx);
    gr_poly_init(g, ctx);
    gr_poly_preinv_init(P, ctx);

    status |= gr_ctx_fq_order(q, ctx);
    status |= gr_poly_make_monic(v, f, ctx);
    status |= gr_poly_preinv_set(P, v, ctx);
    status |= gr_poly_gen(x, ctx);
    status |= gr_poly_preinv_powmod_x_fmpz(h, q, P, ctx);

    for (i = 1; i <= k && 2 * i <= v->length - 1 && status == GR_SUCCESS; i++)
    {
        if (i > 1)
            status |= gr_poly_preinv_powmod_fmpz_sliding(h, h, q, 0, P, ctx);
        status |= gr_poly_sub(g, h, x, ctx);
        status |= gr_poly_gcd(g, v, g, ctx);

        if (status != GR_SUCCESS)
            break;

        if (g->length != 1)
        {
            res = (g->length > 1 && gr_is_zero(gr_poly_coeff_srcptr(g, g->length - 1, ctx), ctx) == T_FALSE) ? T_TRUE : T_UNKNOWN;
            break;
        }
    }

    if (status != GR_SUCCESS)
        res = T_UNKNOWN;

    fmpz_clear(q);
    gr_poly_clear(v, ctx);
    gr_poly_clear(x, ctx);
    gr_poly_clear(h, ctx);
    gr_poly_clear(g, ctx);
    gr_poly_preinv_clear(P, ctx);

    return res;
}

/* Rabin's test (a few compositions of the Frobenius) beats the distinct
   degree test (about 2 sqrt(n) compositions) for inputs without small
   factors from about this degree on. */
#define GR_POLY_IS_IRREDUCIBLE_RABIN_CUTOFF 600

truth_t
gr_poly_is_irreducible(const gr_poly_t f, gr_ctx_t ctx)
{
    slong n = f->length - 1;
    truth_t res;

    if (n < 1)
        return (n < 0 || gr_is_zero(f->coeffs, ctx) == T_FALSE) ? T_FALSE : T_UNKNOWN;

    if (gr_is_zero(gr_poly_coeff_srcptr(f, n, ctx), ctx) != T_FALSE)
        return T_UNKNOWN;

    if (n == 1)
        return T_TRUE;

    if (_gr_poly_factor_ff_info(NULL, NULL, NULL, ctx) != GR_SUCCESS)
        return T_UNKNOWN;

    /* For small p, trial division by linear factors quickly filters
       out many candidates when testing random polynomials. */
    if (_gr_poly_has_root_small_prime_field(f, ctx) == T_TRUE)
        return T_FALSE;

    if (n < GR_POLY_IS_IRREDUCIBLE_RABIN_CUTOFF)
        return gr_poly_is_irreducible_ddf(f, ctx);

    /* Large degree: filter factors of small degree (the likely
       obstruction for random inputs), then Rabin's test. */
    res = gr_poly_is_squarefree(f, ctx);
    if (res != T_TRUE)
        return res;

    res = _gr_poly_has_small_factor(f, 4, ctx);
    if (res == T_TRUE)
        return T_FALSE;
    if (res == T_UNKNOWN)
        return T_UNKNOWN;

    return gr_poly_is_irreducible_rabin(f, ctx);
}
