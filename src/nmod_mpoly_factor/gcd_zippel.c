/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mat.h"
#include "nmod_mpoly_factor.h"

/* set up mock to point to the coefficients of A, which are not owned by mock */
static void nmod_mpoly_mock_eval_coeff(
    n_polyun_t mock,
    const nmod_mpoly_t A,
    const n_polyun_t Aeh_inc,
    const nmod_mpoly_ctx_t ctx)
{
    slong i, k;

    if (mock->alloc < Aeh_inc->length)
    {
        mock->alloc = FLINT_MAX(Aeh_inc->length, mock->alloc + mock->alloc/2);
        mock->coeffs = FLINT_ARRAY_REALLOC(mock->coeffs, mock->alloc,
                                                                n_poly_struct);
    }

    mock->length = Aeh_inc->length;

    k = 0;
    for (i = 0; i < Aeh_inc->length; i++)
    {
        slong l = Aeh_inc->coeffs[i].length;
        mock->coeffs[i].coeffs = A->coeffs + k;
        mock->coeffs[i].alloc = l;
        mock->coeffs[i].length = l;
        k += l;
    }

    FLINT_ASSERT(k == A->length);
}

static void nmod_mpoly_monomial_evals1(
    n_polyun_t E,
    const nmod_mpoly_t A, const ulong * Amarks, slong Amarkslen,
    n_poly_struct * betas,
    slong m,
    const nmod_mpoly_ctx_t ctx)
{
    mpoly1_monomial_evals_nmod(E, A->exps, A->bits, Amarks,
                                Amarkslen, betas, m, ctx->minfo, ctx->mod);
}

static void n_polyu1n_mod_zip_eval_cur_inc_coeff(
    n_poly_t E,
    n_polyun_t Acur,
    const n_polyun_t Ainc,
    const n_polyun_t Acoeff,
    const nmod_t ctx)
{
    slong i;
    mp_limb_t c;

    FLINT_ASSERT(Acur->length > 0);
    FLINT_ASSERT(Acur->length == Ainc->length);
    FLINT_ASSERT(Acur->length == Acoeff->length);

    n_poly_zero(E);

    for (i = 0; i < Acur->length; i++)
    {
        slong this_len = Acur->coeffs[i].length;
        FLINT_ASSERT(this_len == Ainc->coeffs[i].length);
        FLINT_ASSERT(this_len == Acoeff->coeffs[i].length);
        c = _nmod_zip_eval_step(Acur->coeffs[i].coeffs, Ainc->coeffs[i].coeffs,
                                Acoeff->coeffs[i].coeffs, this_len, ctx);
        n_poly_set_coeff(E, Acur->exps[i], c);
    }
}

static int n_poly_add_zip_must_match(
    n_polyun_t Z,
    n_poly_t A,
    slong cur_length)
{
    slong i, ai;
    slong Alen = A->length;
    ulong * Zexps = Z->exps;
    n_poly_struct * Zcoeffs = Z->coeffs;
    mp_limb_t * Acoeffs = A->coeffs;

    ai = Alen - 1;

    for (i = 0; i < Z->length; i++)
    {
        if (ai >= 0 && Zexps[i] == ai)
        {
            /* Z present, A present */
            Zcoeffs[i].coeffs[cur_length] = Acoeffs[ai];
            Zcoeffs[i].length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[ai] == 0);
        }
        else if (ai < 0 || Zexps[i] > ai)
        {
            /* Z present, A missing */
            Zcoeffs[i].coeffs[cur_length] = 0;
            Zcoeffs[i].length = cur_length + 1;
        }
        else
        {
            /* Z missing, A present */
            return 0;
        }
    }

    return ai < 0;
}


static mp_limb_t * nmod_mat_row_ref(nmod_mat_t M, slong i)
{
    return M->rows[i];
}

static void _nmod_vec_mul(mp_limb_t * a, mp_limb_t * b, mp_limb_t * c,
                                                           slong n, nmod_t ctx)
{
    for (n--; n >= 0; n--)
        a[n] = nmod_mul(b[n], c[n], ctx);
}


/*
    Try to set G to the gcd of A and B given the form f of G.

    f = sum_i X^i sum_j c_ij x^e_j

    images are formed from evaluation at x = alpha^(k+1), 0 <= k < l

    Assume all image gcds have the same degree and are

        sum_i X^i g_ik, 0 <= k < l

    Then, for each 0 <= k < l there is a scalar factor s_k such that

        sum_j c_ij alpha^((k+1)e_j) = s_k g_ik  for 0 <= i <= degree(gcd)

*/
int nmod_mpolyl_gcds_zippel(
    nmod_mpoly_t G, const ulong * Gmarks, slong Gmarkslen,
    nmod_mpoly_t A,
    nmod_mpoly_t B,
    slong *perm,
    slong l,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state,
    slong * Gdegbound,
    n_poly_t Amarks,    /* temps */
    n_poly_t Bmarks)
{
    int success;
    int betas_tries_left, underdetermined_left, exceeded_left;
    slong i, s, S, n, cur_zip_image;
    slong Adeg, Bdeg, Gdeg;
    flint_bitcnt_t bits = A->bits;
    n_poly_t Aev, Bev, Gev;
    n_polyun_t Aeh_inc, Aeh_cur, Aeh_coeff_mock;
    n_polyun_t Beh_inc, Beh_cur, Beh_coeff_mock;
    n_polyun_t HG, MG, ZG;
    mp_limb_t * betas;
    n_poly_struct * beta_caches;
    nmod_mat_struct * ML;
    nmod_mat_t MF, Msol, MFtemp, Mwindow;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(1 < var && var <= ctx->minfo->nvars);

    FLINT_ASSERT(bits <= FLINT_BITS);
    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);

    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(G->length > 0);
    FLINT_ASSERT(Gmarkslen > 0);
    FLINT_ASSERT(*Gdegbound == mpoly_degree_si(G->exps, G->length,
                                                         bits, 0, ctx->minfo));
    if (Gmarkslen < 2)
    {
        G->coeffs[0] = 1;
        return Gmarks[1] - Gmarks[0] == 1;
    }

    betas = FLINT_ARRAY_ALLOC(var, mp_limb_t);
    beta_caches = FLINT_ARRAY_ALLOC(3*var, n_poly_struct);
    for (i = 0; i < var; i++)
    {
        n_poly_init(beta_caches + 3*i + 0);
        n_poly_init(beta_caches + 3*i + 1);
        n_poly_init(beta_caches + 3*i + 2);
    }

    ML = FLINT_ARRAY_ALLOC(Gmarkslen, nmod_mat_struct);
    for (i = 0; i < Gmarkslen; i++)
        nmod_mat_init(ML + i, 0, 0, nmod_mpoly_ctx_modulus(ctx));

    nmod_mat_init(MF, 0, l, nmod_mpoly_ctx_modulus(ctx));
    nmod_mat_init(Msol, l, 1, nmod_mpoly_ctx_modulus(ctx));

    n_poly_init(Aev);
    n_poly_init(Bev);
    n_poly_init(Gev);

    n_polyun_init(Aeh_inc);
    n_polyun_init(Aeh_cur);
    n_polyun_init(Beh_inc);
    n_polyun_init(Beh_cur);
    Aeh_coeff_mock->exps = NULL;
    Aeh_coeff_mock->coeffs = NULL;
    Aeh_coeff_mock->length = 0;
    Aeh_coeff_mock->alloc = 0;
    Beh_coeff_mock->exps = NULL;
    Beh_coeff_mock->coeffs = NULL;
    Beh_coeff_mock->length = 0;
    Beh_coeff_mock->alloc = 0;

    n_polyun_init(HG);
    n_polyun_init(MG);
    n_polyun_init(ZG);

    Adeg = nmod_mpoly_degree_si(A, 0, ctx);
    Bdeg = nmod_mpoly_degree_si(B, 0, ctx);

    mpoly1_fill_marks(&Amarks->coeffs, &Amarks->length, &Amarks->alloc,
                                         A->exps, A->length, bits, ctx->minfo);

    mpoly1_fill_marks(&Bmarks->coeffs, &Bmarks->length, &Bmarks->alloc,
                                         B->exps, B->length, bits, ctx->minfo);

    betas_tries_left = 10;
    underdetermined_left = 3;
    exceeded_left = 3;

next_betas:

    if (--betas_tries_left < 0)
    {
        success = -1;
        goto cleanup;
    }

    for (i = 1; i < var; i++)
    {
        betas[i] = 1 + n_randint(state, nmod_mpoly_ctx_modulus(ctx) - 1);
        nmod_pow_cache_start(betas[i], beta_caches + 3*i + 0,
                                       beta_caches + 3*i + 1,
                                       beta_caches + 3*i + 2);
    }

    nmod_mpoly_monomial_evals1(Aeh_inc, A, Amarks->coeffs, Amarks->length,
                                                  beta_caches + 3*1, var, ctx);
    nmod_mpoly_monomial_evals1(Beh_inc, B, Bmarks->coeffs, Bmarks->length,
                                                  beta_caches + 3*1, var, ctx);

    n_polyun_set(Aeh_cur, Aeh_inc);
    n_polyun_set(Beh_cur, Beh_inc);
    nmod_mpoly_mock_eval_coeff(Aeh_coeff_mock, A, Aeh_inc, ctx);
    nmod_mpoly_mock_eval_coeff(Beh_coeff_mock, B, Beh_inc, ctx);

    nmod_mpoly_monomial_evals1(HG, G, Gmarks, Gmarkslen,
                                                  beta_caches + 3*1, var, ctx);
    n_polyun_zip_start(ZG, HG, l);
    n = n_polyun_product_roots(MG, HG, ctx->mod);
    FLINT_ASSERT(n <= l);

    for (cur_zip_image = 0; cur_zip_image < l; cur_zip_image++)
    {
        n_polyu1n_mod_zip_eval_cur_inc_coeff(Aev, Aeh_cur, Aeh_inc,
                                                     Aeh_coeff_mock, ctx->mod);
        n_polyu1n_mod_zip_eval_cur_inc_coeff(Bev, Beh_cur, Beh_inc,
                                                     Beh_coeff_mock, ctx->mod);

        if (n_poly_degree(Aev) != Adeg)
            goto next_betas;
        if (n_poly_degree(Bev) != Bdeg)
            goto next_betas;

        n_poly_mod_gcd(Gev, Aev, Bev, ctx->mod);
        Gdeg = n_poly_degree(Gev);

        if (Gdeg > *Gdegbound)
        {
            if (--exceeded_left >= 0)
                goto next_betas;

            success = -1;
            goto cleanup;
        }

        if (Gdeg < *Gdegbound)
        {
            *Gdegbound = Gdeg;
            success = 0;
            goto cleanup;
        }

        if (!n_poly_add_zip_must_match(ZG, Gev, cur_zip_image))
            goto next_betas;
    }

    nmod_mat_clear(Msol);
    nmod_mat_init(Msol, 1, l, nmod_mpoly_ctx_modulus(ctx));

    s = perm[0];
    if (Gmarks[s + 1] - Gmarks[s] == 1)
    {
        /* monic case */
        mp_limb_t temp = 1;

        for (i = 0; i < l; i++)
        {
            temp = nmod_mul(temp, HG->coeffs[s].coeffs[0], ctx->mod);

            if (ZG->coeffs[s].coeffs[i] == 0)
                goto general_case;

            nmod_mat_entry(Msol, 0, i) = nmod_div(temp,
                                            ZG->coeffs[s].coeffs[i], ctx->mod);
        }

        goto try_it;
    }

general_case:

    nmod_mat_clear(Msol);
    nmod_mat_init(Msol, 0, l, nmod_mpoly_ctx_modulus(ctx));

    nmod_mat_clear(MF);
    nmod_mat_init(MF, 0, l, nmod_mpoly_ctx_modulus(ctx));

    for (S = 0; S < Gmarkslen; S++)
    {
        s = perm[S];
        n = Gmarks[s + 1] - Gmarks[s];

        FLINT_ASSERT(n <= l);

        if (nmod_mat_nrows(ML + s) != l ||
            nmod_mat_ncols(ML + s) != l + n)
        {
            nmod_mat_clear(ML + s);
            nmod_mat_init(ML + s, l, l + n, nmod_mpoly_ctx_modulus(ctx));
        }

        _nmod_vec_set(nmod_mat_row_ref(ML + s, 0), HG->coeffs[s].coeffs, n);
        for (i = 1; i < l; i++)
        {
            _nmod_vec_mul(nmod_mat_row_ref(ML + s, i),
                          nmod_mat_row_ref(ML + s, i - 1),
                                      HG->coeffs[s].coeffs, n, ctx->mod);
        }

        for (i = 0; i < l; i++)
        {
            _nmod_vec_zero(nmod_mat_row_ref(ML + s, i) + n, l);
            nmod_mat_entry(ML + s, i, n + i) = ZG->coeffs[s].coeffs[i];
        }

        /*
            l x (n + l) matrix ML[s]

                     n                    l
            | monomial evals 1 | image coeffs      0 |
          l |     ...          |         on the      |
            | monomial evals l | 0          diagonal |

        */

        nmod_mat_rref(ML + s);

        for (i = 0; i < n; i++)
            if (nmod_mat_entry(ML + s, i, i) != 1)
                goto next_betas;

        /* delete zero rows from MF, matrix interface makes this fun */
        i = nmod_mat_nrows(MF);
        while (i > 1 && _nmod_vec_is_zero(nmod_mat_row_ref(MF, i - 1), l))
            i--;

        if (i < nmod_mat_nrows(MF))
        {
            nmod_mat_window_init(Mwindow, MF, 0, 0, i, l);
            nmod_mat_init(MFtemp, i, l, nmod_mpoly_ctx_modulus(ctx));
            nmod_mat_set(MFtemp, Mwindow);
            nmod_mat_swap(MF, MFtemp);
            nmod_mat_clear(MFtemp);
            nmod_mat_window_clear(Mwindow);
        }

        /* appends rows to MF */
        nmod_mat_window_init(Mwindow, ML + s, n, n, l, n + l);
        nmod_mat_init(MFtemp, i + l - n, l, nmod_mpoly_ctx_modulus(ctx));
        nmod_mat_concat_vertical(MFtemp, MF, Mwindow);
        nmod_mat_swap(MF, MFtemp);
        nmod_mat_clear(MFtemp);
        nmod_mat_window_clear(Mwindow);

        nmod_mat_clear(Msol);
        nmod_mat_init_nullspace_tr(Msol, MF);

        if (nmod_mat_nrows(Msol) < 1)
        {
            success = 0;
            goto cleanup;
        }

        if (nmod_mat_nrows(Msol) == 1)
            break;
    }

    if (nmod_mat_nrows(Msol) != 1)
    {
        if (--underdetermined_left >= 0)
            goto next_betas;

        success = -1;
        goto cleanup;
    }

try_it:

    /* the first row of Msol has the scales */

    for (s = 0; s < Gmarkslen; s++)
    {
        _nmod_vec_mul(ZG->coeffs[s].coeffs, ZG->coeffs[s].coeffs,
                                       nmod_mat_row_ref(Msol, 0), l, ctx->mod);
    }

    success = n_polyun_zip_solve(G, ZG, HG, MG, ctx);

    if (success < 0)
        goto next_betas;

cleanup:

    n_poly_clear(Aev);
    n_poly_clear(Bev);
    n_poly_clear(Gev);

    for (i = 0; i < var; i++)
    {
        n_poly_clear(beta_caches + 3*i + 0);
        n_poly_clear(beta_caches + 3*i + 1);
        n_poly_clear(beta_caches + 3*i + 2);
    }
    flint_free(beta_caches);
    flint_free(betas);

    for (i = 0; i < Gmarkslen; i++)
        nmod_mat_clear(ML + i);
    flint_free(ML);

    nmod_mat_clear(MF);
    nmod_mat_clear(Msol);

    n_polyun_clear(Aeh_inc);
    n_polyun_clear(Aeh_cur);
    n_polyun_clear(Beh_inc);
    n_polyun_clear(Beh_cur);
    flint_free(Aeh_coeff_mock->coeffs);
    flint_free(Beh_coeff_mock->coeffs);

    n_polyun_clear(HG);
    n_polyun_clear(MG);
    n_polyun_clear(ZG);

    return success;
}


static int _do_bivar_or_univar(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    nmod_mpoly_t A,
    nmod_mpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    if (var == 1)
    {
        int success;
        n_poly_t c;
        n_polyun_t Aev, Bev, Gev, Abarev, Bbarev;
        n_poly_polyun_stack_t St;

        n_poly_stack_init(St->poly_stack);
        n_polyun_stack_init(St->polyun_stack);
        n_polyun_init(Aev);
        n_polyun_init(Bev);
        n_polyun_init(Gev);
        n_polyun_init(Abarev);
        n_polyun_init(Bbarev);
        n_poly_init(c);

        nmod_mpoly_get_polyu1n(Aev, A, 0, 1, ctx);
        nmod_mpoly_get_polyu1n(Bev, B, 0, 1, ctx);

        success = n_polyu1n_mod_gcd_brown_smprime(Gev, Abarev, Bbarev,
                                                       Aev, Bev, ctx->mod, St);
        if (success)
        {
            _n_poly_vec_mod_content(c, Gev->coeffs, Gev->length, ctx->mod);
            success = n_poly_is_one(c);
            nmod_mpoly_set_polyu1n(G, Gev, 0, 1, ctx);
        }

        n_poly_clear(c);
        n_polyun_clear(Aev);
        n_polyun_clear(Bev);
        n_polyun_clear(Gev);
        n_polyun_clear(Abarev);
        n_polyun_clear(Bbarev);
        n_poly_stack_clear(St->poly_stack);
        n_polyun_stack_clear(St->polyun_stack);

        return success;
    }
    else
    {
        n_poly_t a, b, c;
        n_poly_init(a);
        n_poly_init(b);
        n_poly_init(c);
        nmod_mpoly_get_n_poly(a, A, 0, ctx);
        nmod_mpoly_get_n_poly(b, B, 0, ctx);
        n_poly_mod_gcd(c, a, b, ctx->mod);
        _nmod_mpoly_set_nmod_poly(G, G->bits, c->coeffs, c->length, 0, ctx);
        n_poly_clear(a);
        n_poly_clear(b);
        n_poly_clear(c);
        return 1;
    }
}


/*
    A and B depend only on x_0, ..., x_var.
    Assume content_x0(G) = 1, and fail otherwise

    when var = 1, this is checked explicitly.

    when var = 2,

        The possibility of f(x2) | G(x0,x1,x2) for nonconstant f(x2) is ruled
        out by the explicit check.

        suppose f(x1,x2) | G(x0,x1,x2) where f is not in Fp[x2]

        the evaluation point x2 = a2 is chosen so that lc_{x0,x1}(G(x0,x1,x2))
        does not vanish at x2 = a2. This means that lc_x1(f(x1, x2)) does
        not vanish at x2 = a2, and so f(x1, a2) is a nontrivial divisor of
        G(x0,x1,a2) which would be caught by the recursive call.

    ...
*/
int nmod_mpolyl_gcdp_zippel_smprime(
    nmod_mpoly_t G,
    nmod_mpoly_t Abar,
    nmod_mpoly_t Bbar,
    nmod_mpoly_t A,
    nmod_mpoly_t B,
    slong var,
    const nmod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    flint_bitcnt_t bits = A->bits;
    slong i, j, N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong Adeg, Bdeg, Alastdeg, Blastdeg, Gdeg;
    slong bound, Gdegbound, lastdeg, req_zip_images;
    int success, changed, have_enough;
    mp_limb_t alpha, start_alpha, gammaeval, temp;
    n_poly_t a, b, c, gamma, modulus, alphapow;
    nmod_mpoly_t Ac, Bc, Aeval, Beval, Geval, Abareval, Bbareval;
    nmod_mpolyn_t H, T;
    n_poly_t Amarks, Bmarks, Gmarks;
    slong * perm = NULL;

    FLINT_ASSERT(ctx->minfo->ord == ORD_LEX);
    FLINT_ASSERT(0 <= var && var < ctx->minfo->nvars);

    FLINT_ASSERT(A->bits == bits);
    FLINT_ASSERT(B->bits == bits);
    FLINT_ASSERT(G->bits == bits);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);

    if (var < 2)
        return _do_bivar_or_univar(G, Abar, Bbar, A, B, var, ctx, state);

    nmod_mpoly_init3(Ac, 0, bits, ctx);
    nmod_mpoly_init3(Bc, 0, bits, ctx);
    nmod_mpoly_init3(Aeval, 0, bits, ctx);
    nmod_mpoly_init3(Beval, 0, bits, ctx);
    nmod_mpoly_init3(Geval, 0, bits, ctx);
    nmod_mpoly_init3(Abareval, 0, bits, ctx);
    nmod_mpoly_init3(Bbareval, 0, bits, ctx);
    n_poly_init(a);
    n_poly_init(b);
    n_poly_init(c);
    n_poly_init(gamma);
    n_poly_init(modulus);
    n_poly_init2(alphapow, 3);
    nmod_mpolyn_init(H, bits, ctx);
    nmod_mpolyn_init(T, bits, ctx);
    n_poly_init(Amarks);
    n_poly_init(Bmarks);
    n_poly_init(Gmarks);

    nmod_mpolyl_content(Ac, A, var, ctx);
    if (!nmod_mpoly_is_one(Ac, ctx))
    {
        success = nmod_mpoly_divides(A, A, Ac, ctx);
        FLINT_ASSERT(success);
    }

    nmod_mpolyl_content(Bc, B, var, ctx);
    if (!nmod_mpoly_is_one(Bc, ctx))
    {
        success = nmod_mpoly_divides(B, B, Bc, ctx);
        FLINT_ASSERT(success);
    }

    nmod_mpoly_get_n_poly(a, Ac, var, ctx);
    nmod_mpoly_get_n_poly(b, Bc, var, ctx);
    n_poly_mod_gcd(c, a, b, ctx->mod);
    success = n_poly_is_one(c);
    if (!success)
        goto cleanup;

    nmod_mpolyl_lead_coeff(Ac, A, var, ctx);
    nmod_mpolyl_lead_coeff(Bc, B, var, ctx);
    nmod_mpoly_get_n_poly(a, Ac, var, ctx);
    nmod_mpoly_get_n_poly(b, Bc, var, ctx);
    n_poly_mod_gcd(gamma, a, b, ctx->mod);

    /* degree bound on the gcd */
    Adeg = nmod_mpoly_degree_si(A, 0, ctx);
    Bdeg = nmod_mpoly_degree_si(B, 0, ctx);
    Gdegbound = FLINT_MIN(Adeg, Bdeg);

    /* bound of the number of images required */
    Alastdeg = nmod_mpoly_degree_si(A, var, ctx);
    Blastdeg = nmod_mpoly_degree_si(B, var, ctx);
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg) + n_poly_degree(gamma);

    n_poly_one(modulus);

    start_alpha = 1 + n_randint(state, ctx->mod.n - 1);
    alpha = start_alpha;

outer_loop:

    /* get new evaluation point */
    if (alpha < 2)
        alpha = ctx->mod.n;
    alpha -= 1;
    success = (alpha != start_alpha);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;
    alphapow->length = 2;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    gammaeval = n_poly_mod_eval_pow(gamma, alphapow, ctx->mod);
    if (gammaeval == 0)
        goto outer_loop;

    /* make sure evaluation point does not kill either A or B */
    nmod_mpoly_evaluate_one_ui(Aeval, A, var, alpha, ctx);
    nmod_mpoly_evaluate_one_ui(Beval, B, var, alpha, ctx);
    if (Aeval->length == 0 || Beval->length == 0)
        goto outer_loop;

    nmod_mpoly_repack_bits_inplace(Aeval, bits, ctx);
    nmod_mpoly_repack_bits_inplace(Beval, bits, ctx);

    success = nmod_mpolyl_gcdp_zippel_smprime(Geval, Abareval, Bbareval,
                                            Aeval, Beval, var - 1, ctx, state);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(Geval->coeffs[0] == 1);

    if (nmod_mpoly_is_one(Geval, ctx))
    {
        nmod_mpoly_one(G, ctx);
        nmod_mpoly_set(Abar, A, ctx);
        nmod_mpoly_set(Bbar, B, ctx);
        success = 1;
        goto cleanup;
    }

    Gdeg = nmod_mpoly_degree_si(Geval, 0, ctx);

    if (Gdeg > Gdegbound)
        goto outer_loop;

    if (Gdeg < Gdegbound)
        n_poly_one(modulus);

    Gdegbound = Gdeg;

    /* update interpolant H */

    nmod_mpoly_scalar_mul_nmod_invertible(Geval, Geval, gammaeval, ctx);

    if (n_poly_degree(modulus) > 0)
    {
        temp = n_poly_mod_eval_pow(modulus, alphapow, ctx->mod);
        temp = nmod_inv(temp, ctx->mod);
        _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx->mod);
        changed = nmod_mpolyn_interp_crt_sm_mpoly(&lastdeg, H, T, Geval,
                                                       modulus, alphapow->coeffs[1], ctx);
        if (!changed)
        {
            _n_poly_vec_mod_remove_content(c, H->coeffs, H->length, ctx->mod);
            nmod_mpoly_cvtfrom_mpolyn(G, H, var, ctx);
            nmod_mpoly_make_monic(G, G, ctx);
            success = nmod_mpoly_divides(Abar, A, G, ctx) &&
                      nmod_mpoly_divides(Bbar, B, G, ctx);
            if (success)
                goto cleanup;
            /* restore H */
            _n_poly_vec_mod_mul_poly(H->coeffs, H->length, c, ctx->mod);
            goto outer_loop;
        }
    }
    else
    {
        nmod_mpolyn_interp_lift_sm_mpoly(H, Geval, ctx);
    }

    temp = nmod_neg(alpha, ctx->mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 1, temp, ctx->mod);


    nmod_mpoly_fit_length(G, H->length, ctx);
    mpoly_copy_monomials(G->exps, H->exps, H->length, N);
    G->length = H->length;

    mpoly1_fill_marks(&Gmarks->coeffs, &Gmarks->length, &Gmarks->alloc,
                                         G->exps, G->length, bits, ctx->minfo);

    perm = FLINT_ARRAY_REALLOC(perm, Gmarks->length, slong);
    for (i = 0; i < Gmarks->length; i++)
        perm[i] = i;

#define length(k) Gmarks->coeffs[(k)+1] - Gmarks->coeffs[k]

    for (i = 1; i < Gmarks->length; i++)
        for (j = i; j > 0 && length(perm[j-1]) > length(perm[j]); j--)
            FLINT_SWAP(slong, perm[j-1], perm[j]);

    req_zip_images = Gmarks->length - 2;
    j = 0;
    for (i = 0; i < Gmarks->length; i++)
    {
        req_zip_images += length(i);
        j = FLINT_MAX(j, length(i));
    }

    if (Gmarks->length > 1)
        req_zip_images = req_zip_images / (Gmarks->length - 1);

    req_zip_images = FLINT_MAX(req_zip_images, j);
    req_zip_images += 1;

inner_loop:

    /* get new evaluation point */
    if (alpha < 2)
        alpha = ctx->mod.n;
    alpha -= 1;
    success = (alpha != start_alpha);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;
    alphapow->length = 2;

    /* make sure evaluation does not kill both lc(A) and lc(B) */
    gammaeval = n_poly_mod_eval_pow(gamma, alphapow, ctx->mod);
    if (gammaeval == 0)
        goto inner_loop;

    /* make sure evaluation does not kill either A or B */
    nmod_mpoly_evaluate_one_ui(Aeval, A, var, alpha, ctx);
    nmod_mpoly_evaluate_one_ui(Beval, B, var, alpha, ctx);
    if (Aeval->length == 0 || Beval->length == 0)
        goto outer_loop;

    nmod_mpoly_repack_bits_inplace(Aeval, bits, ctx);
    nmod_mpoly_repack_bits_inplace(Beval, bits, ctx);

    success = nmod_mpolyl_gcds_zippel(Geval, Gmarks->coeffs, Gmarks->length,
                                Aeval, Beval, perm, req_zip_images, var, ctx,
                                            state, &Gdegbound, Amarks, Bmarks);
    if (success == 0)
    {
        n_poly_one(modulus);
        goto outer_loop;
    }

    if (success < 0 || Geval->coeffs[0] == 0)
        goto inner_loop;

    /* update interpolant H */
    temp = nmod_div(gammaeval, Geval->coeffs[0], ctx->mod);
    nmod_mpoly_scalar_mul_nmod_invertible(Geval, Geval, temp, ctx);

    FLINT_ASSERT(n_poly_degree(modulus) > 0);

    temp = n_poly_mod_eval_pow(modulus, alphapow, ctx->mod);
    temp = nmod_inv(temp, ctx->mod);
    _n_poly_mod_scalar_mul_nmod_inplace(modulus, temp, ctx->mod);

    changed = nmod_mpolyn_interp_mcrt_sm_mpoly(&lastdeg, H, Geval,
                                                       modulus, alphapow, ctx);
    temp = nmod_neg(alpha, ctx->mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 1, temp, ctx->mod);

    have_enough = n_poly_degree(modulus) >= bound;

    if (changed && !have_enough)
        goto inner_loop;

    if (!changed || have_enough)
    {
        /*
            since H started from a correct modular image and was scaled
            wrt to the univariate gamma, if H is the correct gcd modulo
            content, the only content is univariate in the last variable
        */
        _n_poly_vec_mod_remove_content(c, H->coeffs, H->length, ctx->mod);

        nmod_mpoly_cvtfrom_mpolyn(G, H, var, ctx);
        nmod_mpoly_make_monic(G, G, ctx);

        success = nmod_mpoly_divides(Abar, A, G, ctx) &&
                  nmod_mpoly_divides(Bbar, B, G, ctx);
        if (success)
            goto cleanup;

        /* restore H */
        _n_poly_vec_mod_mul_poly(H->coeffs, H->length, c, ctx->mod);
    }

    if (have_enough)
    {
        n_poly_one(modulus);
        goto outer_loop;
    }

    goto inner_loop;

cleanup:

    flint_free(perm);

    nmod_mpoly_clear(Ac, ctx);
    nmod_mpoly_clear(Bc, ctx);
    nmod_mpoly_clear(Aeval, ctx);
    nmod_mpoly_clear(Beval, ctx);
    nmod_mpoly_clear(Geval, ctx);
    nmod_mpoly_clear(Abareval, ctx);
    nmod_mpoly_clear(Bbareval, ctx);
    n_poly_clear(a);
    n_poly_clear(b);
    n_poly_clear(c);
    n_poly_clear(gamma);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    nmod_mpolyn_clear(H, ctx);
    nmod_mpolyn_clear(T, ctx);
    n_poly_clear(Amarks);
    n_poly_clear(Bmarks);
    n_poly_clear(Gmarks);

    if (success)
    {
        FLINT_ASSERT(G->bits == bits);
    }

    return success;
}
