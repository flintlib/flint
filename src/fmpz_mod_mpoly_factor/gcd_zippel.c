/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_vec.h"
#include "fmpz_mod_vec.h"
#include "fmpz_mod_mat.h"
#include "fmpz_mod_mpoly_factor.h"
#include "n_poly.h"

static void fmpz_mod_mpoly_monomial_evals1(
    fmpz_mod_polyun_t E,
    const fmpz_mod_mpoly_t A, const ulong * Amarks, slong Amarkslen,
    fmpz_mod_poly_struct * betas,
    slong m,
    const fmpz_mod_mpoly_ctx_t ctx) /* temp space */
{
    mpoly1_monomial_evals_fmpz_mod(E, A->exps, A->bits, Amarks,
                                Amarkslen, betas, m, ctx->minfo, ctx->ffinfo);
}


static void fmpz_mod_polyu1n_zip_eval_cur_inc_coeff(
    fmpz_mod_poly_t E,
    fmpz_mod_polyun_t Acur,
    const fmpz_mod_polyun_t Ainc,
    const fmpz_mod_polyun_t Acoeff,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    fmpz_t c;

    FLINT_ASSERT(Acur->length > 0);
    FLINT_ASSERT(Acur->length == Ainc->length);
    FLINT_ASSERT(Acur->length == Acoeff->length);

    fmpz_init(c);

    fmpz_mod_poly_zero(E, ctx);

    for (i = 0; i < Acur->length; i++)
    {
        slong this_len = Acur->coeffs[i].length;
        FLINT_ASSERT(this_len == Ainc->coeffs[i].length);
        FLINT_ASSERT(this_len == Acoeff->coeffs[i].length);
        _fmpz_mod_zip_eval_step(c, Acur->coeffs[i].coeffs,
              Ainc->coeffs[i].coeffs, Acoeff->coeffs[i].coeffs, this_len, ctx);
        fmpz_mod_poly_set_coeff_fmpz(E, Acur->exps[i], c, ctx);
    }

    fmpz_clear(c);
}

static int fmpz_mod_poly_add_zip_must_match(
    fmpz_mod_polyun_t Z,
    fmpz_mod_poly_t A,
    slong cur_length)
{
    slong i, ai;
    slong Alen = A->length;
    ulong * Zexps = Z->exps;
    fmpz_mod_poly_struct * Zcoeffs = Z->coeffs;
    fmpz * Acoeffs = A->coeffs;

    ai = Alen - 1;

    for (i = 0; i < Z->length; i++)
    {
        if (ai >= 0 && Zexps[i] == ai)
        {
            /* Z present, A present */
            fmpz_set(Zcoeffs[i].coeffs + cur_length, Acoeffs + ai);
            Zcoeffs[i].length = cur_length + 1;
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs + ai));
        }
        else if (ai < 0 || Zexps[i] > ai)
        {
            /* Z present, A missing */
            fmpz_zero(Zcoeffs[i].coeffs + cur_length);
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


static fmpz * fmpz_mod_mat_row_ref(fmpz_mod_mat_t M, slong i)
{
    return M->mat->rows[i];
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
int fmpz_mod_mpolyl_gcds_zippel(
    fmpz_mod_mpoly_t G, const ulong * Gmarks, slong Gmarkslen,
    fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_t B,
    slong *perm,
    slong l,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state,
    slong * Gdegbound,
    n_poly_t Amarks,
    n_poly_t Bmarks)
{
    int success;
    int betas_tries_left, underdetermined_left, exceeded_left;
    slong i, s, S, n, cur_zip_image;
    slong Adeg, Bdeg, Gdeg;
    flint_bitcnt_t bits = A->bits;
    fmpz_mod_poly_t Aev, Bev, Gev;
    fmpz_mod_polyun_t Aeh_inc, Aeh_cur, Aeh_coeff_mock;
    fmpz_mod_polyun_t Beh_inc, Beh_cur, Beh_coeff_mock;
    fmpz_mod_polyun_t HG, MG, ZG;
    fmpz * betas;
    fmpz_mod_poly_struct * beta_caches;
    fmpz_mod_mat_struct * ML;
    fmpz_mod_mat_t MF, Msol, MFtemp, Mwindow;

    FLINT_ASSERT(1 < var && var < ctx->minfo->nvars);
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
        fmpz_one(G->coeffs + 0);
        return Gmarks[1] - Gmarks[0] == 1;
    }

    betas = _fmpz_vec_init(var);
    beta_caches = FLINT_ARRAY_ALLOC(var, fmpz_mod_poly_struct);
    for (i = 0; i < var; i++)
        fmpz_mod_poly_init(beta_caches + i, ctx->ffinfo);

    ML = FLINT_ARRAY_ALLOC(Gmarkslen, fmpz_mod_mat_struct);
    for (i = 0; i < Gmarkslen; i++)
        fmpz_mod_mat_init(ML + i, 0, 0, fmpz_mod_ctx_modulus(ctx->ffinfo));

    fmpz_mod_mat_init(MF, 0, l, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_mod_mat_init(Msol, l, 1, fmpz_mod_ctx_modulus(ctx->ffinfo));

    fmpz_mod_poly_init(Aev, ctx->ffinfo);
    fmpz_mod_poly_init(Bev, ctx->ffinfo);
    fmpz_mod_poly_init(Gev, ctx->ffinfo);

    fmpz_mod_polyun_init(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_init(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_inc, ctx->ffinfo);
    fmpz_mod_polyun_init(Beh_cur, ctx->ffinfo);
    Aeh_coeff_mock->exps = NULL;
    Aeh_coeff_mock->coeffs = NULL;
    Aeh_coeff_mock->length = 0;
    Aeh_coeff_mock->alloc = 0;
    Beh_coeff_mock->exps = NULL;
    Beh_coeff_mock->coeffs = NULL;
    Beh_coeff_mock->length = 0;
    Beh_coeff_mock->alloc = 0;

    fmpz_mod_polyun_init(HG, ctx->ffinfo);
    fmpz_mod_polyun_init(MG, ctx->ffinfo);
    fmpz_mod_polyun_init(ZG, ctx->ffinfo);

    Adeg = fmpz_mod_mpoly_degree_si(A, 0, ctx);
    Bdeg = fmpz_mod_mpoly_degree_si(B, 0, ctx);

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
        fmpz_mod_rand_not_zero(betas + i, state, ctx->ffinfo);
        fmpz_mod_pow_cache_start(betas + i, beta_caches + i, ctx->ffinfo);
    }

    fmpz_mod_mpoly_monomial_evals1(Aeh_inc, A, Amarks->coeffs, Amarks->length, beta_caches + 1, var, ctx);
    fmpz_mod_mpoly_monomial_evals1(Beh_inc, B, Bmarks->coeffs, Bmarks->length, beta_caches + 1, var, ctx);

    fmpz_mod_polyun_set(Aeh_cur, Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_set(Beh_cur, Beh_inc, ctx->ffinfo);
    fmpz_mod_mpoly_mock_eval_coeff(Aeh_coeff_mock, A, Aeh_inc, ctx);
    fmpz_mod_mpoly_mock_eval_coeff(Beh_coeff_mock, B, Beh_inc, ctx);

    fmpz_mod_mpoly_monomial_evals1(HG, G, Gmarks, Gmarkslen, beta_caches + 1, var, ctx);

    fmpz_mod_polyun_zip_start(ZG, HG, l, ctx->ffinfo);
    n = fmpz_mod_polyun_product_roots(MG, HG, ctx->ffinfo);
    FLINT_ASSERT(n <= l);

    for (cur_zip_image = 0; cur_zip_image < l; cur_zip_image++)
    {
        fmpz_mod_polyu1n_zip_eval_cur_inc_coeff(Aev, Aeh_cur, Aeh_inc, Aeh_coeff_mock, ctx->ffinfo);
        fmpz_mod_polyu1n_zip_eval_cur_inc_coeff(Bev, Beh_cur, Beh_inc, Beh_coeff_mock, ctx->ffinfo);

        if (_fmpz_mod_poly_degree(Aev) != Adeg)
            goto next_betas;
        if (_fmpz_mod_poly_degree(Bev) != Bdeg)
            goto next_betas;

        fmpz_mod_poly_gcd(Gev, Aev, Bev, ctx->ffinfo);
        Gdeg = _fmpz_mod_poly_degree(Gev);

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

        if (!fmpz_mod_poly_add_zip_must_match(ZG, Gev, cur_zip_image))
            goto next_betas;
    }

    fmpz_mod_mat_clear(Msol);
    fmpz_mod_mat_init(Msol, 1, l, fmpz_mod_ctx_modulus(ctx->ffinfo));

    s = perm[0];
    if (Gmarks[s + 1] - Gmarks[s] == 1)
    {
        /* monic case */
        fmpz * temp = betas; /* use betas[0] for tmp space */

        fmpz_one(temp);
        for (i = 0; i < l; i++)
        {
            fmpz_mod_mul(temp, temp, HG->coeffs[s].coeffs + 0, ctx->ffinfo);
            if (!fmpz_mod_divides(fmpz_mod_mat_entry(Msol, 0, i), temp,
                                        ZG->coeffs[s].coeffs + i, ctx->ffinfo))
            {
                goto general_case;
            }
        }

        goto try_it;
    }

general_case:

    fmpz_mod_mat_clear(Msol);
    fmpz_mod_mat_init(Msol, 0, l, fmpz_mod_ctx_modulus(ctx->ffinfo));

    fmpz_mod_mat_clear(MF);
    fmpz_mod_mat_init(MF, 0, l, fmpz_mod_ctx_modulus(ctx->ffinfo));

    for (S = 0; S < Gmarkslen; S++)
    {
        s = perm[S];
        n = Gmarks[s + 1] - Gmarks[s];

        FLINT_ASSERT(n <= l);

        if (fmpz_mod_mat_nrows(ML + s) != l ||
            fmpz_mod_mat_ncols(ML + s) != l + n)
        {
            fmpz_mod_mat_clear(ML + s);
            fmpz_mod_mat_init(ML + s, l, l + n, fmpz_mod_ctx_modulus(ctx->ffinfo));
        }

        _fmpz_vec_set(fmpz_mod_mat_row_ref(ML + s, 0), HG->coeffs[s].coeffs, n);
        for (i = 1; i < l; i++)
        {
            _fmpz_mod_vec_mul(fmpz_mod_mat_row_ref(ML + s, i),
                              fmpz_mod_mat_row_ref(ML + s, i - 1),
                                         HG->coeffs[s].coeffs, n, ctx->ffinfo);
        }

        for (i = 0; i < l; i++)
        {
            _fmpz_vec_zero(fmpz_mod_mat_row_ref(ML + s, i) + n, l);
            fmpz_set(fmpz_mod_mat_row_ref(ML + s, i) + n + i, ZG->coeffs[s].coeffs + i);
        }

        /*
            l x (n + l) matrix ML[s]

                     n                    l
            | monomial evals 1 | image coeffs      0 |
          l |     ...          |         on the      |
            | monomial evals l | 0          diagonal |

        */

        fmpz_mod_mat_rref(NULL, ML + s);

        for (i = 0; i < n; i++)
            if (!fmpz_is_one(fmpz_mod_mat_entry(ML + s, i, i)))
                goto next_betas;

        /* delete zero rows from MF, matrix interface makes this fun */
        i = fmpz_mod_mat_nrows(MF);
        while (i > 1 && _fmpz_vec_is_zero(fmpz_mod_mat_row_ref(MF, i - 1), l))
            i--;

        if (i < fmpz_mod_mat_nrows(MF))
        {
            fmpz_mod_mat_window_init(Mwindow, MF, 0, 0, i, l);
            fmpz_mod_mat_init(MFtemp, i, l, fmpz_mod_ctx_modulus(ctx->ffinfo));
            fmpz_mod_mat_set(MFtemp, Mwindow);
            fmpz_mod_mat_swap(MF, MFtemp);
            fmpz_mod_mat_clear(MFtemp);
            fmpz_mod_mat_window_clear(Mwindow);
        }

        /* appends rows to MF */
        fmpz_mod_mat_window_init(Mwindow, ML + s, n, n, l, n + l);
        fmpz_mod_mat_init(MFtemp, i + l - n, l, fmpz_mod_ctx_modulus(ctx->ffinfo));
        fmpz_mod_mat_concat_vertical(MFtemp, MF, Mwindow);
        fmpz_mod_mat_swap(MF, MFtemp);
        fmpz_mod_mat_clear(MFtemp);
        fmpz_mod_mat_window_clear(Mwindow);

        fmpz_mod_mat_clear(Msol);
        fmpz_mod_mat_init_nullspace_tr(Msol, MF, ctx->ffinfo);

        if (fmpz_mod_mat_nrows(Msol) < 1)
        {
            success = 0;
            goto cleanup;
        }

        if (fmpz_mod_mat_nrows(Msol) == 1)
            break;
    }

    if (fmpz_mod_mat_nrows(Msol) != 1)
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
        _fmpz_mod_vec_mul(ZG->coeffs[s].coeffs, ZG->coeffs[s].coeffs,
                                fmpz_mod_mat_row_ref(Msol, 0), l, ctx->ffinfo);
    }

    success = fmpz_mod_polyun_zip_solve(G, ZG, HG, MG, ctx);

    if (success < 0)
        goto next_betas;

cleanup:

    fmpz_mod_poly_clear(Aev, ctx->ffinfo);
    fmpz_mod_poly_clear(Bev, ctx->ffinfo);
    fmpz_mod_poly_clear(Gev, ctx->ffinfo);

    _fmpz_vec_clear(betas, var);
    for (i = 0; i < var; i++)
        fmpz_mod_poly_clear(beta_caches + i, ctx->ffinfo);
    flint_free(beta_caches);

    for (i = 0; i < Gmarkslen; i++)
        fmpz_mod_mat_clear(ML + i);
    flint_free(ML);

    fmpz_mod_mat_clear(MF);
    fmpz_mod_mat_clear(Msol);

    fmpz_mod_polyun_clear(Aeh_inc, ctx->ffinfo);
    fmpz_mod_polyun_clear(Aeh_cur, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_inc, ctx->ffinfo);
    fmpz_mod_polyun_clear(Beh_cur, ctx->ffinfo);
    flint_free(Aeh_coeff_mock->coeffs);
    flint_free(Beh_coeff_mock->coeffs);

    fmpz_mod_polyun_clear(HG, ctx->ffinfo);
    fmpz_mod_polyun_clear(MG, ctx->ffinfo);
    fmpz_mod_polyun_clear(ZG, ctx->ffinfo);

    return success;
}


static int _do_bivar_or_univar(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    if (var == 1)
    {
        int success;
        fmpz_mod_poly_t c;
        fmpz_mod_polyun_t Aev, Bev, Gev, Abarev, Bbarev;
        fmpz_mod_poly_polyun_stack_t St;

        fmpz_mod_poly_stack_init(St->poly_stack);
        fmpz_mod_polyun_stack_init(St->polyun_stack);
        fmpz_mod_polyun_init(Aev, ctx->ffinfo);
        fmpz_mod_polyun_init(Bev, ctx->ffinfo);
        fmpz_mod_polyun_init(Gev, ctx->ffinfo);
        fmpz_mod_polyun_init(Abarev, ctx->ffinfo);
        fmpz_mod_polyun_init(Bbarev, ctx->ffinfo);
        fmpz_mod_poly_init(c, ctx->ffinfo);

        fmpz_mod_mpoly_get_polyu1n(Aev, A, 0, 1, ctx);
        fmpz_mod_mpoly_get_polyu1n(Bev, B, 0, 1, ctx);

        success = fmpz_mod_polyu1n_gcd_brown_smprime(Gev, Abarev, Bbarev,
                      Aev, Bev, ctx->ffinfo, St->poly_stack, St->polyun_stack);

        if (success)
        {
            _fmpz_mod_poly_vec_content(c, Gev->coeffs, Gev->length, ctx->ffinfo);
            success = fmpz_mod_poly_is_one(c, ctx->ffinfo);
            fmpz_mod_mpoly_set_polyu1n(G, Gev, 0, 1, ctx);
        }

        fmpz_mod_poly_clear(c, ctx->ffinfo);
        fmpz_mod_polyun_clear(Aev, ctx->ffinfo);
        fmpz_mod_polyun_clear(Bev, ctx->ffinfo);
        fmpz_mod_polyun_clear(Gev, ctx->ffinfo);
        fmpz_mod_polyun_clear(Abarev, ctx->ffinfo);
        fmpz_mod_polyun_clear(Bbarev, ctx->ffinfo);
        fmpz_mod_poly_stack_clear(St->poly_stack);
        fmpz_mod_polyun_stack_clear(St->polyun_stack);

        return success;
    }
    else
    {
        fmpz_mod_poly_t a, b, c;
        fmpz_mod_poly_init(a, ctx->ffinfo);
        fmpz_mod_poly_init(b, ctx->ffinfo);
        fmpz_mod_poly_init(c, ctx->ffinfo);
        fmpz_mod_mpoly_get_fmpz_mod_poly(a, A, 0, ctx);
        fmpz_mod_mpoly_get_fmpz_mod_poly(b, B, 0, ctx);
        fmpz_mod_poly_gcd(c, a, b, ctx->ffinfo);
        _fmpz_mod_mpoly_set_fmpz_mod_poly(G, G->bits, c->coeffs, c->length, 0, ctx);
        fmpz_mod_poly_clear(a, ctx->ffinfo);
        fmpz_mod_poly_clear(b, ctx->ffinfo);
        fmpz_mod_poly_clear(c, ctx->ffinfo);
        return 1;
    }
}

/*
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

    when var = 3,

        The possibility of f(x3) | G(x0,x1,x2,x3) for nonconstant f(x3) is ruled
        out by the explicit check.

        suppose f(x1,x2,x3) | G(x0,x1,x2,x3) where f is not in Fp[x3]

        the evaluation point x3 = a3 is chosen so that lc_{x0,x1,X2}(G(x0,x1,x2,X3))
        does not vanish at x3 = a3. This means that G(x0, x1, x2, a3) has
        content that would be caught be the recursive call.

    ...
*/

int fmpz_mod_mpolyl_gcdp_zippel(
    fmpz_mod_mpoly_t G,
    fmpz_mod_mpoly_t Abar,
    fmpz_mod_mpoly_t Bbar,
    fmpz_mod_mpoly_t A,
    fmpz_mod_mpoly_t B,
    slong var,
    const fmpz_mod_mpoly_ctx_t ctx,
    flint_rand_t state)
{
    flint_bitcnt_t bits = A->bits;
    slong i, j, N = mpoly_words_per_exp_sp(bits, ctx->minfo);
    slong Adeg, Bdeg, Alastdeg, Blastdeg, Gdeg;
    slong bound, Gdegbound, lastdeg, req_zip_images;
    int success, changed, have_enough;
    fmpz_t alpha, start_alpha, gammaeval, temp;
    fmpz_mod_poly_t a, b, c, gamma, modulus, alphapow;
    fmpz_mod_mpoly_t Ac, Bc, Aeval, Beval, Geval, Abareval, Bbareval;
    fmpz_mod_mpolyn_t H, T;
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

    fmpz_init(alpha);
    fmpz_init(start_alpha);
    fmpz_init(gammaeval);
    fmpz_init(temp);
    fmpz_mod_mpoly_init3(Ac, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Bc, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Aeval, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Beval, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Geval, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Abareval, 0, bits, ctx);
    fmpz_mod_mpoly_init3(Bbareval, 0, bits, ctx);
    fmpz_mod_poly_init(a, ctx->ffinfo);
    fmpz_mod_poly_init(b, ctx->ffinfo);
    fmpz_mod_poly_init(c, ctx->ffinfo);
    fmpz_mod_poly_init(gamma, ctx->ffinfo);
    fmpz_mod_poly_init(modulus, ctx->ffinfo);
    fmpz_mod_poly_init2(alphapow, 3, ctx->ffinfo);
    fmpz_mod_mpolyn_init(H, bits, ctx);
    fmpz_mod_mpolyn_init(T, bits, ctx);
    n_poly_init(Amarks);
    n_poly_init(Bmarks);
    n_poly_init(Gmarks);

    fmpz_mod_mpolyl_content(Ac, A, var, ctx);
    if (!fmpz_mod_mpoly_is_one(Ac, ctx))
    {
        success = fmpz_mod_mpoly_divides(A, A, Ac, ctx);
        FLINT_ASSERT(success);
    }

    fmpz_mod_mpolyl_content(Bc, B, var, ctx);
    if (!fmpz_mod_mpoly_is_one(Bc, ctx))
    {
        success = fmpz_mod_mpoly_divides(B, B, Bc, ctx);
        FLINT_ASSERT(success);
    }

    fmpz_mod_mpoly_get_fmpz_mod_poly(a, Ac, var, ctx);
    fmpz_mod_mpoly_get_fmpz_mod_poly(b, Bc, var, ctx);
    fmpz_mod_poly_gcd(c, a, b, ctx->ffinfo);
    success = fmpz_mod_poly_is_one(c, ctx->ffinfo);
    if (!success)
        goto cleanup;

    fmpz_mod_mpolyl_lead_coeff(Ac, A, var, ctx);
    fmpz_mod_mpolyl_lead_coeff(Bc, B, var, ctx);
    fmpz_mod_mpoly_get_fmpz_mod_poly(a, Ac, var, ctx);
    fmpz_mod_mpoly_get_fmpz_mod_poly(b, Bc, var, ctx);
    fmpz_mod_poly_gcd(gamma, a, b, ctx->ffinfo);

    /* degree bound on the gcd */
    Adeg = fmpz_mod_mpoly_degree_si(A, 0, ctx);
    Bdeg = fmpz_mod_mpoly_degree_si(B, 0, ctx);
    Gdegbound = FLINT_MIN(Adeg, Bdeg);

    /* bound of the number of images required */
    Alastdeg = fmpz_mod_mpoly_degree_si(A, var, ctx);
    Blastdeg = fmpz_mod_mpoly_degree_si(B, var, ctx);
    bound = 1 + FLINT_MIN(Alastdeg, Blastdeg) + _fmpz_mod_poly_degree(gamma);

    fmpz_mod_poly_one(modulus, ctx->ffinfo);

    fmpz_mod_rand_not_zero(start_alpha, state, ctx->ffinfo);
    fmpz_set(alpha, start_alpha);

outer_loop:

    /* get new evaluation point */
    if (fmpz_cmp_ui(alpha, 2) < 0)
        fmpz_set(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_sub_ui(alpha, alpha, 1);
    success = !fmpz_equal(alpha, start_alpha);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(alphapow->alloc >= 2);
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);
    alphapow->length = 2;

    /* make sure evaluation point does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_eval_pow(gammaeval, gamma, alphapow, ctx->ffinfo);
    if (fmpz_is_zero(gammaeval))
        goto outer_loop;

    /* make sure evaluation point does not kill either A or B */
    fmpz_mod_mpoly_evaluate_one_fmpz(Aeval, A, var, alpha, ctx);
    fmpz_mod_mpoly_evaluate_one_fmpz(Beval, B, var, alpha, ctx);
    if (Aeval->length == 0 || Beval->length == 0)
        goto outer_loop;

    fmpz_mod_mpoly_repack_bits_inplace(Aeval, bits, ctx);
    fmpz_mod_mpoly_repack_bits_inplace(Beval, bits, ctx);

    success = fmpz_mod_mpolyl_gcdp_zippel(Geval, Abareval, Bbareval,
                                            Aeval, Beval, var - 1, ctx, state);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(Geval->length > 0);
    FLINT_ASSERT(fmpz_is_one(Geval->coeffs + 0));

    if (fmpz_mod_mpoly_is_one(Geval, ctx))
        goto gcd_is_trivial;

    Gdeg = fmpz_mod_mpoly_degree_si(Geval, 0, ctx);

    if (Gdeg > Gdegbound)
        goto outer_loop;

    if (Gdeg < Gdegbound)
        fmpz_mod_poly_one(modulus, ctx->ffinfo);

    Gdegbound = Gdeg;

    /* update interpolant H */

    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(Geval, Geval, gammaeval, ctx);

    if (_fmpz_mod_poly_degree(modulus) > 0)
    {
        fmpz_mod_poly_eval_pow(temp, modulus, alphapow, ctx->ffinfo);
        fmpz_mod_poly_scalar_div_fmpz(modulus, modulus, temp, ctx->ffinfo);
        changed = fmpz_mod_mpolyn_interp_crt_sm_mpoly(&lastdeg, H, T, Geval,
                                                       modulus, alphapow, ctx);
        if (!changed)
        {
            _fmpz_mod_poly_vec_remove_content(c, H->coeffs, H->length, ctx->ffinfo);
            fmpz_mod_mpoly_cvtfrom_mpolyn(G, H, var, ctx);
            fmpz_mod_mpoly_make_monic(G, G, ctx);
            success = fmpz_mod_mpoly_divides(Abar, A, G, ctx) &&
                      fmpz_mod_mpoly_divides(Bbar, B, G, ctx);
            if (success)
                goto cleanup;
            /* restore H */
            _fmpz_mod_poly_vec_mul_poly(H->coeffs, H->length, c, ctx->ffinfo);
            goto outer_loop;
        }
    }
    else
    {
        fmpz_mod_mpolyn_interp_lift_sm_mpoly(H, Geval, ctx);
    }

    fmpz_mod_neg(temp, alpha, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, temp, ctx->ffinfo);


    fmpz_mod_mpoly_fit_length(G, H->length, ctx);
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
    if (fmpz_cmp_ui(alpha, 2) < 0)
        fmpz_set(alpha, fmpz_mod_ctx_modulus(ctx->ffinfo));
    fmpz_sub_ui(alpha, alpha, 1);
    success = !fmpz_equal(alpha, start_alpha);
    if (!success)
        goto cleanup;

    FLINT_ASSERT(alphapow->alloc >= 2);
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);
    alphapow->length = 2;

    /* make sure evaluation does not kill both lc(A) and lc(B) */
    fmpz_mod_poly_eval_pow(gammaeval, gamma, alphapow, ctx->ffinfo);
    if (fmpz_is_zero(gammaeval))
        goto inner_loop;

    /* make sure evaluation does not kill either A or B */
    fmpz_mod_mpoly_evaluate_one_fmpz(Aeval, A, var, alpha, ctx);
    fmpz_mod_mpoly_evaluate_one_fmpz(Beval, B, var, alpha, ctx);
    if (Aeval->length == 0 || Beval->length == 0)
        goto outer_loop;

    fmpz_mod_mpoly_repack_bits_inplace(Aeval, bits, ctx);
    fmpz_mod_mpoly_repack_bits_inplace(Beval, bits, ctx);

    success = fmpz_mod_mpolyl_gcds_zippel(Geval, Gmarks->coeffs, Gmarks->length,
                                Aeval, Beval, perm, req_zip_images, var, ctx,
                                            state, &Gdegbound, Amarks, Bmarks);
    if (success == 0)
    {
        fmpz_mod_poly_one(modulus, ctx->ffinfo);
        goto outer_loop;
    }

    if (success < 0 || fmpz_is_zero(Geval->coeffs + 0))
        goto inner_loop;

    /* update interpolant H */
    fmpz_mod_divides(temp, gammaeval, Geval->coeffs + 0, ctx->ffinfo);
    fmpz_mod_mpoly_scalar_mul_fmpz_mod_invertible(Geval, Geval, temp, ctx);

    FLINT_ASSERT(_fmpz_mod_poly_degree(modulus) > 0);

    fmpz_mod_poly_eval_pow(temp, modulus, alphapow, ctx->ffinfo);
    fmpz_mod_poly_scalar_div_fmpz(modulus, modulus, temp, ctx->ffinfo);

    changed = fmpz_mod_mpolyn_interp_mcrt_sm_mpoly(&lastdeg, H, Geval,
                                                       modulus, alphapow, ctx);
    fmpz_mod_neg(temp, alpha, ctx->ffinfo);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 1, temp, ctx->ffinfo);

    have_enough = _fmpz_mod_poly_degree(modulus) >= bound;

    if (changed && !have_enough)
        goto inner_loop;

    if (!changed || have_enough)
    {
        /*
            since H started from a correct modular image and was scaled
            wrt to the univariate gamma, if H is the correct gcd modulo
            content, the only content is univariate in the last variable
        */
        _fmpz_mod_poly_vec_remove_content(c, H->coeffs, H->length, ctx->ffinfo);

        fmpz_mod_mpoly_cvtfrom_mpolyn(G, H, var, ctx);
        fmpz_mod_mpoly_make_monic(G, G, ctx);

        success = fmpz_mod_mpoly_divides(Abar, A, G, ctx) &&
                  fmpz_mod_mpoly_divides(Bbar, B, G, ctx);
        if (success)
            goto cleanup;

        /* restore H */
        _fmpz_mod_poly_vec_mul_poly(H->coeffs, H->length, c, ctx->ffinfo);
    }

    if (have_enough)
    {
        fmpz_mod_poly_one(modulus, ctx->ffinfo);
        goto outer_loop;
    }

    goto inner_loop;

cleanup:

    flint_free(perm);

    fmpz_clear(alpha);
    fmpz_clear(start_alpha);
    fmpz_clear(gammaeval);
    fmpz_clear(temp);
    fmpz_mod_mpoly_clear(Ac, ctx);
    fmpz_mod_mpoly_clear(Bc, ctx);
    fmpz_mod_mpoly_clear(Aeval, ctx);
    fmpz_mod_mpoly_clear(Beval, ctx);
    fmpz_mod_mpoly_clear(Geval, ctx);
    fmpz_mod_mpoly_clear(Abareval, ctx);
    fmpz_mod_mpoly_clear(Bbareval, ctx);
    fmpz_mod_poly_clear(a, ctx->ffinfo);
    fmpz_mod_poly_clear(b, ctx->ffinfo);
    fmpz_mod_poly_clear(c, ctx->ffinfo);
    fmpz_mod_poly_clear(gamma, ctx->ffinfo);
    fmpz_mod_poly_clear(modulus, ctx->ffinfo);
    fmpz_mod_poly_clear(alphapow, ctx->ffinfo);
    fmpz_mod_mpolyn_clear(H, ctx);
    fmpz_mod_mpolyn_clear(T, ctx);
    n_poly_clear(Amarks);
    n_poly_clear(Bmarks);
    n_poly_clear(Gmarks);

    FLINT_ASSERT(success || fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx->ffinfo)));

    return success;

gcd_is_trivial:

    fmpz_mod_mpoly_one(G, ctx);
    fmpz_mod_mpoly_set(Abar, A, ctx);
    fmpz_mod_mpoly_set(Bbar, B, ctx);
    success = 1;
    goto cleanup;
}
