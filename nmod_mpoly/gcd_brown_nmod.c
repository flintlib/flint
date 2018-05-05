/*
    Copyright (C) 2018 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <http://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly.h"

/*
    permute the variable ording in the dense represenation
    to something favourable for gcd_brown
*/
int nmod_mpolyd_ctx_set_for_gcd(nmod_mpolyd_ctx_t dctx,
                            const nmod_mpoly_t A, const nmod_mpoly_t B,
                                                    const nmod_mpoly_ctx_t ctx)
{
    int success;
    slong i, j, degb_prod;
    slong * Aexps, * Bexps, * deg_bounds;
    slong nvars = ctx->minfo->nvars;
    slong * perm = dctx->perm;
    TMP_INIT;

    success = 0;

    TMP_START;
    Aexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    Bexps = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    if (A->bits > FLINT_BITS || B->bits > FLINT_BITS)
        goto cleanup;
    nmod_mpoly_degrees_si(Aexps, A, ctx);
    nmod_mpoly_degrees_si(Bexps, B, ctx);

    deg_bounds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (i = 0; i < nvars; i++)
    {
        dctx->perm[i] = i;
    }

    degb_prod = 1;
    for (i = 0; i < nvars; i++)
    {
        ulong hi;
        deg_bounds[i] = FLINT_MAX(Aexps[i] + 1, Bexps[i] + 1);
        umul_ppmm(hi, degb_prod, degb_prod, deg_bounds[i]);
        if (hi != WORD(0) || degb_prod < 0)
            goto cleanup;
    }

    success = 1;
    for (i = 1; i < nvars; i++)
    {
        for (j = i; (j > 0) && (deg_bounds[j-1] < deg_bounds[j-0]); j--)
        {
            slong t1, t2;
            t1 = deg_bounds[j-1];
            t2 = deg_bounds[j-0];
            deg_bounds[j-0] = t1;
            deg_bounds[j-1] = t2;
            t1 = perm[j-1];
            t2 = perm[j-0];
            perm[j-0] = t1;
            perm[j-1] = t2;
        }
    }

cleanup:
    TMP_END;
    return success;
}

/*
    Set G = gcd(A,B), Abar = A/G, Bbar = B/G
*/
void nmod_mpolyd_gcd_brown_univar(nmod_mpolyd_t G,
                                 nmod_mpolyd_t Abar,       nmod_mpolyd_t Bbar,
                           const nmod_mpolyd_t A   , const nmod_mpolyd_t B,
                                                        const nmodf_ctx_t fctx)
{
    slong Alen, Blen;

    Alen = A->deg_bounds[0];
    while (A->coeffs[Alen-1] == 0) {
        Alen --;
        if (Alen == 0)
            break;
    }

    Blen = B->deg_bounds[0];
    while (B->coeffs[Blen-1] == 0) {
        Blen --;
        if (Blen == 0)
            break;
    }

    nmod_mpolyd_set_nvars(G, 1);
    nmod_mpolyd_set_nvars(Abar, 1);
    nmod_mpolyd_set_nvars(Bbar, 1);
    if (Alen == 0)
    {
        if (Blen == 0)
        {
            G->coeffs[0] = 0;
            G->deg_bounds[0] = 1;
            Abar->coeffs[0] = 0;
            Abar->deg_bounds[0] = 1;
            Bbar->coeffs[0] = 0;
            Bbar->deg_bounds[0] = 1;
        } else
        {
            nmod_mpolyd_fit_length(G, Blen);
            _nmod_poly_make_monic(G->coeffs, B->coeffs, Blen, fctx->mod);
            G->deg_bounds[0] = Blen;
            Abar->coeffs[0] = 0;
            Abar->deg_bounds[0] = 1;
            Bbar->coeffs[0] = 1;
            Bbar->deg_bounds[0] = 1;
        }

    } else {
        if (Blen == 0) {
            nmod_mpolyd_fit_length(G, Alen);
            _nmod_poly_make_monic(G->coeffs, A->coeffs, Alen, fctx->mod);
            G->deg_bounds[0] = Alen;
            Abar->coeffs[0] = 1;
            Abar->deg_bounds[0] = 1;
            Bbar->coeffs[0] = 0;
            Bbar->deg_bounds[0] = 1;
        } else {
            if (Alen >= Blen) {
                nmod_mpolyd_fit_length(G, Blen);
                G->deg_bounds[0] = _nmod_poly_gcd_euclidean(G->coeffs,
                                         A->coeffs, Alen, B->coeffs, Blen,
                                                                fctx->mod);
            } else {
                nmod_mpolyd_fit_length(G, Alen);
                G->deg_bounds[0] = _nmod_poly_gcd_euclidean(G->coeffs,
                                         B->coeffs, Blen, A->coeffs, Alen,
                                                                fctx->mod);
            }
            if (G->deg_bounds[0] <= 1)
            {
                G->coeffs[0] = 1;
            } else
            {
                _nmod_poly_make_monic(G->coeffs,
                                   G->coeffs, G->deg_bounds[0], fctx->mod);
            }

            Abar->deg_bounds[0] = Alen - G->deg_bounds[0] + 1;
            nmod_mpolyd_fit_length(Abar, Alen - G->deg_bounds[0] + 1);
            _nmod_poly_div(Abar->coeffs, A->coeffs, Alen,
                                   G->coeffs, G->deg_bounds[0], fctx->mod);

            Bbar->deg_bounds[0] = Blen - G->deg_bounds[0] + 1;
            nmod_mpolyd_fit_length(Bbar, Blen - G->deg_bounds[0] + 1);
            _nmod_poly_div(Bbar->coeffs, B->coeffs, Blen,
                                   G->coeffs, G->deg_bounds[0], fctx->mod);
        }
    }
    return;
}

/*
    Set cont to the content of A in the last (innermost) variable only
*/
void nmod_mpolyd_last_content(nmod_poly_t cont, const nmod_mpolyd_t A,
                                                        const nmodf_ctx_t fctx) 
{
    mp_limb_t * temp;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
    slong i, j, Plen;
    slong degb_prod, degb_last=0;
    TMP_INIT;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    TMP_START;
    nmod_poly_zero(cont);
    nmod_poly_fit_length(cont, degb_last);
    temp = (mp_limb_t *) TMP_ALLOC(degb_last*sizeof(mp_limb_t));

    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        if (Plen != 0)
        {
            if (cont->length == 0)
            {
                _nmod_poly_make_monic(cont->coeffs, P, Plen, fctx->mod);
                cont->length = Plen;
            } else
            {
                if (cont->length < Plen) {
                    cont->length = _nmod_poly_gcd_euclidean(temp,
                                         P, Plen, cont->coeffs, cont->length,
                                                                    fctx->mod);
                } else {
                    cont->length = _nmod_poly_gcd_euclidean(temp,
                                         cont->coeffs, cont->length, P, Plen,
                                                                    fctx->mod);
                }
                if (cont->length <= 1)
                {
                    cont->coeffs[0] = 1;
                } else
                {
                    _nmod_poly_make_monic(cont->coeffs, temp, cont->length, fctx->mod);
                }
            }
        }
    }
    TMP_END;
}

/*
    Divide A by b in the last (innermost) variable only
*/
void nmod_mpolyd_div_last_poly(nmod_mpolyd_t A, nmod_poly_t b,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_last=0, new_degb_last;
    mp_limb_t * temp;
    TMP_INIT;

    FLINT_ASSERT(b->length != 0);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }
    new_degb_last = degb_last - (b->length - 1);

    TMP_START;
    temp = (mp_limb_t *) TMP_ALLOC(new_degb_last*sizeof(mp_limb_t));

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        slong Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            FLINT_ASSERT(Plen > b->length - 1);
            _nmod_poly_div(temp, P, Plen, b->coeffs, b->length, fctx->mod);
            while (j < Plen - (b->length - 1))
            {
                A->coeffs[k++] = temp[j];
                j++;
            }
        }
        while (j < new_degb_last)
        {
            A->coeffs[k++] = WORD(0);
            j++;
        }
    }

    A->deg_bounds[A->nvars - 1] = new_degb_last;

    TMP_END;
}

/*
    Set M = A * b in the last (innermost) variable only
*/

void nmod_mpolyd_mul_last_poly(nmod_mpolyd_t M,
                   nmod_mpolyd_t A, nmod_poly_t b, const nmodf_ctx_t fctx)
{
    slong i, j, k;
    slong degb_prod, degb_last, new_degb_last;

    FLINT_ASSERT(M != A);
    FLINT_ASSERT(b->length != 0);

    nmod_mpolyd_set_nvars(M, A->nvars);

    new_degb_last = nmod_mpolyd_last_degree(A, fctx) + 1;
    new_degb_last += (b->length - 1);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars - 1; j++) {
        M->deg_bounds[j] = A->deg_bounds[j];
        degb_prod *= A->deg_bounds[j];
    }
    j = A->nvars - 1;
        M->deg_bounds[j] = new_degb_last;
        nmod_mpolyd_fit_length(M, degb_prod*new_degb_last);

        degb_prod *= A->deg_bounds[j];
        degb_last = A->deg_bounds[j];

    k = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t * P = A->coeffs + i;
        slong Plen = degb_last;
        while (P[Plen-1] == 0)
        {
            Plen --;
            if (Plen == 0)
                break;
        }

        j = 0;
        if (Plen != 0)
        {
            _nmod_poly_mul_classical(M->coeffs + k, P, Plen,
                                              b->coeffs, b->length, fctx->mod);
            j +=  Plen + (b->length - 1);
        }
        FLINT_ASSERT(j <= new_degb_last);
        while (j < new_degb_last)
        {

            M->coeffs[k+j] = WORD(0);
            j++;
        }

        k += new_degb_last;
    }
}


void nmod_mpolyd_mul_scalar(nmod_mpolyd_t A, mp_limb_t b, const nmodf_ctx_t fctx)
{
    slong j, degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_prod *= A->deg_bounds[j];
    }

    _nmod_vec_scalar_mul_nmod(A->coeffs, A->coeffs, degb_prod, b, fctx->mod);
}

/*
    Set lc to the leading coefficient of A when considered as a polynomial
    in all but the last variable
*/
void nmod_mpolyd_last_lc(nmod_poly_t lc, const nmod_mpolyd_t A,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, lc_len;
    slong degb_prod, degb_last=0;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++) {
        degb_last = A->deg_bounds[j];
        degb_prod *= degb_last;
    }

    nmod_poly_zero(lc);

    for (i = degb_prod - degb_last; i >= 0; i -= degb_last)
    {
        lc_len = degb_last;
        while (lc_len > 0)
        {
            if (A->coeffs[i+lc_len-1] != 0)
            {
                nmod_poly_fit_length(lc, lc_len);
                flint_mpn_copyi(lc->coeffs, A->coeffs + i, lc_len);
                lc->length = lc_len;
                return;
            }
            lc_len--;
        }
    }

    return;
}

/*
    set v = P(alpha) given precomputed powers of alpha
*/
void _nmod_mpolyd_eval_uni(mp_limb_t * v, mp_limb_t * Pcoeffs, slong Plen,
                              mp_limb_t * alpha_powers, const nmodf_ctx_t fctx)
{
    mp_limb_t pp1, pp0, ac0, ac1, ac2;
    slong k;
    FLINT_ASSERT(Plen >= WORD(1));

    ac0 = Pcoeffs[0];
    ac1 = 0;
    ac2 = 0;

    for (k = WORD(1); k < Plen; k++)
    {
        umul_ppmm(pp1, pp0, Pcoeffs[k], alpha_powers[k]);
        add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
    }
    NMOD_RED3(v[0], ac2, ac1, ac0, fctx->mod);
}

/*
    set v{p|m} = P(+-alpha) given precomputed powers of alpha
*/
void _nmod_mpolyd_eval2_uni(mp_limb_t * vp, mp_limb_t * vm,
                                      mp_limb_t * Pcoeffs, slong Plen,
                              mp_limb_t * alpha_powers, const nmodf_ctx_t fctx)
{
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;
    FLINT_ASSERT(Plen >= WORD(1));

    a0 = a1 = a2 = UWORD(0);
    b0 = b1 = b2 = UWORD(0);

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, WORD(0), q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, WORD(0), p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, fctx->mod);
    NMOD_RED3(q0, b2, b1, b0, fctx->mod);

    vp[0] = nmod_add(p0, q0, fctx->mod);
    vm[0] = nmod_sub(p0, q0, fctx->mod);
}

/*
    Set E to A evaluated at last_var = alpha. E will have one fewer variable.
*/
void nmod_mpolyd_eval_last(nmod_mpolyd_t E, const nmod_mpolyd_t A,
                              mp_limb_t * alpha_powers, const nmodf_ctx_t fctx)
{

    slong i, j, k;
    slong degb_prod, degb_last=0, degb_prod_last=0;

    nmod_mpolyd_set_nvars(E, A->nvars - 1);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod_last = degb_prod;
        degb_prod *= degb_last;
        if (j < E->nvars)
        {
            E->deg_bounds[j] = A->deg_bounds[j];
        }
    }

    nmod_mpolyd_fit_length(E, degb_prod_last);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        mp_limb_t pp0, pp1, ac0=A->coeffs[i + 0], ac1=0, ac2=0, v = WORD(0);
        for (k = 1; k < degb_last; k++)
        {
            umul_ppmm(pp1, pp0, A->coeffs[i + k], alpha_powers[k]);
            add_sssaaaaaa(ac2, ac1, ac0, ac2, ac1, ac0, WORD(0), pp1, pp0);
        }
        NMOD_RED3(v, ac2, ac1, ac0, fctx->mod);
        E->coeffs[j++] = v;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}

/*
    Set E{p|m} to A evaluated at last_var = +-alpha.
*/
void nmod_mpolyd_eval2_last(nmod_mpolyd_t Ep, nmod_mpolyd_t Em,
                           const nmod_mpolyd_t A,  mp_limb_t * alpha_powers,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j;
    slong degb_prod, degb_last=0, degb_prod_last=0;

    nmod_mpolyd_set_nvars(Ep, A->nvars - 1);
    nmod_mpolyd_set_nvars(Em, A->nvars - 1);

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_last = A->deg_bounds[j];
        degb_prod_last = degb_prod;
        degb_prod *= degb_last;
        if (j < Ep->nvars)
        {
            Ep->deg_bounds[j] = A->deg_bounds[j];
            Em->deg_bounds[j] = A->deg_bounds[j];
        }
    }

    nmod_mpolyd_fit_length(Ep, degb_prod_last);
    nmod_mpolyd_fit_length(Em, degb_prod_last);

    j = 0;
    for (i = 0; i < degb_prod; i += degb_last)
    {
        _nmod_mpolyd_eval2_uni(Ep->coeffs + j, Em->coeffs + j,
                                 A->coeffs + i, degb_last, alpha_powers, fctx);
        j++;
    }

    FLINT_ASSERT(j == degb_prod_last);

    return;
}


slong nmod_mpolyd_leadmon(slong * exps, const nmod_mpolyd_t A)
{

    slong i, j, k;
    slong degb_prod;

    degb_prod = WORD(1);
    for (j = 0; j < A->nvars; j++)
    {
        degb_prod *= A->deg_bounds[j];
    }

    for (i = degb_prod-1; i >= 0; i--)
    {
        if (A->coeffs[i] != 0)
            break;
    }

    FLINT_ASSERT(i>=0);

    k = i;
    for (j = A->nvars - 1; j >= 0; j--) 
    {
        ulong m = A->deg_bounds[j];
        ulong e = k % m;
        k = k / m;
        exps[j] = e;
    }
    FLINT_ASSERT(k == 0);

    return i;
}


void nmod_mpolyd_set_ui(nmod_mpolyd_t A, mp_limb_t v)
{
    slong j;

    for (j = 0; j < A->nvars; j++)
        A->deg_bounds[j] = WORD(1);

    A->coeffs[0] = v;
}

void nmod_mpolyd_set(nmod_mpolyd_t A, const nmod_mpolyd_t B)
{
    slong j;
    slong degb_prod;

    nmod_mpolyd_set_nvars(A, B->nvars);

    degb_prod = WORD(1);
    for (j = 0; j < B->nvars; j++)
    {
        A->deg_bounds[j] = B->deg_bounds[j];
        degb_prod *= B->deg_bounds[j];
    }

    nmod_mpolyd_fit_length(A, degb_prod);
    flint_mpn_copyi(A->coeffs, B->coeffs, degb_prod);
}



/*
    F = F + modulus*(N - F(alpha))
*/
void nmod_mpolyd_addinterp(nmod_mpolyd_t F, nmod_mpolyd_t T,
             nmod_mpolyd_t N, nmod_poly_t modulus, mp_limb_t * alpha_powers,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, k, degb_prod;
    slong nvars = F->nvars;
    mp_limb_t v, u, Nvalue, changed = 0;
    int carry, Fok, Nok;
    slong * inds, Find, Nind, Tind;
    TMP_INIT;

    FLINT_ASSERT(N->nvars == nvars - 1);
    FLINT_ASSERT(modulus->length > 0);

    nmod_mpolyd_set_nvars(T, nvars);

    degb_prod = 1;
    for (j = 0; j < nvars - 1; j++)
    {
        T->deg_bounds[j] = FLINT_MAX(F->deg_bounds[j], N->deg_bounds[j]);
        degb_prod *= T->deg_bounds[j];
    }
        T->deg_bounds[j] = FLINT_MAX(1 + nmod_mpolyd_last_degree(F, fctx),
                                                  modulus->length);
    
    nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1]);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars - 1; j++)
        inds[j] = 0;

    Fok = 1;
    Nok = 1;
    Find = 0;
    Nind = 0;
    Tind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        for (k = 0; k < T->deg_bounds[nvars-1]; k++) {
            T->coeffs[Tind+k] = UWORD(0);
        }

        Nvalue = 0;
        if (Nok)
        {
            Nvalue = N->coeffs[Nind];
        }

        if (Fok)
        {
            for (k = 0; k < T->deg_bounds[nvars-1]
                      && k < F->deg_bounds[nvars-1]; k++)
            {
                T->coeffs[Tind+k] = F->coeffs[Find+k];
            }

            _nmod_mpolyd_eval_uni(&v, F->coeffs + Find, F->deg_bounds[nvars-1],
                                                           alpha_powers, fctx);

            v = nmod_sub(Nvalue, v, fctx->mod);

            if (v != UWORD(0)) 
            {
                changed = 1;
                for (k = 0; k < modulus->length; k++)
                {
                    u = nmod_mul(modulus->coeffs[k], v, fctx->mod);
                    T->coeffs[Tind + k] = nmod_add(T->coeffs[Tind + k], u,
                                                                    fctx->mod);
                }
            }

        } else
        {
            if (Nvalue != UWORD(0))
            {
                changed = 1;
                for (k = 0; k < modulus->length; k++)
                {
                    T->coeffs[Tind + k] = nmod_mul(modulus->coeffs[k],  Nvalue,
                                                                    fctx->mod);
                }
            }
        }

        /* move indices to next chunk */

        carry = 1;
        for (j = nvars - 2; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < T->deg_bounds[j])
            {
                carry = 0;
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }

        Tind += T->deg_bounds[nvars-1];

        Find = 0;
        Nind = 0;
        Fok = 1;
        Nok = 1;
        for (j = 0; j < nvars - 1; j++)
        {
            Fok = Fok && (inds[j] < F->deg_bounds[j]);
            Nok = Nok && (inds[j] < N->deg_bounds[j]);
            Find = inds[j] + F->deg_bounds[j]*Find;
            Nind = inds[j] + N->deg_bounds[j]*Nind;
        }
        Find *= F->deg_bounds[nvars - 1];
    }

    if (changed) 
    {
        nmod_mpolyd_swap(F, T);
    }
    TMP_END;
}


int nmod_mpolyd_addinterp2(nmod_mpolyd_t F, nmod_mpolyd_t T,
             nmod_mpolyd_t P, nmod_mpolyd_t M, nmod_poly_t modulus,
                                 mp_limb_t alpha, mp_limb_t * alpha_powers,
                                                        const nmodf_ctx_t fctx)
{
    slong i, j, k, degb_prod;
    slong nvars = F->nvars;
    mp_limb_t Pvalue, Mvalue, Fvalue, Fvaluem, v, u, changed = 0;
    int carry, Fok, Pok, Mok;
    slong * inds, Find, Pind, Mind, Tind;
    TMP_INIT;

    FLINT_ASSERT((modulus->length & WORD(1)) == WORD(1));
    FLINT_ASSERT(P->nvars == nvars - 1);
    FLINT_ASSERT(M->nvars == nvars - 1);
    FLINT_ASSERT(modulus->length > 0);

    nmod_mpolyd_set_nvars(T, nvars);

    degb_prod = 1;
    for (j = 0; j < nvars - 1; j++)
    {
        T->deg_bounds[j] = FLINT_MAX(F->deg_bounds[j],
                                FLINT_MAX(P->deg_bounds[j], M->deg_bounds[j]));
        degb_prod *= T->deg_bounds[j];
    }
        T->deg_bounds[j] = FLINT_MAX(1 + nmod_mpolyd_last_degree(F, fctx),
                                                          modulus->length + 1);
    
    nmod_mpolyd_fit_length(T, degb_prod*T->deg_bounds[nvars-1]);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars - 1; j++)
        inds[j] = 0;

    Fok = 1;
    Pok = 1;
    Mok = 1;
    Find = 0;
    Pind = 0;
    Mind = 0;
    Tind = 0;
    for (i = 0; i < degb_prod; i++)
    {
        for (k = 0; k < T->deg_bounds[nvars-1]; k++) {
            T->coeffs[Tind+k] = UWORD(0);
        }

        Pvalue = 0;
        if (Pok)
        {
            Pvalue = P->coeffs[Pind];
        }

        Mvalue = 0;
        if (Mok)
        {
            Mvalue = M->coeffs[Mind];
        }

        /* f(x) */
        Fvalue = 0;
        Fvaluem = 0;
        if (Fok)
        {
            for (k = 0; k < T->deg_bounds[nvars-1]
                       && k < F->deg_bounds[nvars-1]; k++)
            {
                T->coeffs[Tind+k] = F->coeffs[Find+k];
            }
            _nmod_mpolyd_eval2_uni(&Fvalue, &Fvaluem,
                 F->coeffs + Find, F->deg_bounds[nvars-1], alpha_powers, fctx);
        }

        Fvalue = nmod_sub(Fvalue, Pvalue, fctx->mod);
        Fvaluem = nmod_sub(Fvaluem, Mvalue, fctx->mod);
        
        changed |= Fvalue | Fvaluem;

        u = nmod_sub(Fvalue, Fvaluem, fctx->mod);
        v = nmod_mul(alpha, nmod_add(Fvalue, Fvaluem, fctx->mod), fctx->mod);

        /* + m(x) * ([F(a)-f(a)] - [F(-a)-f(-a)]) * x */
        if (v != UWORD(0))
        {
            for (k = 0; k < modulus->length; k += 2)
            {
                if (k > 0)
                    FLINT_ASSERT(modulus->coeffs[k-1] == 0);
                T->coeffs[Tind + k] = nmod_sub(T->coeffs[Tind + k],
                                     nmod_mul(modulus->coeffs[k], v, fctx->mod),
                                      fctx->mod);
            }
        }

        /* + m(x) * ([F(a)-f(a)] + [F(-a)-f(-a)]) * a */
        if (u != UWORD(0))
        {
            for (k = 0; k < modulus->length; k += 2)
            {
                if (k > 0)
                    FLINT_ASSERT(modulus->coeffs[k-1] == 0);
                T->coeffs[Tind + k + 1] = nmod_sub(T->coeffs[Tind + k + 1],
                                     nmod_mul(modulus->coeffs[k], u, fctx->mod),
                                          fctx->mod);
            }
        }

        /* move indices to next chunk */

        carry = 1;
        for (j = nvars - 2; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < T->deg_bounds[j])
            {
                carry = 0;
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }

        Tind += T->deg_bounds[nvars-1];

        Find = 0;
        Pind = 0;
        Mind = 0;
        Fok = 1;
        Pok = 1;
        Mok = 1;
        for (j = 0; j < nvars - 1; j++)
        {
            Fok = Fok && (inds[j] < F->deg_bounds[j]);
            Pok = Pok && (inds[j] < P->deg_bounds[j]);
            Mok = Mok && (inds[j] < M->deg_bounds[j]);
            Find = inds[j] + F->deg_bounds[j]*Find;
            Pind = inds[j] + P->deg_bounds[j]*Pind;
            Mind = inds[j] + M->deg_bounds[j]*Mind;
        }
        Find *= F->deg_bounds[nvars - 1];
    }

    if (changed != WORD(0)) 
    {
        nmod_mpolyd_swap(F, T);
    }
    TMP_END;

    return changed == WORD(0);
}


void nmod_mpolyd_startinterp(nmod_mpolyd_t fxn,
             nmod_mpolyd_t newvalue)
{
    slong nvars = newvalue->nvars;

    nmod_mpolyd_set_nvars(newvalue, nvars+1);
    newvalue->deg_bounds[nvars] = 1;

    nmod_mpolyd_set(fxn, newvalue);
}

void nmod_mpolyd_startinterp2(nmod_mpolyd_t F, const nmod_mpolyd_t P,
                                 const nmod_mpolyd_t M, mp_limb_t alpha,
                                                        const nmodf_ctx_t fctx)
{
    int Pok, Mok, carry;
    slong Pind, Mind;
    ulong c0, c1, d0, d1;
    slong i, j;
    slong * inds;
    slong nvars = P->nvars;
    slong degb_prod;
    TMP_INIT;

    FLINT_ASSERT(F != P);
    FLINT_ASSERT(F != M);
    FLINT_ASSERT(P->nvars == M->nvars);

    nmod_mpolyd_set_nvars(F, nvars + 1);

    degb_prod = 1;
    for (j = 0; j < nvars; j++)
    {
        F->deg_bounds[j] = FLINT_MAX(P->deg_bounds[j], M->deg_bounds[j]);
        degb_prod *= F->deg_bounds[j];
    }
        F->deg_bounds[j] = 2;
        degb_prod *= F->deg_bounds[j];


    nmod_mpolyd_fit_length(F, degb_prod);

    TMP_START;
    inds = (slong *) TMP_ALLOC(nvars*sizeof(slong));
    for (j = 0; j < nvars; j++)
        inds[j] = 0;
    Pok = 1;
    Mok = 1;
    Pind = 0;
    Mind = 0;

    d0 = n_invmod(UWORD(2), fctx->mod.n);
    d1 = n_invmod(nmod_add(alpha, alpha, fctx->mod), fctx->mod.n);

    for (i = 0; i < degb_prod; i+= 2)
    {
        c0 = c1 = WORD(0);
        carry = 1;

        if (Pok) {
            c0 += P->coeffs[Pind];
            c1 += P->coeffs[Pind];
            Pind++;
        }
        if (Mok)
        {
            c0 = nmod_add(c0, M->coeffs[Mind], fctx->mod);
            c1 = nmod_sub(c1, M->coeffs[Mind], fctx->mod);
            Mind++;
        }

        F->coeffs[i + 0] = nmod_mul(c0, d0, fctx->mod);
        F->coeffs[i + 1] = nmod_mul(c1, d1, fctx->mod);

        Pok = 1;
        Mok = 1;
        for (j = nvars - 1; j >= 0; j--)
        {
            inds[j] += carry;
            if (inds[j] < F->deg_bounds[j])
            {
                carry = 0;
                Pok = Pok && (inds[j] < P->deg_bounds[j]);
                Mok = Mok && (inds[j] < M->deg_bounds[j]);
            } else
            {
                carry = 1;
                inds[j] = 0;
            }
        }
    }

    TMP_END;
}


/*
    Check divisibility of A by D and expection on
    the degree and leading coeff of D
*/
int nmod_mpolyd_divides_univar(nmod_mpolyd_t Q, nmod_mpolyd_t A, nmod_mpolyd_t D,
                 slong expected_deg, mp_limb_t expected_lc, const nmodf_ctx_t fctx)
{
    int success = 0;
    slong i;
    nmod_poly_t q, r;
    nmod_poly_t a, d;

    nmod_poly_init(q, fctx->mod.n);
    nmod_poly_init(r, fctx->mod.n);

    FLINT_ASSERT(A->nvars == WORD(1));

    a->coeffs = A->coeffs;
    a->length = A->deg_bounds[0];
    a->alloc = A->deg_bounds[0];
    a->mod.n = fctx->mod.n;
    a->mod.ninv = fctx->mod.ninv;
    a->mod.norm = fctx->mod.norm;
    while ((a->length > 0) && (a->coeffs[a->length - 1] == 0))
        a->length--;


    d->coeffs = D->coeffs;
    d->length = D->deg_bounds[0];
    d->alloc = D->deg_bounds[0];
    d->mod.n = fctx->mod.n;
    d->mod.ninv = fctx->mod.ninv;
    d->mod.norm = fctx->mod.norm;
    while ((d->length > 0) && (d->coeffs[d->length - 1] == 0))
        d->length--;


    if (d->length != expected_deg + 1
            || d->coeffs[d->length - 1] != expected_lc)
    {
        goto clean_up;
    }

    nmod_poly_divrem_basecase(q, r, a, d);

    if (r->length != 0)
    {
        goto clean_up;
    }

    nmod_mpolyd_set_nvars(Q, 1);
    nmod_mpolyd_fit_length(Q, q->length);
    Q->deg_bounds[0] = q->length;
    for (i = 0; i < q->length; i++)
    {
        Q->coeffs[i] = q->coeffs[i];
    }

    success = 1;

clean_up:

    nmod_poly_clear(q);
    nmod_poly_clear(r);

    return success;
}


int nmod_mpolyd_gcd_brown_smprime_bivar(nmod_mpolyd_t G,
             nmod_mpolyd_t Abar, nmod_mpolyd_t Bbar,
                   nmod_mpolyd_t A, nmod_mpolyd_t B, const nmodf_ctx_t fctx)
{
    int success;
    slong j, bound;
    slong nvars = A->nvars;
    mp_limb_t alpha, gamma_eval, gammam_eval;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus, modulus2;
    nmod_mpolyd_t T, Gs, Abars, Bbars, phiA, phiB, phiAm, phiBm,
                                         gs, abars, bbars, gsm, abarsm, bbarsm;
    slong leadmon_gs_idx, leadmon_gsm_idx;
    slong * leadmon_gs, * leadmon_gsm, * leadmon_Gs;
    slong deggamma, degGs, degA, degB, degAbars, degBbars;
    slong ABlenmax;
    int gstab, astab, bstab, use_stab;
    mp_limb_t * alpha_powers;
    mp_limb_t * alpham_powers;

    FLINT_ASSERT(G != A);
    FLINT_ASSERT(G != B);
    FLINT_ASSERT(A->nvars == B->nvars);
    FLINT_ASSERT(A->nvars == WORD(2));

    leadmon_gs = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_gsm = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_Gs = (slong *) flint_malloc(nvars*sizeof(slong));

    ABlenmax = 1 + FLINT_MAX(A->deg_bounds[nvars-1],B->deg_bounds[nvars-1]);
    alpha_powers = (mp_limb_t *) flint_malloc(ABlenmax*sizeof(mp_limb_t));
    alpham_powers = (mp_limb_t *) flint_malloc(ABlenmax*sizeof(mp_limb_t));

    nmod_poly_init(cA, fctx->mod.n);
    nmod_poly_init(cB, fctx->mod.n);
    nmod_mpolyd_last_content(cA, A, fctx);
    nmod_mpolyd_last_content(cB, B, fctx);

    nmod_mpolyd_div_last_poly(A, cA, fctx);
    nmod_mpolyd_div_last_poly(B, cB, fctx);

    nmod_poly_init(cG, fctx->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, fctx->mod.n);
    nmod_poly_init(cBbar, fctx->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(lcA, fctx->mod.n);
    nmod_poly_init(lcB, fctx->mod.n);
    nmod_mpolyd_last_lc(lcA, A, fctx);
    nmod_mpolyd_last_lc(lcB, B, fctx);

    nmod_poly_init(gamma, fctx->mod.n);
    nmod_poly_gcd_euclidean(gamma, lcA, lcB);

    bound = 1 + nmod_poly_degree(gamma)
               + FLINT_MAX(nmod_mpolyd_last_degree(A, fctx),
                           nmod_mpolyd_last_degree(B, fctx));

    nmod_mpolyd_init(T, nvars);

    nmod_mpolyd_init(Gs, nvars);
    nmod_mpolyd_init(Abars, nvars);
    nmod_mpolyd_init(Bbars, nvars);

    nmod_mpolyd_init(phiA, nvars - 1);
    nmod_mpolyd_init(phiB, nvars - 1);
    nmod_mpolyd_init(phiAm, nvars - 1);
    nmod_mpolyd_init(phiBm, nvars - 1);

    nmod_mpolyd_init(gs, nvars - 1);
    nmod_mpolyd_init(abars, nvars - 1);
    nmod_mpolyd_init(bbars, nvars - 1);
    nmod_mpolyd_init(gsm, nvars - 1);
    nmod_mpolyd_init(abarsm, nvars - 1);
    nmod_mpolyd_init(bbarsm, nvars - 1);

    nmod_poly_init(modulus, fctx->mod.n);
    nmod_poly_init(modulus2, fctx->mod.n);

    if ((fctx->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    use_stab = 1;
    gstab = bstab = astab = 0;

    nmod_poly_one(modulus);

    for (alpha = (fctx->mod.n - UWORD(1))/UWORD(2); alpha != UWORD(0); alpha--)
    {
        mp_limb_t alpha2 = nmod_mul(alpha, alpha, fctx->mod);
        mp_limb_t temp;

        alpha_powers[0] = UWORD(1);
        alpham_powers[0] = UWORD(1);
        for (j = 1; j < ABlenmax; j++)
        {
            alpha_powers[j] = nmod_mul(alpha_powers[j - 1], alpha, fctx->mod);
            alpham_powers[j] = (j & 1) ? fctx->mod.n - alpha_powers[j]
                                       : alpha_powers[j];
        }

        gamma_eval = nmod_poly_evaluate_nmod(gamma, alpha);
        gammam_eval = nmod_poly_evaluate_nmod(gamma, fctx->mod.n - alpha);
        if (gamma_eval == WORD(0) || gammam_eval == WORD(0))
            goto break_continue;

        nmod_mpolyd_eval2_last(phiA, phiAm, A, alpha_powers, fctx);
        nmod_mpolyd_eval2_last(phiB, phiBm, B, alpha_powers, fctx);

        if (use_stab && gstab)
        {
            nmod_mpolyd_eval2_last(gs, gsm, Gs, alpha_powers, fctx);
            success = 1;
            success = success && nmod_mpolyd_divides_univar(abars,  phiA,  gs,
                                              leadmon_gs[0], gamma_eval, fctx);
            success = success && nmod_mpolyd_divides_univar(abarsm, phiAm, gsm,
                                             leadmon_gs[0], gammam_eval, fctx);
            success = success && nmod_mpolyd_divides_univar(bbars,  phiB,  gs,
                                              leadmon_gs[0], gamma_eval, fctx);
            success = success && nmod_mpolyd_divides_univar(bbarsm, phiBm, gsm,
                                             leadmon_gs[0], gammam_eval, fctx);
            if (!success) {
                /*
                    the true gcd does not match the interpolated one
                    this can happen mod 41 where

                    true   gcd = 38*x1^0*x0^0 + 33*x1^1*x0^3 + 6*x1^10*x0^4 + ...
                    interp gcd = 38*x1^0*x0^0 + 33*x1^1*x0^3 + 13*x1^0*x0^4 + ...

                    6*x1^10 == 13 for the values x1 = +-12, +-13
                */
                use_stab = 0;
                nmod_poly_one(modulus);
                alpha = (fctx->mod.n - UWORD(1))/UWORD(2);
                goto break_continue;
            }   

            nmod_mpolyd_mul_scalar(abars, gamma_eval, fctx);
            nmod_mpolyd_mul_scalar(abarsm, gammam_eval, fctx);
            nmod_mpolyd_mul_scalar(bbars, gamma_eval, fctx);
            nmod_mpolyd_mul_scalar(bbarsm, gammam_eval, fctx);

        } else
        {
            success = 1;
            success = success && nmod_mpolyd_gcd_brown_smprime(gs, abars, bbars,
                                                             phiA, phiB, fctx);
            success = success && nmod_mpolyd_gcd_brown_smprime(gsm, abarsm, bbarsm,
                                                           phiAm, phiBm, fctx);
            if (success == 0)
                goto break_continue;

            gstab = astab = bstab = 0;
        }

        leadmon_gs_idx = nmod_mpolyd_leadmon(leadmon_gs, gs);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            goto successful;        
        }

        leadmon_gsm_idx = nmod_mpolyd_leadmon(leadmon_gsm, gs);

        if (leadmon_gsm_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gsm_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            goto successful;        
        }


        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_mpolyd_leadmon(leadmon_Gs, Gs);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] != leadmon_gsm[j])
                    goto break_continue;

                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;
                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    nmod_mpolyd_zero(Gs);
                    nmod_mpolyd_zero(Abars);
                    nmod_mpolyd_zero(Bbars);
                    nmod_poly_one(modulus);
                }
            }
        }

        nmod_mpolyd_mul_scalar(gs, gamma_eval, fctx);
        nmod_mpolyd_mul_scalar(gsm, gammam_eval, fctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(modulus, alpha);
            FLINT_ASSERT(temp ==
                        nmod_poly_evaluate_nmod(modulus, fctx->mod.n - alpha));
            temp = nmod_mul(temp, alpha, fctx->mod);
            temp = nmod_add(temp, temp, fctx->mod);
            nmod_poly_scalar_mul_nmod(modulus, modulus,
                                                  n_invmod(temp, fctx->mod.n));

            if (!gstab)
            {
                gstab = nmod_mpolyd_addinterp2(Gs, T, gs, gsm, modulus, alpha,
                                                           alpha_powers, fctx);
            }
            nmod_mpolyd_addinterp2(Abars, T, abars, abarsm, modulus, alpha,
                                                           alpha_powers, fctx);
            nmod_mpolyd_addinterp2(Bbars, T, bbars, bbarsm, modulus, alpha,
                                                           alpha_powers, fctx);            
            nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha2);
            nmod_poly_shift_left(modulus, modulus, 2);
            nmod_poly_sub(modulus, modulus, modulus2);

        } else
        {
            nmod_poly_one(modulus);
            nmod_mpolyd_startinterp2(Gs, gs, gsm, alpha, fctx);
            nmod_mpolyd_startinterp2(Abars, abars, abarsm, alpha, fctx);
            nmod_mpolyd_startinterp2(Bbars, bbars, bbarsm, alpha, fctx);
            nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha2);
            nmod_poly_shift_left(modulus, modulus, 2);
            nmod_poly_sub(modulus, modulus, modulus2);

            gstab = astab = bstab = 0;
        }

        if (nmod_poly_degree(modulus) < bound)
            continue;

        deggamma = nmod_poly_degree(gamma);
        degGs = nmod_mpolyd_last_degree(Gs, fctx);
        degA = nmod_mpolyd_last_degree(A, fctx);
        degB = nmod_mpolyd_last_degree(B, fctx);
        degAbars = nmod_mpolyd_last_degree(Abars, fctx);
        degBbars = nmod_mpolyd_last_degree(Bbars, fctx);

        if (   deggamma + degA == degGs + degAbars
            && deggamma + degB == degGs + degBbars
           )
        {
            goto successful;
        } else
        {
            nmod_mpolyd_zero(Gs);
            nmod_mpolyd_zero(Abars);
            nmod_mpolyd_zero(Bbars);
            nmod_poly_one(modulus);
            continue;
        }

break_continue:
        (void)(NULL);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_gsm);
    flint_free(leadmon_Gs);
    flint_free(alpha_powers);
    flint_free(alpham_powers);

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);

    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(lcA);
    nmod_poly_clear(lcB);

    nmod_poly_clear(gamma);

    nmod_mpolyd_clear(T);
    nmod_mpolyd_clear(Gs);
    nmod_mpolyd_clear(Abars);
    nmod_mpolyd_clear(Bbars);

    nmod_mpolyd_clear(phiA);
    nmod_mpolyd_clear(phiB);
    nmod_mpolyd_clear(phiAm);
    nmod_mpolyd_clear(phiBm);

    nmod_mpolyd_clear(gs);
    nmod_mpolyd_clear(abars);
    nmod_mpolyd_clear(bbars);
    nmod_mpolyd_clear(gsm);
    nmod_mpolyd_clear(abarsm);
    nmod_mpolyd_clear(bbarsm);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;

successful:

    nmod_poly_init(cGs, fctx->mod.n);
    nmod_poly_init(cAbars, fctx->mod.n);
    nmod_poly_init(cBbars, fctx->mod.n);
    nmod_mpolyd_last_content(cGs, Gs, fctx);
    nmod_mpolyd_last_content(cAbars, Abars, fctx);
    nmod_mpolyd_last_content(cBbars, Bbars, fctx);

    nmod_mpolyd_div_last_poly(Gs, cGs, fctx);
    nmod_mpolyd_div_last_poly(Abars, cAbars, fctx);
    nmod_mpolyd_div_last_poly(Bbars, cBbars, fctx);

    nmod_mpolyd_mul_last_poly(G, Gs, cG, fctx);
    nmod_mpolyd_mul_last_poly(Abar, Abars, cAbar, fctx);
    nmod_mpolyd_mul_last_poly(Bbar, Bbars, cBbar, fctx);

    nmod_poly_clear(cGs);
    nmod_poly_clear(cAbars);
    nmod_poly_clear(cBbars);

    success = 1;
    goto cleanup;
}

int nmod_mpolyd_gcd_brown_smprime(nmod_mpolyd_t G,
                nmod_mpolyd_t Abar, nmod_mpolyd_t Bbar,
                     nmod_mpolyd_t A, nmod_mpolyd_t B,  const nmodf_ctx_t fctx)
{
    int success;
    slong j, bound;
    slong nvars = A->nvars;
    mp_limb_t alpha, gamma_eval, gammam_eval;
    nmod_poly_t cA, cB, cG, cAbar, cBbar, lcA, lcB, gamma;
    nmod_poly_t cGs, cAbars, cBbars, modulus, modulus2;
    nmod_mpolyd_t T, Gs, Abars, Bbars, phiA, phiB, phiAm, phiBm,
                                         gs, abars, bbars, gsm, abarsm, bbarsm;
    slong leadmon_gs_idx, leadmon_gsm_idx;
    slong * leadmon_gs, * leadmon_gsm, * leadmon_Gs;
    slong deggamma, degGs, degA, degB, degAbars, degBbars;
    slong ABlenmax;
    mp_limb_t * alpha_powers;
    mp_limb_t * alpham_powers;

    FLINT_ASSERT(G != A);
    FLINT_ASSERT(G != B);
    FLINT_ASSERT(A->nvars == B->nvars);

    if (A->nvars == 1) {
        nmod_mpolyd_gcd_brown_univar(G, Abar, Bbar, A, B, fctx);
        return 1;
    }
    if (A->nvars == 2) {
        return nmod_mpolyd_gcd_brown_smprime_bivar(G, Abar, Bbar, A, B, fctx);
    }

    leadmon_gs = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_gsm = (slong *) flint_malloc(nvars*sizeof(slong));
    leadmon_Gs = (slong *) flint_malloc(nvars*sizeof(slong));

    ABlenmax = 1 + FLINT_MAX(A->deg_bounds[nvars - 1], B->deg_bounds[nvars - 1]);
    alpha_powers = (mp_limb_t *) flint_malloc(ABlenmax*sizeof(mp_limb_t));
    alpham_powers = (mp_limb_t *) flint_malloc(ABlenmax*sizeof(mp_limb_t));

    nmod_poly_init(cA, fctx->mod.n);
    nmod_poly_init(cB, fctx->mod.n);
    nmod_mpolyd_last_content(cA, A, fctx);
    nmod_mpolyd_last_content(cB, B, fctx);

    nmod_mpolyd_div_last_poly(A, cA, fctx);
    nmod_mpolyd_div_last_poly(B, cB, fctx);

    nmod_poly_init(cG, fctx->mod.n);
    nmod_poly_gcd_euclidean(cG, cA, cB);

    nmod_poly_init(cAbar, fctx->mod.n);
    nmod_poly_init(cBbar, fctx->mod.n);
    nmod_poly_div(cAbar, cA, cG);
    nmod_poly_div(cBbar, cB, cG);

    nmod_poly_init(lcA, fctx->mod.n);
    nmod_poly_init(lcB, fctx->mod.n);
    nmod_mpolyd_last_lc(lcA, A, fctx);
    nmod_mpolyd_last_lc(lcB, B, fctx);

    nmod_poly_init(gamma, fctx->mod.n);
    nmod_poly_gcd_euclidean(gamma, lcA, lcB);

    bound = 1 + nmod_poly_degree(gamma)
               + FLINT_MAX(nmod_mpolyd_last_degree(A, fctx),
                           nmod_mpolyd_last_degree(B, fctx));

    nmod_mpolyd_init(T, nvars);

    nmod_mpolyd_init(Gs, nvars);
    nmod_mpolyd_init(Abars, nvars);
    nmod_mpolyd_init(Bbars, nvars);

    nmod_mpolyd_init(phiA, nvars - 1);
    nmod_mpolyd_init(phiB, nvars - 1);
    nmod_mpolyd_init(phiAm, nvars - 1);
    nmod_mpolyd_init(phiBm, nvars - 1);

    nmod_mpolyd_init(gs, nvars - 1);
    nmod_mpolyd_init(abars, nvars - 1);
    nmod_mpolyd_init(bbars, nvars - 1);
    nmod_mpolyd_init(gsm, nvars - 1);
    nmod_mpolyd_init(abarsm, nvars - 1);
    nmod_mpolyd_init(bbarsm, nvars - 1);

    nmod_poly_init(modulus, fctx->mod.n);
    nmod_poly_init(modulus2, fctx->mod.n);

    nmod_poly_one(modulus);

    if ((fctx->mod.n & UWORD(1)) == UWORD(0))
    {
        success = 0;
        goto cleanup;
    }

    for (alpha = (fctx->mod.n - UWORD(1))/UWORD(2); alpha != UWORD(0); alpha--)
    {
        mp_limb_t alpha2 = nmod_mul(alpha, alpha, fctx->mod);
        mp_limb_t temp;

        alpha_powers[0] = UWORD(1);
        alpham_powers[0] = UWORD(1);
        for (j = 1; j < ABlenmax; j++)
        {
            alpha_powers[j] = nmod_mul(alpha_powers[j - 1], alpha, fctx->mod);
            alpham_powers[j] = (j & 1) ? fctx->mod.n - alpha_powers[j]
                                       : alpha_powers[j];
        }

        gamma_eval = nmod_poly_evaluate_nmod(gamma, alpha);
        gammam_eval = nmod_poly_evaluate_nmod(gamma, fctx->mod.n - alpha);
        if (gamma_eval == WORD(0) || gammam_eval == WORD(0))
            goto break_continue;

        nmod_mpolyd_eval2_last(phiA, phiAm, A, alpha_powers, fctx);
        nmod_mpolyd_eval2_last(phiB, phiBm, B, alpha_powers, fctx);

        success = 1;
        success = success && nmod_mpolyd_gcd_brown_smprime(gs, abars, bbars,
                                                             phiA, phiB, fctx);
        success = success && nmod_mpolyd_gcd_brown_smprime(gsm, abarsm, bbarsm,
                                                           phiAm, phiBm, fctx);
        if (success == 0)
            goto break_continue;

        leadmon_gs_idx = nmod_mpolyd_leadmon(leadmon_gs, gs);

        if (leadmon_gs_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gs_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            goto successful;        
        }

        leadmon_gsm_idx = nmod_mpolyd_leadmon(leadmon_gsm, gs);

        if (leadmon_gsm_idx <= 0)
        {
            FLINT_ASSERT(leadmon_gsm_idx == 0);
            nmod_mpolyd_set_ui(Gs, WORD(1));
            nmod_mpolyd_set(Abars, A);
            nmod_mpolyd_set(Bbars, B);
            goto successful;        
        }

        if (nmod_poly_degree(modulus) > 0)
        {
            nmod_mpolyd_leadmon(leadmon_Gs, Gs);
            for (j = 0; j < nvars - 1; j++)
            {
                if (leadmon_gs[j] != leadmon_gsm[j])
                    goto break_continue;

                if (leadmon_gs[j] > leadmon_Gs[j])
                {
                    goto break_continue;
                } else if (leadmon_gs[j] < leadmon_Gs[j])
                {
                    nmod_mpolyd_zero(Gs);
                    nmod_mpolyd_zero(Abars);
                    nmod_mpolyd_zero(Bbars);
                    nmod_poly_one(modulus);
                }
            }
        }

        nmod_mpolyd_mul_scalar(gs, gamma_eval, fctx);
        nmod_mpolyd_mul_scalar(gsm, gammam_eval, fctx);

        if (nmod_poly_degree(modulus) > 0)
        {
            temp = nmod_poly_evaluate_nmod(modulus, alpha);
            FLINT_ASSERT(temp == 
                        nmod_poly_evaluate_nmod(modulus, fctx->mod.n - alpha));
            temp = nmod_mul(temp, alpha, fctx->mod);
            temp = nmod_add(temp, temp, fctx->mod);
            nmod_poly_scalar_mul_nmod(modulus, modulus,
                                                  n_invmod(temp, fctx->mod.n));

            nmod_mpolyd_addinterp2(Gs, T, gs, gsm, modulus, alpha,
                                                           alpha_powers, fctx);
            nmod_mpolyd_addinterp2(Abars, T, abars, abarsm, modulus, alpha,
                                                           alpha_powers, fctx);
            nmod_mpolyd_addinterp2(Bbars, T, bbars, bbarsm, modulus, alpha,
                                                           alpha_powers, fctx);
        } else
        {
            nmod_poly_one(modulus);
            nmod_mpolyd_startinterp2(Gs, gs, gsm, alpha, fctx);
            nmod_mpolyd_startinterp2(Abars, abars, abarsm, alpha, fctx);
            nmod_mpolyd_startinterp2(Bbars, bbars, bbarsm, alpha, fctx);
        }

        nmod_poly_scalar_mul_nmod(modulus2, modulus, alpha2);
        nmod_poly_shift_left(modulus, modulus, 2);
        nmod_poly_sub(modulus, modulus, modulus2);

        if (nmod_poly_degree(modulus) < bound)
            continue;

        deggamma = nmod_poly_degree(gamma);
        degGs = nmod_mpolyd_last_degree(Gs, fctx);
        degA = nmod_mpolyd_last_degree(A, fctx);
        degB = nmod_mpolyd_last_degree(B, fctx);
        degAbars = nmod_mpolyd_last_degree(Abars, fctx);
        degBbars = nmod_mpolyd_last_degree(Bbars, fctx);
        if (   deggamma + degA == degGs + degAbars
            && deggamma + degB == degGs + degBbars
           )
        {
            goto successful;
        } else
        {
            nmod_mpolyd_zero(Gs);
            nmod_mpolyd_zero(Abars);
            nmod_mpolyd_zero(Bbars);
            nmod_poly_one(modulus);
            continue;
        }

break_continue:
        (void)(NULL);
    }

    success = 0;

cleanup:

    flint_free(leadmon_gs);
    flint_free(leadmon_gsm);
    flint_free(leadmon_Gs);
    flint_free(alpha_powers);
    flint_free(alpham_powers);

    nmod_poly_clear(cA);
    nmod_poly_clear(cB);
    nmod_poly_clear(cG);

    nmod_poly_clear(cAbar);
    nmod_poly_clear(cBbar);

    nmod_poly_clear(lcA);
    nmod_poly_clear(lcB);

    nmod_poly_clear(gamma);

    nmod_mpolyd_clear(T);
    nmod_mpolyd_clear(Gs);
    nmod_mpolyd_clear(Abars);
    nmod_mpolyd_clear(Bbars);

    nmod_mpolyd_clear(phiA);
    nmod_mpolyd_clear(phiB);
    nmod_mpolyd_clear(phiAm);
    nmod_mpolyd_clear(phiBm);

    nmod_mpolyd_clear(gs);
    nmod_mpolyd_clear(abars);
    nmod_mpolyd_clear(bbars);
    nmod_mpolyd_clear(gsm);
    nmod_mpolyd_clear(abarsm);
    nmod_mpolyd_clear(bbarsm);

    nmod_poly_clear(modulus);
    nmod_poly_clear(modulus2);

    return success;

successful:

    nmod_poly_init(cGs, fctx->mod.n);
    nmod_poly_init(cAbars, fctx->mod.n);
    nmod_poly_init(cBbars, fctx->mod.n);
    nmod_mpolyd_last_content(cGs, Gs, fctx);
    nmod_mpolyd_last_content(cAbars, Abars, fctx);
    nmod_mpolyd_last_content(cBbars, Bbars, fctx);

    nmod_mpolyd_div_last_poly(Gs, cGs, fctx);
    nmod_mpolyd_div_last_poly(Abars, cAbars, fctx);
    nmod_mpolyd_div_last_poly(Bbars, cBbars, fctx);

    nmod_mpolyd_mul_last_poly(G, Gs, cG, fctx);
    nmod_mpolyd_mul_last_poly(Abar, Abars, cAbar, fctx);
    nmod_mpolyd_mul_last_poly(Bbar, Bbars, cBbar, fctx);

    nmod_poly_clear(cGs);
    nmod_poly_clear(cAbars);
    nmod_poly_clear(cBbars);

    success = 1;
    goto cleanup;
}
