/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "n_poly.h"
#include "mpn_extras.h"

int n_poly_mod_is_canonical(const n_poly_t A, nmod_t mod)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i] >= mod.n)
            return 0;
        if (A->coeffs[i] == 0 && i + 1 == A->length)
            return 0;
    }
    return 1;
}

void n_poly_mod_set_coeff_ui(
    n_poly_t poly,
    slong j,
    ulong c,
    nmod_t ctx)
{
    if (c >= ctx.n)
        NMOD_RED(c, c, ctx);

    n_poly_set_coeff(poly, j, c);
}

void n_poly_mod_add_ui(n_poly_t res, const n_poly_t poly, ulong c, nmod_t ctx)
{
    if (c >= ctx.n)
        NMOD_RED(c, c, ctx);

    if (poly->length < 1)
    {
        n_poly_set_ui(res, c);
    }
    else
    {
        n_poly_set(res, poly);
        res->coeffs[0] = nmod_add(res->coeffs[0], c, ctx);
        _n_poly_normalise(res);
   }
}

mp_limb_t n_poly_mod_div_root(n_poly_t Q,
                                     const n_poly_t A, mp_limb_t c, nmod_t ctx)
{
    mp_limb_t rem;

    slong len = A->length;

    if (len < 2)
    {
        if (len == 1)
        {
            rem = A->coeffs[0];
            n_poly_zero(Q);
            return rem;
        }

        n_poly_zero(Q);
        return 0;
    }

    n_poly_fit_length(Q, len - 1);
    rem = _nmod_poly_div_root(Q->coeffs, A->coeffs, len, c, ctx);
    Q->length = len - 1;
    return rem;
}

void n_poly_mod_pow(n_poly_t res, const n_poly_t poly, ulong e, nmod_t ctx)
{
    const slong len = poly->length;
    slong rlen;

    if ((len < 2) | (e < UWORD(3)))
    {
        if (len == 0)
        {
            if (e == 0)
                n_poly_one(res);
            else
                n_poly_zero(res);
        }
        else if (len == 1)
        {
            n_poly_set_ui(res,
                     n_powmod2_ui_preinv(poly->coeffs[0], e, ctx.n, ctx.ninv));
        }
        else if (e == 0)
        {
            n_poly_one(res);
        }
        else if (e == 1)
            n_poly_set(res, poly);
        else  /* e == UWORD(2) */
            n_poly_mod_mul(res, poly, poly, ctx);

        return;
    }

    rlen = (slong) e * (len - 1) + 1;

    if (res != poly)
    {
        n_poly_fit_length(res, rlen);
        _nmod_poly_pow(res->coeffs, poly->coeffs, len, e, ctx);
    }
    else
    {
        n_poly_t t;
        n_poly_init2(t, rlen);
        _nmod_poly_pow(t->coeffs, poly->coeffs, len, e, ctx);
        n_poly_swap(res, t);
        n_poly_clear(t);
    }

    res->length = rlen;
    _n_poly_normalise(res);
}


void n_poly_mod_mul(n_poly_t res, const n_poly_t poly1,
                 const n_poly_t poly2, nmod_t ctx)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    if (len1 == 0 || len2 == 0)
    {
        n_poly_zero(res);

        return;
    }

    len_out = poly1->length + poly2->length - 1;

    if (res == poly1 || res == poly2)
    {
        n_poly_t temp;

        n_poly_init2(temp, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, ctx);
        else
            _nmod_poly_mul(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, ctx);

        n_poly_swap(temp, res);
        n_poly_clear(temp);
    }
    else
    {
        n_poly_fit_length(res, len_out);

        if (len1 >= len2)
            _nmod_poly_mul(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, ctx);
        else
            _nmod_poly_mul(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, ctx);
    }

    res->length = len_out;
    _n_poly_normalise(res);
}

void n_poly_mod_mullow(
    n_poly_t res,
    const n_poly_t poly1,
    const n_poly_t poly2,
    slong trunc,
    nmod_t ctx)
{
    slong len1, len2, len_out;

    len1 = poly1->length;
    len2 = poly2->length;

    len_out = poly1->length + poly2->length - 1;
    if (trunc > len_out)
        trunc = len_out;

    if (len1 <= 0 || len2 <= 0 || trunc <= 0)
    {
        n_poly_zero(res);
        return;
    }

    if (res == poly1 || res == poly2)
    {
        n_poly_t temp;

        n_poly_init2(temp, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(temp->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, ctx);
        else
            _nmod_poly_mullow(temp->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, ctx);

        n_poly_swap(temp, res);
        n_poly_clear(temp);
    }
    else
    {
        n_poly_fit_length(res, trunc);

        if (len1 >= len2)
            _nmod_poly_mullow(res->coeffs, poly1->coeffs, len1,
                           poly2->coeffs, len2, trunc, ctx);
        else
            _nmod_poly_mullow(res->coeffs, poly2->coeffs, len2,
                           poly1->coeffs, len1, trunc, ctx);
    }

    res->length = trunc;
    _n_poly_normalise(res);
}


void n_poly_mod_mulmod(n_poly_t res, const n_poly_t poly1,
                            const n_poly_t poly2, const n_poly_t f, nmod_t ctx)
{
    slong len1, len2, lenf;
    mp_ptr fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_mulmod). Divide by zero.\n");
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        n_poly_zero(res);
        return;
    }

    if (len1 + len2 - lenf > 0)
    {
        if (f == res)
        {
            fcoeffs = flint_malloc(sizeof(mp_limb_t) * lenf);
            _nmod_vec_set(fcoeffs, f->coeffs, lenf);
        }
        else
            fcoeffs = f->coeffs;

        n_poly_fit_length(res, lenf - 1);
        _nmod_poly_mulmod(res->coeffs, poly1->coeffs, len1,
                                       poly2->coeffs, len2,
                                       fcoeffs, lenf,
                                       ctx);
        if (f == res)
            flint_free(fcoeffs);

        res->length = lenf - 1;
        _n_poly_normalise(res);
    }
    else
    {
        n_poly_mod_mul(res, poly1, poly2, ctx);
    }
}


void n_poly_mod_div(n_poly_t Q, const n_poly_t A, const n_poly_t B, nmod_t ctx)
{
    n_poly_t tQ;
    mp_ptr q;
    slong A_len, B_len;

    B_len = B->length;

    if (B_len == 0)
    {
        if (ctx.n == 1)
        {
            n_poly_set(Q, A);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (n_poly_mod_div). Division by zero.\n");
        }
    }

    A_len = A->length;

    if (A_len < B_len)
    {
        n_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        n_poly_init2(tQ, A_len - B_len + 1);
        q = tQ->coeffs;
    }
    else
    {
        n_poly_fit_length(Q, A_len - B_len + 1);
        q = Q->coeffs;
    }

    _nmod_poly_div(q, A->coeffs, A_len, B->coeffs, B_len, ctx);

    if (Q == A || Q == B)
    {
        n_poly_swap(tQ, Q);
        n_poly_clear(tQ);
    }

    Q->length = A_len - B_len + 1;
}

void n_poly_mod_rem(n_poly_t R, const n_poly_t A, const n_poly_t B, nmod_t ctx)
{
    const slong lenA = A->length, lenB = B->length;
    n_poly_t tR;
    mp_ptr r;

    if (lenB == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (nmod_poly_rem). Division by zero.\n");
    }
    if (lenA < lenB)
    {
        n_poly_set(R, A);
        return;
    }

    if (R == A || R == B)
    {
        n_poly_init2(tR, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        n_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_rem(r, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (R == A || R == B)
    {
        n_poly_swap(R, tR);
        n_poly_clear(tR);
    }

    R->length = lenB - 1;
    _n_poly_normalise(R);
}


void n_poly_mod_divrem(n_poly_t Q, n_poly_t R,
                                const n_poly_t A, const n_poly_t B, nmod_t ctx)
{
    const slong lenA = A->length, lenB = B->length;
    n_poly_t tQ, tR;
    mp_ptr q, r;

    if (lenB == 0)
    {
        if (ctx.n == 1)
        {
            n_poly_set(Q, A);
            n_poly_zero(R);
            return;
        }
        else
        {
            flint_throw(FLINT_ERROR, "Exception (n_poly_mod_divrem). Division by zero.");
        }
    }

    if (lenA < lenB)
    {
        n_poly_set(R, A);
        n_poly_zero(Q);
        return;
    }

    if (Q == A || Q == B)
    {
        n_poly_init2(tQ, lenA - lenB + 1);
        q = tQ->coeffs;
    }
    else
    {
        n_poly_fit_length(Q, lenA - lenB + 1);
        q = Q->coeffs;
    }

    if (R == A || R == B)
    {
        n_poly_fit_length(tR, lenB - 1);
        r = tR->coeffs;
    }
    else
    {
        n_poly_fit_length(R, lenB - 1);
        r = R->coeffs;
    }

    _nmod_poly_divrem(q, r, A->coeffs, lenA, B->coeffs, lenB, ctx);

    if (Q == A || Q == B)
    {
        n_poly_swap(Q, tQ);
        n_poly_clear(tQ);
    }
    if (R == A || R == B)
    {
        n_poly_swap(R, tR);
        n_poly_clear(tR);
    }

    Q->length = lenA - lenB + 1;
    R->length = lenB - 1;

    _n_poly_normalise(R);
}


int n_poly_mod_invmod(n_poly_t A, const n_poly_t B, const n_poly_t P, nmod_t ctx)
{
    const slong lenB = B->length, lenP = P->length;
    mp_limb_t * a;
    n_poly_t tA;
    int ans;

    if (lenP < 2)
        flint_throw(FLINT_ERROR, "lenP < 2 in %s\n", __func__);

    if (lenB == 0)
    {
        n_poly_zero(A);
        return 0;
    }
    if (lenB >= lenP)
    {
        n_poly_t T;

        n_poly_init(T);
        n_poly_mod_rem(T, B, P, ctx);
        ans = n_poly_mod_invmod(A, T, P, ctx);
        n_poly_clear(T);
        return ans;
    }

    if (A != B && A != P)
    {
        n_poly_fit_length(A, lenP - 1);
        a = A->coeffs;
    }
    else
    {
        n_poly_init2(tA, lenP - 1);
        a = tA->coeffs;
    }

    ans = _nmod_poly_invmod(a, B->coeffs, lenB, P->coeffs, lenP, ctx);

    if (A == B || A == P)
    {
        n_poly_swap(A, tA);
        n_poly_clear(tA);
    }

    A->length = lenP - 1;
    _n_poly_normalise(A);
    return ans;
}


void n_poly_mod_gcd(n_poly_t G, const n_poly_t A, const n_poly_t B, nmod_t ctx)
{
    if (A->length < B->length)
    {
        n_poly_mod_gcd(G, B, A, ctx);
    }
    else /* lenA >= lenB >= 0 */
    {
        slong lenA = A->length, lenB = B->length, lenG;
        n_poly_t tG;
        mp_ptr g;

        if (lenA == 0) /* lenA = lenB = 0 */
        {
            n_poly_zero(G);
        }
        else if (lenB == 0) /* lenA > lenB = 0 */
        {
            n_poly_mod_make_monic(G, A, ctx);
        }
        else /* lenA >= lenB >= 1 */
        {
            if (G == A || G == B)
            {
                n_poly_init2(tG, FLINT_MIN(lenA, lenB));
                g = tG->coeffs;
            }
            else
            {
                n_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }

            lenG = _nmod_poly_gcd(g, A->coeffs, lenA, B->coeffs, lenB, ctx);

            if (G == A || G == B)
            {
                n_poly_swap(tG, G);
                n_poly_clear(tG);
            }
            G->length = lenG;

            if (G->length == 1)
                G->coeffs[0] = 1;
            else
                n_poly_mod_make_monic(G, G, ctx);
        }
    }
}

void n_poly_mod_xgcd(
    n_poly_t G,
    n_poly_t S,
    n_poly_t T,
    const n_poly_t A,
    const n_poly_t B,
    nmod_t ctx)
{
    if (A->length < B->length)
    {
        n_poly_mod_xgcd(G, T, S, B, A, ctx);
    }
    else  /* lenA >= lenB >= 0 */
    {
        const slong lenA = A->length, lenB = B->length;
        mp_limb_t inv;

        if (lenA == 0)  /* lenA = lenB = 0 */
        {
            n_poly_zero(G);
            n_poly_zero(S);
            n_poly_zero(T);
        }
        else if (lenB == 0)  /* lenA > lenB = 0 */
        {
            inv = n_invmod(A->coeffs[lenA - 1], ctx.n);
            _n_poly_mod_scalar_mul_nmod(G, A, inv, ctx);
            n_poly_zero(T);
            n_poly_set_coeff(S, 0, inv);
            S->length = 1;
        }
        else if (lenB == 1)  /* lenA >= lenB = 1 */
        {
            n_poly_fit_length(T, 1);
            T->length = 1;
            T->coeffs[0] = n_invmod(B->coeffs[0], ctx.n);
            n_poly_one(G);
            n_poly_zero(S);
        }
        else  /* lenA >= lenB >= 2 */
        {
            mp_ptr g, s, t;
            slong lenG;

            if (G == A || G == B)
            {
                g = _nmod_vec_init(FLINT_MIN(lenA, lenB));
            }
            else
            {
                n_poly_fit_length(G, FLINT_MIN(lenA, lenB));
                g = G->coeffs;
            }
            if (S == A || S == B)
            {
                s = _nmod_vec_init(lenB - 1);
            }
            else
            {
                n_poly_fit_length(S, lenB - 1);
                s = S->coeffs;
            }
            if (T == A || T == B)
            {
                t = _nmod_vec_init(lenA - 1);
            }
            else
            {
                n_poly_fit_length(T, lenA - 1);
                t = T->coeffs;
            }

            if (lenA >= lenB)
                lenG = _nmod_poly_xgcd(g, s, t, A->coeffs, lenA,
                                                         B->coeffs, lenB, ctx);
            else
                lenG = _nmod_poly_xgcd(g, t, s, B->coeffs, lenB,
                                                         A->coeffs, lenA, ctx);

            if (G == A || G == B)
            {
                flint_free(G->coeffs);
                G->coeffs = g;
                G->alloc  = FLINT_MIN(lenA, lenB);
            }
            if (S == A || S == B)
            {
                flint_free(S->coeffs);
                S->coeffs = s;
                S->alloc  = lenB - 1;
            }
            if (T == A || T == B)
            {
                flint_free(T->coeffs);
                T->coeffs = t;
                T->alloc  = lenA - 1;
            }

            G->length = lenG;
            S->length = FLINT_MAX(lenB - lenG, 1);
            T->length = FLINT_MAX(lenA - lenG, 1);
            MPN_NORM(S->coeffs, S->length);
            MPN_NORM(T->coeffs, T->length);

            if (G->coeffs[lenG - 1] != 1)
            {
                inv = nmod_inv(G->coeffs[lenG - 1], ctx);
                _n_poly_mod_scalar_mul_nmod(G, G, inv, ctx);
                _n_poly_mod_scalar_mul_nmod(S, S, inv, ctx);
                _n_poly_mod_scalar_mul_nmod(T, T, inv, ctx);
            }
        }
    }
}


void n_poly_reverse(n_poly_t output, const n_poly_t input, slong m)
{
    n_poly_fit_length(output, m);
    _nmod_poly_reverse(output->coeffs, input->coeffs, input->length, m);
    output->length = m;
    _n_poly_normalise(output);
}


void n_poly_mod_mulmod_preinv(
    n_poly_t res,
    const n_poly_t poly1,
    const n_poly_t poly2,
    const n_poly_t f,
    const n_poly_t finv,
    nmod_t ctx)
{
    slong len1, len2, lenf;
    mp_ptr fcoeffs;

    lenf = f->length;
    len1 = poly1->length;
    len2 = poly2->length;

    if (lenf <= len1 || lenf <= len2)
    {
        flint_throw(FLINT_ERROR, "n_poly_mod_mulmod_preinv: Input is larger than modulus.");
    }

    if (lenf == 1 || len1 == 0 || len2 == 0)
    {
        n_poly_zero(res);
        return;
    }

    if (len1 + len2 > lenf)
    {
        if (f == res)
        {
            fcoeffs = flint_malloc(sizeof(mp_limb_t) * lenf);
            _nmod_vec_set(fcoeffs, f->coeffs, lenf);
        }
        else
        {
            fcoeffs = f->coeffs;
        }

        n_poly_fit_length(res, lenf - 1);
        _nmod_poly_mulmod_preinv(res->coeffs, poly1->coeffs, len1,
                                       poly2->coeffs, len2,
                                       fcoeffs, lenf,
                                       finv->coeffs, finv->length, ctx);
        if (f == res)
            flint_free(fcoeffs);

        res->length = lenf - 1;
        _n_poly_normalise(res);
    }
    else
    {
        n_poly_mod_mul(res, poly1, poly2, ctx);
    }
}

void n_poly_mod_inv_series(n_poly_t Qinv, const n_poly_t Q, slong n, nmod_t ctx)
{
    slong Qlen = Q->length;

    Qlen = FLINT_MIN(Qlen, n);

    if (Qlen == 0)
    {
        flint_throw(FLINT_ERROR, "n_poly_mod_inv_series_newton: Division by zero.");
    }

    if (Qinv != Q)
    {
        n_poly_fit_length(Qinv, n);
        _nmod_poly_inv_series_newton(Qinv->coeffs, Q->coeffs, Qlen, n, ctx);
    }
    else
    {
        n_poly_t t;
        n_poly_init2(t, n);
        _nmod_poly_inv_series_newton(t->coeffs, Q->coeffs, Qlen, n, ctx);
        n_poly_swap(Qinv, t);
        n_poly_clear(t);
    }

    Qinv->length = n;
    _n_poly_normalise(Qinv);
}


void n_poly_mod_div_series(n_poly_t Q, const n_poly_t A, const n_poly_t B,
                                                       slong order, nmod_t ctx)
{
    slong Blen = B->length;
    slong Alen = A->length;

    if (order < 1 || Blen == 0 || B->coeffs[0] == 0)
    {
        flint_throw(FLINT_ERROR, "Exception (n_poly_div_series). Division by zero.\n");
    }

    if (Alen == 0)
    {
        n_poly_zero(Q);
        return;
    }

    if (Q != A && Q != B)
    {
        n_poly_fit_length(Q, order);
        _nmod_poly_div_series(Q->coeffs, A->coeffs, Alen, B->coeffs, Blen, order, ctx);
    }
    else
    {
        n_poly_t t;
        n_poly_init(t);
        _nmod_poly_div_series(t->coeffs, A->coeffs, Alen, B->coeffs, Blen, order, ctx);
        n_poly_swap(Q, t);
        n_poly_clear(t);
    }

    Q->length = order;
    _n_poly_normalise(Q);
}

void n_poly_mod_scalar_mul_ui(n_poly_t A, const n_poly_t B, mp_limb_t c, nmod_t ctx)
{
    if (c >= ctx.n)
    {
        NMOD_RED(c, c, ctx);
    }

    if (c == 0 || B->length < 1)
    {
        n_poly_zero(A);
        return;
    }

    _n_poly_mod_scalar_mul_nmod(A, B, c, ctx);
    _n_poly_normalise(A);
}

/* multiply A by (x^k + c) */
void n_poly_mod_shift_left_scalar_addmul(n_poly_t A, slong k, mp_limb_t c,
                                                                    nmod_t ctx)
{
    mp_limb_t * Acoeffs;
    slong i;
    slong Alen = A->length;

    n_poly_fit_length(A, Alen + k);

    Acoeffs = A->coeffs;

    flint_mpn_copyd(Acoeffs + k, Acoeffs, Alen);
    flint_mpn_zero(Acoeffs, k);

    for (i = 0; i < A->length; i++)
        Acoeffs[i] = nmod_addmul(Acoeffs[i], c, Acoeffs[i + k], ctx);

    A->length = Alen + k;
}

/* A = B + C*(d1*x+d0) */
void n_poly_mod_addmul_linear(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    mp_limb_t d1, mp_limb_t d0,
    nmod_t ctx)
{
    slong i;
    mp_limb_t * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(B->length, C->length + 1);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_poly_fit_length(A, Alen);
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Ccoeffs = C->coeffs;

    for (i = 0; i < Alen; i++)
    {
        ulong p1, p0, t0 = 0, t1 = 0, t2 = 0;

        if (i < Blen)
        {
            t0 = Bcoeffs[i];
        }
        if (i < Clen)
        {
            umul_ppmm(p1, p0, Ccoeffs[i], d0);
            add_ssaaaa(t1, t0, t1, t0, p1, p0);
        }
        if (0 < i && i - 1 < Clen)
        {
            umul_ppmm(p1, p0, Ccoeffs[i - 1], d1);
            add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        }
        NMOD_RED3(Acoeffs[i], t2, t1, t0, ctx);
    }

    A->length = Alen;
    _n_poly_normalise(A);
}

/* A = B + C*d0 */
void n_poly_mod_scalar_addmul_nmod(
    n_poly_t A,
    const n_poly_t B,
    const n_poly_t C,
    mp_limb_t d0,
    nmod_t ctx)
{
    slong i;
    mp_limb_t t0, t1;
    mp_limb_t * Acoeffs, * Bcoeffs, * Ccoeffs;
    slong Blen = B->length;
    slong Clen = C->length;
    slong Alen = FLINT_MAX(B->length, C->length);

    n_poly_fit_length(A, Alen);
    Acoeffs = A->coeffs;
    Bcoeffs = B->coeffs;
    Ccoeffs = C->coeffs;

    if (ctx.norm >= (FLINT_BITS + 1)/2)
    {
        for (i = 0; i + 2 <= FLINT_MIN(Blen, Clen); i += 2)
        {
            NMOD_RED2(t0, 0, Bcoeffs[i + 0] + d0*Ccoeffs[i + 0], ctx);
            NMOD_RED2(t1, 0, Bcoeffs[i + 1] + d0*Ccoeffs[i + 1], ctx);
            Acoeffs[i + 0] = t0;
            Acoeffs[i + 1] = t1;
        }
        for ( ; i < FLINT_MIN(Blen, Clen); i++)
        {
            NMOD_RED2(Acoeffs[i], 0, Bcoeffs[i] + d0*Ccoeffs[i], ctx);
        }

        for ( ; i + 2 <= Clen; i += 2)
        {
            NMOD_RED2(t0, 0, d0*Ccoeffs[i + 0], ctx);
            NMOD_RED2(t1, 0, d0*Ccoeffs[i + 1], ctx);
            Acoeffs[i + 0] = t0;
            Acoeffs[i + 1] = t1;
        }
        for ( ; i < Clen; i++)
        {
            NMOD_RED2(Acoeffs[i], 0, d0*Ccoeffs[i], ctx);
        }
    }
    else
    {
        for (i = 0; i < FLINT_MIN(Blen, Clen); i++)
        {
            t0 = Bcoeffs[i];
            NMOD_ADDMUL(t0, d0, Ccoeffs[i], ctx);
            Acoeffs[i] = t0;
        }

        while (i < Clen)
        {
            Acoeffs[i] = nmod_mul(d0, Ccoeffs[i], ctx);
            i++;
        }
    }

    while (i < Blen)
    {
        Acoeffs[i] = Bcoeffs[i];
        i++;
    }

    A->length = Alen;
    _n_poly_normalise(A);
}


ulong n_poly_mod_remove(n_poly_t f, const n_poly_t p, nmod_t ctx)
{
    n_poly_t q, r;
    ulong i = 0;

    n_poly_init(q);
    n_poly_init(r);

    while (1)
    {
        if (f->length < p->length)
            break;
        n_poly_mod_divrem(q, r, f, p, ctx);
        if (r->length == 0)
            n_poly_swap(q, f);
        else
            break;
        i++;
    }

    n_poly_clear(q);
    n_poly_clear(r);

    return i;
}

mp_limb_t _n_poly_eval_pow(n_poly_t P, n_poly_t alphapow, int nlimbs, nmod_t ctx)
{
    mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t res;
    slong k;

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
        for (k = oldlength; k < Plen; k++)
            alpha_powers[k] = nmod_mul(alpha_powers[k - 1], alpha_powers[1], ctx);
    }

    NMOD_VEC_DOT(res, k, Plen, Pcoeffs[k], alpha_powers[k], ctx, nlimbs);

    return res;
}

mp_limb_t n_poly_mod_eval_pow(n_poly_t P, n_poly_t alphapow, nmod_t ctx)
{
    int nlimbs = _nmod_vec_dot_bound_limbs(P->length, ctx);
    return _n_poly_eval_pow(P, alphapow, nlimbs, ctx);
}

void n_poly_mod_eval2_pow(
    mp_limb_t * vp,
    mp_limb_t * vm,
    const n_poly_t P,
    n_poly_t alphapow,
    nmod_t ctx)
{
    const mp_limb_t * Pcoeffs = P->coeffs;
    slong Plen = P->length;
    mp_limb_t * alpha_powers = alphapow->coeffs;
    mp_limb_t p1, p0, a0, a1, a2, q1, q0, b0, b1, b2;
    slong k;

    a0 = a1 = a2 = 0;
    b0 = b1 = b2 = 0;

    if (Plen > alphapow->length)
    {
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, Plen);
        for (k = oldlength; k < Plen; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                               alphapow->coeffs[1], ctx);
        }
        alphapow->length = Plen;
        alpha_powers = alphapow->coeffs;
    }

    for (k = 0; k + 2 <= Plen; k += 2)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        umul_ppmm(q1, q0, Pcoeffs[k + 1], alpha_powers[k + 1]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
        add_sssaaaaaa(b2, b1, b0, b2, b1, b0, 0, q1, q0);
    }

    if (k < Plen)
    {
        umul_ppmm(p1, p0, Pcoeffs[k + 0], alpha_powers[k + 0]);
        add_sssaaaaaa(a2, a1, a0, a2, a1, a0, 0, p1, p0);
        k++;
    }

    FLINT_ASSERT(k == Plen);

    NMOD_RED3(p0, a2, a1, a0, ctx);
    NMOD_RED3(q0, b2, b1, b0, ctx);

    *vp = nmod_add(p0, q0, ctx);
    *vm = nmod_sub(p0, q0, ctx);
}

mp_limb_t n_poly_mod_eval_step2(
    n_poly_t Acur,
    const n_poly_t Ainc,
    nmod_t mod)
{
    slong i, Alen = Acur->length;
    mp_limb_t * cur = Acur->coeffs;
    const mp_limb_t * inc = Ainc->coeffs;
    ulong t0, t1, t2, p0, p1;

    FLINT_ASSERT(2*Alen == Ainc->length);

    t2 = t1 = t0 = 0;
    for (i = 0; i < Alen; i++)
    {
        umul_ppmm(p1, p0, cur[i], inc[2*i + 0]);
        add_sssaaaaaa(t2, t1, t0, t2, t1, t0, 0, p1, p0);
        cur[i] = nmod_mul(cur[i], inc[2*i + 1], mod);
    }
    NMOD_RED3(t0, t2, t1, t0, mod);
    return t0;
}

