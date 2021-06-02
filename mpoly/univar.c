/*
    Copyright (C) 2021 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpoly.h"
#include "fmpz_mod_mpoly.h"

static void _fmpz_max(fmpz_t a, const fmpz_t b, const fmpz_t c)
{
    fmpz_set(a, fmpz_cmp(b, c) > 0 ? b : c);
}

#define COEFF(A, i) ((void*)(A->coeffs + (i)*R->elem_size))
#define IS_ZERO(A) (R->is_zero(A, R->ctx))
#define ZERO(A) (R->zero(A, R->ctx))
#define ONE(A) (R->one(A, R->ctx))
#define SET(A, B) (R->set(A, B, R->ctx))
#define SWAP(A, B) (R->swap(A, B, R->ctx))
#define NEG(A, B) (R->neg(A, B, R->ctx))
#define ADD(A, B, C) (R->add(A, B, C, R->ctx))
#define SUB(A, B, C) (R->sub(A, B, C, R->ctx))
#define MUL_FMPZ(A, B, C) (R->mul_fmpz(A, B, C, R->ctx))
#define MUL(A, B, C) (R->mul(A, B, C, R->ctx))
#define DIVEXACT(A, B, C) (R->divexact(A, B, C, R->ctx))
#define DIVIDES(A, B, C) (R->divides(A, B, C, R->ctx))
#define POW_FMPZ(A, B, C) (R->pow_fmpz(A, B, C, R->ctx))

void * mpoly_void_ring_elem_init(mpoly_void_ring_t R)
{
    void * a = flint_malloc(R->elem_size);
    R->init(a, R->ctx);
    return a;
}

void mpoly_void_ring_elem_clear(void * a, mpoly_void_ring_t R)
{
    R->clear(a, R->ctx);
    flint_free(a);
}

void mpoly_univar_init(mpoly_univar_t A, mpoly_void_ring_t R)
{
    A->coeffs = NULL;
    A->exps = NULL;
    A->alloc = 0;
    A->length = 0;
}

void mpoly_univar_clear(mpoly_univar_t A, mpoly_void_ring_t R)
{
    slong i;

    for (i = 0; i < A->alloc; i++)
    {
        R->clear(COEFF(A, i), R->ctx);
        fmpz_clear(A->exps + i);
    }

    flint_free(A->coeffs);
    flint_free(A->exps);
}

void mpoly_univar_swap(mpoly_univar_t A, mpoly_univar_t B)
{
    mpoly_univar_struct t = *B;
    *B = *A;
    *A = t;
}

void mpoly_univar_fit_length(mpoly_univar_t A, slong len, mpoly_void_ring_t R)
{
    slong i;
    slong old_alloc = A->alloc;
    slong new_alloc = FLINT_MAX(len, 2*A->alloc);

    if (len <= old_alloc)
        return;

    A->exps = FLINT_ARRAY_REALLOC(A->exps, new_alloc, fmpz);
    A->coeffs = flint_realloc(A->coeffs, new_alloc*R->elem_size);

    for (i = old_alloc; i < new_alloc; i++)
    {
        fmpz_init(A->exps + i);
        R->init(COEFF(A, i), R->ctx);
    }

    A->alloc = new_alloc;
}

void mpoly_univar_init2(mpoly_univar_t A, slong len, mpoly_void_ring_t R)
{
    mpoly_univar_init(A, R);
    mpoly_univar_fit_length(A, len, R);
}

/*
    A = prem(A, -B)
    C is used for working space
*/
void mpoly_univar_prem(
    mpoly_univar_t A,
    const mpoly_univar_t B,
    mpoly_univar_t C,
    mpoly_void_ring_t R)
{
    slong i, j;
    fmpz_t z1, delta, delta_org;
    void * u = mpoly_void_ring_elem_init(R);
    void * v = mpoly_void_ring_elem_init(R);

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(B != C);
    FLINT_ASSERT(C != A);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(fmpz_cmp(A->exps + 0, B->exps + 0) >= 0);

    fmpz_init(z1);
    fmpz_init(delta);
    fmpz_init(delta_org);

    fmpz_sub(delta_org, A->exps + 0, B->exps + 0);
    fmpz_add_ui(delta_org, delta_org, 1);

looper:

    if (A->length < 1)
        goto done;

    fmpz_sub(delta, A->exps + 0, B->exps + 0);
    if (fmpz_sgn(delta) < 0)
        goto done;

    i = 1;
    j = 1;
    C->length = 0;
    while (i < A->length || j < B->length)
    {
        mpoly_univar_fit_length(C, C->length + 1, R);

        if (j < B->length)
            fmpz_add(z1, B->exps + j, delta);

        if (i < A->length && j < B->length && fmpz_equal(A->exps + i, z1))
        {
            MUL(u, COEFF(A, i), COEFF(B, 0));
            MUL(v, COEFF(A, 0), COEFF(B, j));
            SUB(COEFF(C, C->length), v, u);
            fmpz_set(C->exps + C->length, A->exps + i);
            i++;
            j++;
        }
        else if (i < A->length && (j >= B->length ||
                                                fmpz_cmp(A->exps + i, z1) > 0))
        {
            MUL(COEFF(C, C->length), COEFF(A, i), COEFF(B, 0));
            NEG(COEFF(C, C->length), COEFF(C, C->length));
            fmpz_set(C->exps + C->length, A->exps + i);
            i++;
        }
        else
        {
            FLINT_ASSERT(j < B->length && (i >= A->length ||
                                               fmpz_cmp(A->exps + i, z1) < 0));

            MUL(COEFF(C, C->length), COEFF(A, 0), COEFF(B, j));
            fmpz_set(C->exps + C->length, z1);
            j++;
        }

        C->length += !IS_ZERO(COEFF(C, C->length));
    }

    mpoly_univar_swap(A, C);
    fmpz_sub_ui(delta_org, delta_org, 1);
    goto looper;

done:

    FLINT_ASSERT(fmpz_sgn(delta_org) >= 0);

    if (!fmpz_is_zero(delta_org))
    {
        NEG(v, COEFF(B, 0));
        POW_FMPZ(u, v, delta_org);
        for (i = 0; i < A->length; i++)
            MUL(COEFF(A, i), COEFF(A, i), u);
    }

    mpoly_void_ring_elem_clear(u, R);
    mpoly_void_ring_elem_clear(v, R);

    fmpz_clear(z1);
    fmpz_clear(delta);
    fmpz_clear(delta_org);
}

/*
    G = pseudo-gcd of B and A
      = last nonzero subresultant polynomial starting with B and A

    B and A are clobbered
*/
int mpoly_univar_pseudo_gcd_ducos(
    mpoly_univar_t G,
    mpoly_univar_t B,
    mpoly_univar_t A,
    mpoly_void_ring_t R)
{
    slong i, j, k, aJ, ae;
    fmpz_t n, d, e, J, z1, alpha;
    int iexists, jexists, kexists;
    void * u, * v, * w, * s;
    mpoly_univar_t C, D, H, T;
    mpoly_univar_struct * last;

    FLINT_ASSERT(B->length > 0);
    FLINT_ASSERT(A->length > 0);
    FLINT_ASSERT(fmpz_cmp(B->exps + 0, A->exps + 0) >= 0);
    FLINT_ASSERT(fmpz_sgn(A->exps + 0) >= 0);

    if (fmpz_is_zero(A->exps + 0))
    {
        mpoly_univar_fit_length(G, 1, R);
        G->length = 1;
        fmpz_zero(G->exps + 0);
        return POW_FMPZ(COEFF(G, 0), COEFF(A, 0), B->exps + 0);
    }

    fmpz_init(n);
    fmpz_init(d);
    fmpz_init(e);
    fmpz_init(J);
    fmpz_init(z1);
    fmpz_init(alpha);

    u = mpoly_void_ring_elem_init(R);
    v = mpoly_void_ring_elem_init(R);
    w = mpoly_void_ring_elem_init(R);
    s = mpoly_void_ring_elem_init(R);

    i = FLINT_MAX(B->length, A->length);
    mpoly_univar_init2(C, i + 1, R);
    mpoly_univar_init2(D, i + 1, R);
    mpoly_univar_init2(H, i + 1, R);
    mpoly_univar_init2(T, i + 1, R);

    last = A;

    fmpz_sub(z1, B->exps + 0, A->exps + 0);
    POW_FMPZ(s, A->coeffs + 0, z1);

    mpoly_univar_prem(B, A, D, R);

looper:

    if (B->length < 1)
        goto done;

    fmpz_set(d, A->exps + 0);
    fmpz_set(e, B->exps + 0);

    last = B;

    fmpz_sub(z1, d, e);
    if (fmpz_is_one(z1))
    {
        if (fmpz_is_zero(e))
            goto done;

        /* D = (B[e]*A - A[e]*B)/A[d] */
        /*           i        j       */
        i = 1;
        j = 1;
        if (A->length > 1 && fmpz_equal(A->exps + 1, e))
            i++;
        else
            j = B->length;
        D->length = 0;
        while (i < A->length || j < B->length)
        {
            mpoly_univar_fit_length(D, D->length + 1, R);

            if (i < A->length && j < B->length &&
                                          fmpz_equal(A->exps + i, B->exps + j))
            {
                MUL(u, COEFF(A, i), COEFF(B, 0));
                MUL(v, COEFF(A, 1), COEFF(B, j));
                SUB(w, u, v);
                DIVEXACT(COEFF(D, D->length), w, COEFF(A, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;                
            }
            else if (i < A->length && (j >= B->length ||
                                       fmpz_cmp(A->exps + i, B->exps + j) > 0))
            {
                MUL(u, COEFF(A, i), COEFF(B, 0));
                DIVEXACT(COEFF(D, D->length), u, COEFF(A, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT((j < B->length && (i >= A->length ||
                                     fmpz_cmp(B->exps + j, A->exps + i) > 0)));

                MUL(v, COEFF(A, 1), COEFF(B, j));
                DIVEXACT(COEFF(D, D->length), v, COEFF(A, 0));
                NEG(COEFF(D, D->length), COEFF(D, D->length));
                fmpz_set(D->exps + D->length, B->exps + j);
                j++;
            }

            D->length += !IS_ZERO(COEFF(D, D->length));
        }

        /* A = (B[e]*(D - B*x) + B[e-1]*B)/s */
        /*            i    j            k    */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (B->length > 1 && fmpz_equal(B->exps + 1, z1))
        {
            j = 2;            
            k = 1;
        }
        else
        {
            j = 1;
            k = B->length;
        }

        A->length = 0;
        while (i < D->length || j < B->length || k < B->length)
        {
            fmpz * exp;

            mpoly_univar_fit_length(A, A->length + 1, R);

            exp = A->exps + A->length;

            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < B->length)
            {
                fmpz_add_ui(z1, B->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length) && fmpz_equal(exp, D->exps + i);
            jexists = (j < B->length) && fmpz_equal(exp, z1);
            kexists = (k < B->length) && fmpz_equal(exp, B->exps + k);

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    SUB(w, COEFF(D, i), COEFF(B, j));
                    MUL(u, COEFF(B, 0), w);
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(D, i));
                }

                if (kexists)
                {
                    MUL(v, COEFF(B, 1), COEFF(B, k));
                    ADD(w, u, v);
                    DIVEXACT(COEFF(A, A->length), w, s);
                }
                else
                {
                    DIVEXACT(COEFF(A, A->length), u, s);
                }
            }
            else
            {
                if (kexists)
                {
                    MUL(u, COEFF(B, 1), COEFF(B, k));
                    if (jexists)
                    {
                        MUL(v, COEFF(B, 0), COEFF(B, j));
                        SUB(w, u, v);
                        DIVEXACT(COEFF(A, A->length), w, s);
                    }
                    else
                    {
                        DIVEXACT(COEFF(A, A->length), u, s);
                    }
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(B, j));
                    DIVEXACT(COEFF(A, A->length), u, s);
                    NEG(COEFF(A, A->length), COEFF(A, A->length));
                }
            }

            A->length += !IS_ZERO(COEFF(A, A->length));

            i += iexists;
            j += jexists;
            k += kexists;
        }

        mpoly_univar_swap(A, B);
        SET(s, COEFF(A, 0));
        last = A;
    }
    else
    {
        fmpz_sub(n, d, e);
        fmpz_sub_ui(n, n, 1);
        fmpz_one(alpha);
        while (fmpz_add(z1, alpha, alpha), fmpz_cmp(z1, n) <= 0)
            fmpz_set(alpha, z1);

        SET(u, COEFF(B, 0));
        fmpz_sub(n, n, alpha);
        while (fmpz_cmp_ui(alpha, 1) > 0)
        {
            fmpz_tdiv_q_2exp(alpha, alpha, 1);
            MUL(v, u, u);
            DIVEXACT(u, v, s);
            if (fmpz_cmp(n, alpha) >= 0)
            {
                MUL(v, u, COEFF(B, 0));
                DIVEXACT(u, v, s);
                fmpz_sub(n, n, alpha);
            }
        }

        mpoly_univar_fit_length(C, B->length, R);
        for (i = 0; i < B->length; i++)
        {
            MUL(v, u, COEFF(B, i));
            DIVEXACT(COEFF(C, i), v, s);
            fmpz_set(C->exps + i, B->exps + i);
        }
        C->length = B->length;

        last = C;

        if (fmpz_is_zero(e))
            goto done;

        /* H = C - C[e]*x^e */
        mpoly_univar_fit_length(H, C->length, R);
        for (i = 1; i < C->length; i++)
        {
            SET(COEFF(H, i - 1), COEFF(C, i));
            fmpz_set(H->exps + i - 1, C->exps + i);
        }
        H->length = C->length - 1;

        /* D = C[e]*A - A[e]*H  (truncated to powers of x < e) */
        i = 0;
        j = H->length;
        ae = A->length;
        while (i < A->length && fmpz_cmp(A->exps + i, e) >= 0)
        {
            if (fmpz_equal(A->exps + i, e))
            {
                j = 0;
                ae = i;
            }
            i++;
        }
        D->length = 0;
        while (i < A->length || j < H->length)
        {
            mpoly_univar_fit_length(D, D->length + 1, R);

            if (i < A->length && j < H->length &&
                                          fmpz_equal(A->exps + i, H->exps + j))
            {
                MUL(u, COEFF(A, i), COEFF(C, 0));
                MUL(v, COEFF(A, ae), COEFF(H, j));
                SUB(COEFF(D, D->length), u, v);
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
                j++;                
            }
            else if (i < A->length && (j >= H->length ||
                                       fmpz_cmp(A->exps + i, H->exps + j) > 0))
            {
                MUL(COEFF(D, D->length), COEFF(A, i), COEFF(C, 0));
                fmpz_set(D->exps + D->length, A->exps + i);
                i++;
            }
            else
            {
                FLINT_ASSERT(j < H->length && (i >= A->length ||
                                      fmpz_cmp(H->exps + j, A->exps + i) > 0));

                MUL(COEFF(D, D->length), COEFF(A, ae), COEFF(H, j));
                NEG(COEFF(D, D->length), COEFF(D, D->length));
                fmpz_set(D->exps + D->length, H->exps + j);
                j++;
            }

            D->length += !IS_ZERO(COEFF(D, D->length));
        }

        for (fmpz_add_ui(J, e, 1); fmpz_cmp(J, d) < 0; fmpz_add_ui(J, J, 1))
        {
            if (H->length < 1)
                break;

            /* H = H*x - H[e-1]*B/B[e] */
            fmpz_sub_ui(z1, e, 1);
            if (fmpz_equal(H->exps + 0, z1))
            {
                i = 1;
                j = 1;
                T->length = 0;
                while (i < H->length || j < B->length)
                {
                    mpoly_univar_fit_length(T, T->length + 1, R);

                    if (i < H->length)
                        fmpz_add_ui(z1, H->exps + i, 1);

                    if (i < H->length && j < B->length &&
                                                   fmpz_equal(z1, B->exps + j))
                    {
                        MUL(u, COEFF(H, 0), COEFF(B, j));
                        DIVEXACT(v, u, COEFF(B, 0));
                        SUB(COEFF(T, T->length), COEFF(H, i), v);
                        fmpz_set(T->exps + T->length, B->exps + j);
                        i++;
                        j++;                
                    }
                    else if (i < H->length && (j >= B->length ||
                                                fmpz_cmp(z1, B->exps + j) > 0))
                    {
                        SET(COEFF(T, T->length), COEFF(H, i));
                        fmpz_set(T->exps + T->length, z1);
                        i++;
                    }
                    else
                    {
                        FLINT_ASSERT(j < B->length && (i >= H->length ||
                                               fmpz_cmp(z1, B->exps + j) < 0));

                        MUL(u, COEFF(H, 0), COEFF(B, j));
                        DIVEXACT(COEFF(T, T->length), u, COEFF(B, 0));
                        NEG(COEFF(T, T->length), COEFF(T, T->length));
                        fmpz_set(T->exps + T->length, B->exps + j);
                        j++;                
                    }

                    T->length += !IS_ZERO(COEFF(T, T->length));
                }

                mpoly_univar_swap(H, T);
            }
            else
            {
                FLINT_ASSERT(fmpz_cmp(H->exps + 0, z1) < 0);
                for (i = 0; i < H->length; i++)
                    fmpz_add_ui(H->exps + i, H->exps + i, 1);
            }

            /* find coefficient of x^J in A */
            aJ = 0;
            while (aJ < A->length && !fmpz_equal(A->exps + aJ, J))
                aJ++;
            if (aJ >= A->length)
                continue;

            /* D = D - A[J]*H */
            i = 0;
            j = 0;
            T->length = 0;
            while (i < D->length || j < H->length)
            {
                mpoly_univar_fit_length(T, T->length + 1, R);

                if (i < D->length && j < H->length &&
                                          fmpz_equal(D->exps + i, H->exps + j))
                {
                    MUL(u, COEFF(H, j), COEFF(A, aJ));
                    SUB(COEFF(T, T->length), COEFF(D, i), u);
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                    j++;                
                }
                else if (i < D->length && (j >= H->length ||
                                       fmpz_cmp(D->exps + i, H->exps + j) > 0))
                {
                    SET(COEFF(T, T->length), COEFF(D, i));
                    fmpz_set(T->exps + T->length, D->exps + i);
                    i++;
                }
                else
                {
                    FLINT_ASSERT(j < H->length && (i >= D->length ||
                                      fmpz_cmp(D->exps + i, H->exps + j) < 0));

                    MUL(COEFF(T, T->length), COEFF(H, j), COEFF(A, aJ));
                    NEG(COEFF(T, T->length), COEFF(T, T->length));
                    fmpz_set(T->exps + T->length, H->exps + j);
                    j++;
                }

                T->length += !IS_ZERO(COEFF(T, T->length));
            }
            mpoly_univar_swap(D, T);
        }

        /* B = (-1)^(d-e+1) * (B[e]*(D/A[d] - H*x) +  H[e-1]*B)/s */
        i = 0;
        fmpz_sub_ui(z1, e, 1);
        if (H->length > 0 && fmpz_equal(H->exps + 0, z1))
        {
            j = 1;
            k = 1;
        }
        else
        {
            j = 0;
            k = B->length;
        }
        T->length = 0;
        while (i < D->length || j < H->length || k < B->length)
        {
            fmpz * exp;

            mpoly_univar_fit_length(T, T->length + 1, R);

            exp = T->exps + T->length;
            fmpz_zero(exp);

            if (i < D->length)
                _fmpz_max(exp, exp, D->exps + i);

            if (j < H->length)
            {
                fmpz_add_ui(z1, H->exps + j, 1);
                _fmpz_max(exp, exp, z1);
            }

            if (k < B->length)
                _fmpz_max(exp, exp, B->exps + k);

            iexists = (i < D->length && fmpz_equal(exp, D->exps + i));
            jexists = (j < H->length && fmpz_equal(exp, z1));
            kexists = (k < B->length && fmpz_equal(exp, B->exps + k));

            FLINT_ASSERT(iexists || jexists || kexists);

            if (iexists)
            {
                if (jexists)
                {
                    DIVEXACT(u, COEFF(D, i), COEFF(A, 0));
                    SUB(w, u, COEFF(H, j));
                    MUL(u, COEFF(B, 0), w);
                }
                else
                {
                    DIVEXACT(u, COEFF(D, i), COEFF(A, 0));
                    MUL(u, COEFF(B, 0), u);
                }
                if (kexists)
                {
                    MUL(v, COEFF(H, 0), COEFF(B, k));
                    ADD(w, u, v);
                    DIVEXACT(COEFF(T, T->length), w, s);
                }
                else
                {
                    DIVEXACT(COEFF(T, T->length), u, s);
                }
            }
            else
            {
                if (kexists)
                {
                    MUL(u, COEFF(H, 0), COEFF(B, k));
                    if (jexists)
                    {
                        MUL(v, COEFF(B, 0), COEFF(H, j));
                        SUB(w, u, v);
                        DIVEXACT(COEFF(T, T->length), w, s);
                    }
                    else
                    {
                        DIVEXACT(COEFF(T, T->length), u, s);
                    }
                }
                else
                {
                    MUL(u, COEFF(B, 0), COEFF(H, j));
                    DIVEXACT(COEFF(T, T->length), u, s);
                    NEG(COEFF(T, T->length), COEFF(T, T->length));
                }
            }

            if (((fmpz_get_ui(d) - fmpz_get_ui(e)) & 1) == 0)
                NEG(COEFF(T, T->length), COEFF(T, T->length));            

            T->length += !IS_ZERO(COEFF(T, T->length));

            i += iexists;
            j += jexists;
            k += kexists;
        }

        mpoly_univar_swap(B, T);
        mpoly_univar_swap(A, C);
        SET(s, COEFF(A, 0));
        last = A;
    }

    goto looper;

done:

    mpoly_univar_swap(G, last);

    fmpz_clear(n);
    fmpz_clear(d);
    fmpz_clear(e);
    fmpz_clear(J);
    fmpz_clear(z1);
    fmpz_clear(alpha);
    mpoly_void_ring_elem_clear(u, R);
    mpoly_void_ring_elem_clear(v, R);
    mpoly_void_ring_elem_clear(w, R);
    mpoly_void_ring_elem_clear(s, R);
    mpoly_univar_clear(C, R);
    mpoly_univar_clear(D, R);
    mpoly_univar_clear(H, R);
    mpoly_univar_clear(T, R);
    return 1;
}



void mpoly_univar_derivative(
    mpoly_univar_t A,
    const mpoly_univar_t B,
    mpoly_void_ring_t R)
{
    slong Ai, Bi;

    mpoly_univar_fit_length(A, B->length, R);

    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    {
        if (fmpz_sgn(B->exps + Bi) <= 0)
            continue;

        MUL_FMPZ(COEFF(A, Ai), COEFF(B, Bi), B->exps + Bi);
        fmpz_sub_ui(A->exps + Ai, B->exps + Bi, 1);
        Ai += !IS_ZERO(COEFF(A, Ai));
    }
    A->length = Ai;
}

/* fx and gx are clobbered */
int mpoly_univar_resultant(
    void * r,
    mpoly_univar_t fx,
    mpoly_univar_t gx,
    mpoly_void_ring_t R)
{
    int success, change_sign;
    mpoly_univar_t rx;
    mpoly_univar_struct * F, * G;

    if (fx->length < 1 || gx->length < 1)
    {
        ZERO(r);
        return 1;
    }

    mpoly_univar_init(rx, R);

    if (fmpz_cmp(fx->exps + 0, gx->exps + 0) < 0)
    {
        change_sign = 1 & fmpz_get_ui(fx->exps + 0) & fmpz_get_ui(gx->exps + 0);
        F = gx;
        G = fx;
    }
    else
    {
        change_sign = 0;
        F = fx;
        G = gx;
    }

    if (fmpz_is_zero(G->exps + 0))
    {
        success = POW_FMPZ(r, COEFF(G, 0), F->exps + 0);
    }
    else
    {
        success = mpoly_univar_pseudo_gcd_ducos(rx, F, G, R);

        if (success && rx->length == 1 && fmpz_is_zero(rx->exps + 0))
            SWAP(r, COEFF(rx, 0));
        else
            ZERO(r);
    }

    if (success && change_sign)
        NEG(r, r);

    mpoly_univar_clear(rx, R);

    return success;
}

/* fx is clobbered */
int mpoly_univar_discriminant(
    void * d,
    mpoly_univar_t fx,
    mpoly_void_ring_t R)
{
    int success;
    mpoly_univar_t rx, fxp;

    if (fx->length < 1 || fmpz_cmp_ui(fx->exps + fx->length - 1, 1) > 0)
    {
        /* the discriminant of the zero polynomial should be zero */
        ZERO(d);
        return 1;
    }

    if (fmpz_is_zero(fx->exps + 0))
    {
        /* the discriminant of the constant polynomial a should be 1/a^2 */
        ONE(d);
        success = DIVIDES(d, d, COEFF(fx, 0));
        if (success)
            MUL(d, d, d);
        return success;
    }

    if (fmpz_is_one(fx->exps + 0))
    {
        /* the discriminant of a linear polynomial should be 1 */
        ONE(d);
        return 1;
    }

    /* the discriminant is (-1)^(n*(n-1)/2)*res(f,f')*a_n^(n-m-2), m=deg(f') */
    mpoly_univar_init(rx, R);
    mpoly_univar_init(fxp, R);
    mpoly_univar_derivative(fxp, fx, R);

    if (fxp->length < 1)
    {
        ZERO(d);
        success = 1;
    }
    else
    {
        int change_sign = fmpz_get_ui(fx->exps + 0) & 2;
        fmpz_t exp_diff;
        void * u;

        fmpz_init(exp_diff);
        fmpz_sub(exp_diff, fx->exps + 0, fxp->exps + 0);
        fmpz_sub_ui(exp_diff, exp_diff, 2);

        u = mpoly_void_ring_elem_init(R);
        SET(u, COEFF(fx, 0));   /* will be clobbered */

        success = mpoly_univar_pseudo_gcd_ducos(rx, fx, fxp, R);

        if (success && rx->length == 1 && fmpz_is_zero(rx->exps + 0))
        {
            if (change_sign != 0)
                NEG(COEFF(rx, 0), COEFF(rx, 0));

            if (fmpz_sgn(exp_diff) < 0)
            {
                FLINT_ASSERT(fmpz_equal_si(exp_diff, -1));
                DIVEXACT(d, COEFF(rx, 0), u);
            }
            else
            {
                success = POW_FMPZ(u, u, exp_diff);
                if (success)
                    MUL(d, COEFF(rx, 0), u);
            }
        }
        else
        {
            ZERO(d);
        }

        fmpz_clear(exp_diff);
        mpoly_void_ring_elem_clear(u, R);
    }

    mpoly_univar_clear(rx, R);
    mpoly_univar_clear(fxp, R);

    return success;
}

