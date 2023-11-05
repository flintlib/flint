/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod_mpoly_factor.h"

#ifdef FLINT_WANT_ASSERT

static void fmpz_mod_polyu_get_fmpz_mod_polyun(
    fmpz_mod_polyu_t A,
    const fmpz_mod_polyun_t B,
    const fmpz_mod_ctx_t ctx)
{
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->coeffs[i].length - 1; j >= 0; j--)
        {
            if (fmpz_is_zero(B->coeffs[i].coeffs + j))
                continue;
            fmpz_mod_polyu_fit_length(A, A->length + 1, ctx);
            fmpz_set(A->coeffs + A->length, B->coeffs[i].coeffs + j);
            A->exps[A->length] = B->exps[i] + j;
            A->length++;
        }
    }
}


static void fmpz_mod_polyu_sort_terms(fmpz_mod_polyu_t A)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->exps[j - 1] < A->exps[j]; j--)
    {
        fmpz_swap(A->coeffs + j - 1, A->coeffs + j);
        FLINT_SWAP(ulong, A->exps[j - 1], A->exps[j]);
    }
    return;
}


int fmpz_mod_polyu_is_canonical(
    const fmpz_mod_polyu_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_mod_is_canonical(A->coeffs + i, ctx) ||
            fmpz_is_zero(A->coeffs + i))
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}

static void fmpz_mod_polyu_combine_like_terms(
    fmpz_mod_polyu_t A,
    const fmpz_mod_ctx_t ctx)
{
    slong in, out;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && A->exps[out] == A->exps[in])
        {
            fmpz_mod_add(A->coeffs + out, A->coeffs + out, A->coeffs + in, ctx);
        }
        else
        {
            if (out < 0 || !fmpz_is_zero(A->coeffs + out))
                out++;

            if (out != in)
            {
                A->exps[out] = A->exps[in];
                fmpz_set(A->coeffs + out, A->coeffs + in);
            }
        }
    }

    if (out < 0 || !fmpz_is_zero(A->coeffs + out))
        out++;

    A->length = out;
}


static int fmpz_mod_polyu_equal(
    fmpz_mod_polyu_t A,
    fmpz_mod_polyu_t B)
{
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (!fmpz_equal(A->coeffs + i, B->coeffs + i))
            return 0;
        if (A->exps[i] != B->exps[i])
            return 0;
    }
    return 1;
}


static void fmpz_mod_polyu_mul(
    fmpz_mod_polyu_t A,
    const fmpz_mod_polyu_t B,
    const fmpz_mod_polyu_t C,
    const fmpz_mod_ctx_t ctx)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    fmpz_mod_polyu_fit_length(A, B->length*C->length, ctx);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->exps[Ai]   = B->exps[Bi] + C->exps[Ci];
        fmpz_mod_mul(A->coeffs + Ai, B->coeffs + Bi, C->coeffs + Ci, ctx);
        Ai++;
    }
    A->length = Ai;
    fmpz_mod_polyu_sort_terms(A);
    fmpz_mod_polyu_combine_like_terms(A, ctx);
}

#endif


void fmpz_mod_poly_fill_powers(
    fmpz_mod_poly_t alphapow,
    slong target,
    const fmpz_mod_ctx_t ctx)
{
    if (target + 1 > alphapow->length)
    {
        slong k;
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        fmpz_mod_poly_fit_length(alphapow, target + 1, ctx);
        for (k = oldlength; k <= target; k++)
        {
            fmpz_mod_mul(alphapow->coeffs + k, alphapow->coeffs + k - 1,
                                                alphapow->coeffs + 1, ctx);
        }
        alphapow->length = target + 1;
    }
}


void fmpz_mod_polyu3_interp_reduce_2sm_bpoly(
    fmpz_mod_bpoly_t Ap,
    fmpz_mod_bpoly_t Am,
    const fmpz_mod_polyu_t A,
    fmpz_mod_poly_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    slong i;
    slong cur0, cur1, e0, e1, e2;
    fmpz_t t, tp, tm;
    const fmpz * Acoeffs = A->coeffs;
    const ulong * Aexps = A->exps;

    fmpz_init(t);
    fmpz_init(tp);
    fmpz_init(tm);

    fmpz_mod_bpoly_zero(Ap, ctx);
    fmpz_mod_bpoly_zero(Am, ctx);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(Aexps[i], 2, 3);
    cur1 = extract_exp(Aexps[i], 1, 3);
    e2   = extract_exp(Aexps[i], 0, 3);

    fmpz_mod_poly_fill_powers(alpha, e2, ctx);

    fmpz_zero(tp);
    fmpz_zero(tm);
    if (e2 & 1)
        fmpz_mod_mul(tm, alpha->coeffs + e2, Acoeffs + i, ctx);
    else
        fmpz_mod_mul(tp, alpha->coeffs + e2, Acoeffs + i, ctx);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(Aexps[i], 2, 3);
        e1 = extract_exp(Aexps[i], 1, 3);
        e2 = extract_exp(Aexps[i], 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            fmpz_mod_add(t, tp, tm, ctx);
            fmpz_mod_bpoly_set_coeff(Ap, cur0, cur1, t, ctx);
            fmpz_mod_sub(t, tp, tm, ctx);
            fmpz_mod_bpoly_set_coeff(Am, cur0, cur1, t, ctx);
            fmpz_zero(tp);
            fmpz_zero(tm);
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        fmpz_mod_poly_fill_powers(alpha, e2, ctx);

    #if WANT_ASSERT
        fmpz_mod_pow_ui(t, alpha->coeffs + 1, e2, ctx);
        FLINT_ASSERT(fmpz_equal(alpha->coeffs + e2, t));
    #endif

        if (e2 & 1)
            fmpz_mod_addmul(tm, tm, alpha->coeffs + e2, Acoeffs + i, ctx);
        else
            fmpz_mod_addmul(tp, tp, alpha->coeffs + e2, Acoeffs + i, ctx);
    }

    fmpz_mod_add(t, tp, tm, ctx);
    fmpz_mod_bpoly_set_coeff(Ap, cur0, cur1, t, ctx);
    fmpz_mod_sub(t, tp, tm, ctx);
    fmpz_mod_bpoly_set_coeff(Am, cur0, cur1, t, ctx);

    fmpz_clear(t);
    fmpz_clear(tp);
    fmpz_clear(tm);
}


/*
    T(x0, x1, x2) is in F[x2][x0, x1]
    A(x0, x1) are B(x0, x1) are in F[x0, x1]
    set T so that
        T(x0, x1, x2) == A(x0, x1) mod (x2 - alpha)
        T(x0, x1, x2) == B(x0, x1) mod (x2 + alpha)
*/

void fmpz_mod_polyu3n_interp_lift_2sm_bpoly(
    slong * lastdeg,
    fmpz_mod_polyun_t T,
    const fmpz_mod_bpoly_t A,
    const fmpz_mod_bpoly_t B,
    const fmpz_t alpha,
    const fmpz_mod_ctx_t ctx)
{
    slong lastlength = 0;
    fmpz_mod_poly_struct * Tcoeffs;
    ulong * Texps;
    slong Ti;
    fmpz_mod_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    fmpz_mod_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    fmpz_t u, v, Avalue, Bvalue;
    fmpz_t d0, d1;

    fmpz_init(d0);
    fmpz_init(d1);
    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(Avalue);
    fmpz_init(Bvalue);

    fmpz_cdiv_q_2exp(d0, fmpz_mod_ctx_modulus(ctx), 1);
    fmpz_mod_add(d1, alpha, alpha, ctx);
    fmpz_mod_inv(d1, d1, ctx);

    fmpz_mod_polyun_fit_length(T, FLINT_MAX(A->length, B->length), ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;

    Ti = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : fmpz_mod_poly_degree(A->coeffs + Ai, ctx);
    bi = (Bi < 0) ? 0 : fmpz_mod_poly_degree(B->coeffs + Bi, ctx);

    while (Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            fmpz_mod_polyun_fit_length(T, Ti + FLINT_MAX(Ai, Bi) + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        FLINT_ASSERT(Ai < 0 || !fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
        FLINT_ASSERT(Bi < 0 || !fmpz_is_zero(Bcoeffs[Bi].coeffs + bi));

        fmpz_zero(Avalue);
        if (Ai >= 0)
        {
            fmpz_set(Avalue, Acoeffs[Ai].coeffs + ai);
            Texps[Ti] = pack_exp3(Ai, ai, 0);
        }

        fmpz_zero(Bvalue);
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            if (fmpz_is_zero(Avalue))
            {
                fmpz_set(Bvalue, Bcoeffs[Bi].coeffs + bi);
                Texps[Ti] = Bexp;
            }
            else
            {
                if (Texps[Ti] <= Bexp)
                {
                    fmpz_set(Bvalue, Bcoeffs[Bi].coeffs + bi);
                }
                if (Texps[Ti] < Bexp)
                {
                    fmpz_zero(Avalue);
                    Texps[Ti] = Bexp;
                }
            }
        }

        fmpz_mod_sub(u, Avalue, Bvalue, ctx);
        fmpz_mod_add(v, Avalue, Bvalue, ctx);
        fmpz_mod_mul(u, u, d1, ctx);
        fmpz_mod_mul(v, v, d0, ctx);

        FLINT_ASSERT(!fmpz_is_zero(u) || !fmpz_is_zero(v));

        fmpz_mod_poly_fit_length(Tcoeffs + Ti, 2, ctx);
        fmpz_set(Tcoeffs[Ti].coeffs + 0, v);
        fmpz_set(Tcoeffs[Ti].coeffs + 1, u);
        Tcoeffs[Ti].length = 1 + !fmpz_is_zero(u);
        lastlength = FLINT_MAX(lastlength, Tcoeffs[Ti].length);
        Ti++;

        if (!fmpz_is_zero(Avalue))
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx);
            }
        }
        if (!fmpz_is_zero(Bvalue))
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero(Bcoeffs[Bi].coeffs + bi));
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = fmpz_mod_poly_degree(Bcoeffs + Bi, ctx);
            }
        }
    }
    T->length = Ti;

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(T, ctx));

    fmpz_clear(d0);
    fmpz_clear(d1);
    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(Avalue);
    fmpz_clear(Bvalue);

    *lastdeg = lastlength - 1;
    return;
}


int fmpz_mod_polyu3n_interp_crt_2sm_bpoly(
    slong * lastdeg,
    fmpz_mod_polyun_t F,
    fmpz_mod_polyun_t T,
    fmpz_mod_bpoly_t A,
    fmpz_mod_bpoly_t B,
    const fmpz_mod_poly_t modulus,
    fmpz_mod_poly_t alphapow,
    const fmpz_mod_ctx_t ctx)
{
    int changed = 0;
    slong lastlength = 0;
    fmpz_t zzero;
    fmpz_mod_poly_t zero;
    fmpz_mod_poly_struct * Tcoeffs;
    ulong * Texps;
    slong Ti;
    fmpz_mod_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    slong Flen = F->length;
    slong Fi;
    fmpz_mod_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    fmpz_mod_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    fmpz_mod_poly_struct * Fvalue;
    fmpz_t u, v, FvalueA, FvalueB;
    const fmpz * Avalue, * Bvalue;
    int texp_set, cmp;

    fmpz_init(u);
    fmpz_init(v);
    fmpz_init(FvalueA);
    fmpz_init(FvalueB);

#ifdef FLINT_WANT_ASSERT
    fmpz_mod_poly_evaluate_fmpz(u, modulus, alphapow->coeffs + 1, ctx);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx);
    fmpz_mod_add(u, u, u, ctx);
    FLINT_ASSERT(fmpz_is_one(u));
    fmpz_mod_neg(u, alphapow->coeffs + 1, ctx);
    fmpz_mod_poly_evaluate_fmpz(u, modulus, u, ctx);
    fmpz_mod_mul(u, u, alphapow->coeffs + 1, ctx);
    fmpz_mod_add(u, u, u, ctx);
    FLINT_ASSERT(fmpz_is_one(u));
#endif

    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(A, ctx));
    FLINT_ASSERT(fmpz_mod_bpoly_is_canonical(B, ctx));
    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(F, ctx));

    fmpz_mod_polyun_fit_length(T, FLINT_MAX(Flen, A->length), ctx);
    Tcoeffs = T->coeffs;
    Texps = T->exps;

    fmpz_init(zzero);
    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    Ti = Fi = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : fmpz_mod_poly_degree(A->coeffs + Ai, ctx);
    bi = (Bi < 0) ? 0 : fmpz_mod_poly_degree(B->coeffs + Bi, ctx);

    while (Fi < Flen || Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            extra = FLINT_MAX(extra, Bi);
            fmpz_mod_polyun_fit_length(T, Ti + extra + 1, ctx);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || Fcoeffs[Fi].length > 0);
        FLINT_ASSERT(Ai < 0 ||!fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
        FLINT_ASSERT(Bi < 0 ||!fmpz_is_zero(Bcoeffs[Bi].coeffs + bi));

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeffs + Fi;
            texp_set = 1;
            Texps[Ti] = Fexps[Fi];
        }

        Avalue = zzero;
        if (Ai >= 0)
        {
            ulong Aexp = pack_exp3(Ai, ai, 0);
            cmp = (!texp_set) ? -1 :
                  Texps[Ti] < Aexp ? -1 :
                  Texps[Ti] > Aexp ? 1 : 0;
            if (cmp <= 0)
            {
                Avalue = Acoeffs[Ai].coeffs + ai;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                Texps[Ti] = Aexp;
            }
        }

        Bvalue = zzero;
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            cmp = (!texp_set) ? -1 :
                  Texps[Ti] < Bexp ? -1 :
                  Texps[Ti] > Bexp ? 1 : 0;
            if (cmp <= 0)
            {
                Bvalue = Bcoeffs[Bi].coeffs + bi;
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = zzero;
                texp_set = 1;
                Texps[Ti] = Bexp;
            }
        }

        FLINT_ASSERT(texp_set);

        fmpz_mod_poly_eval2_pow(FvalueA, FvalueB, Fvalue, alphapow, ctx);
        fmpz_mod_sub(FvalueA, FvalueA, Avalue, ctx);
        fmpz_mod_sub(FvalueB, FvalueB, Bvalue, ctx);
        fmpz_mod_sub(u, FvalueB, FvalueA, ctx);
        fmpz_mod_add(v, FvalueB, FvalueA, ctx);
        fmpz_mod_mul(v, v, alphapow->coeffs + 1, ctx);
        fmpz_mod_neg(v, v, ctx);
        changed |= !fmpz_is_zero(u) || !fmpz_is_zero(v);
        fmpz_mod_poly_addmul_linear(Tcoeffs + Ti, Fvalue, modulus, u, v, ctx);

        FLINT_ASSERT(Tcoeffs[Ti].length >= Fvalue->length);

        Fi += (Fvalue != zero);
        if (Avalue != zzero)
        {
            do {
                ai--;
            } while (ai >= 0 && fmpz_is_zero(Acoeffs[Ai].coeffs + ai));
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = fmpz_mod_poly_degree(Acoeffs + Ai, ctx);
            }
        }
        if (Bvalue != zzero)
        {
            do {
                bi--;
            } while (bi >= 0 && fmpz_is_zero(Bcoeffs[Bi].coeffs + bi));
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = fmpz_mod_poly_degree(Bcoeffs + Bi, ctx);
            }
        }

        FLINT_ASSERT(!fmpz_mod_poly_is_zero(Tcoeffs + Ti, ctx));
        lastlength = FLINT_MAX(lastlength, Tcoeffs[Ti].length);
        Ti++;
    }
    T->length = Ti;

    if (changed)
        fmpz_mod_polyun_swap(T, F);

    FLINT_ASSERT(fmpz_mod_polyun_is_canonical(F, ctx));

    fmpz_clear(zzero);

    fmpz_clear(u);
    fmpz_clear(v);
    fmpz_clear(FvalueA);
    fmpz_clear(FvalueB);

    *lastdeg = lastlength - 1;
    return changed;
}

/*
    Input A, B0, B1, with A(y,x,z) = B0(y,x,z)*B1(y,x,z) mod (y-beta)

    return
       -1: suspect that B0(beta,x,z) & B1(beta,x,z) are not pairwise prime, i.e. could not find an eval point for z making them pairwise prime, or
           A(y,x,z) has wrong degree in x, or
           could not find enough eval points for z with the required degree in x
        0: lift of B0*B1 to true factors is impossible
        1: successfully lifted B0*B1 to true factors BB0*BB1 without changing lc_x
*/
/* r factor version */
int fmpz_mod_polyu3_hlift(
    slong r,
    fmpz_mod_polyun_struct * BB,
    fmpz_mod_polyu_t A,
    fmpz_mod_polyu_struct * B,
    const fmpz_t beta,
    slong degree_inner, /* required degree in x */
    const fmpz_mod_ctx_t ctx)
{
    int success/*, Eok*/;
    slong i, j;
    fmpz_mod_polyun_t T;
    fmpz_mod_bpoly_struct * Bp, * Bm;
    fmpz_mod_bpoly_t Ap, Am;
    fmpz_mod_poly_t modulus, alphapow, t1, t2;
    fmpz_t alpha, c;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;
    fmpz_mod_poly_bpoly_stack_t St;

/*
    if (r < 3)
        return n_polyu3_mod_hlift2(BB + 0, BB + 1, A, B + 0, B + 1,
                                                      beta, degree_inner, ctx);
*/
    FLINT_ASSERT(fmpz_mod_polyu_is_canonical(A, ctx));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(fmpz_mod_polyu_is_canonical(B + i, ctx));

    fmpz_init(alpha);
    fmpz_init(c);

    BBdegZ = FLINT_ARRAY_ALLOC(r, slong);
    Bp = FLINT_ARRAY_ALLOC(r, fmpz_mod_bpoly_struct);
    Bm = FLINT_ARRAY_ALLOC(r, fmpz_mod_bpoly_struct);
    for (i = 0; i < r; i++)
    {
        fmpz_mod_bpoly_init(Bp + i, ctx);
        fmpz_mod_bpoly_init(Bm + i, ctx);
    }
    fmpz_mod_polyun_init(T, ctx);
    fmpz_mod_bpoly_init(Ap, ctx);
    fmpz_mod_bpoly_init(Am, ctx);
    fmpz_mod_poly_init(modulus, ctx);
    fmpz_mod_poly_init(alphapow, ctx);
    fmpz_mod_poly_init(t1, ctx);
    fmpz_mod_poly_init(t2, ctx);
    fmpz_mod_poly_stack_init(St->poly_stack);
    fmpz_mod_bpoly_stack_init(St->bpoly_stack);

    fmpz_mod_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    fmpz_mod_poly_fit_length(alphapow, FLINT_MAX(WORD(3), AdegZ + 2), ctx);
    fmpz_mod_poly_one(modulus, ctx);

    FLINT_ASSERT(fmpz_is_odd(fmpz_mod_ctx_modulus(ctx)));

    fmpz_fdiv_q_2exp(alpha, fmpz_mod_ctx_modulus(ctx), 1);

    bad_primes_left = FLINT_MAX(5, AdegZ);

    goto choose_prime;

choose_prime:

    fmpz_sub_ui(alpha, alpha, 1);
    if (fmpz_sgn(alpha) <= 0)
    {
        success = -1;
        goto cleanup;
    }

    FLINT_ASSERT(fmpz_mod_is_canonical(alpha, ctx));
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    fmpz_one(alphapow->coeffs + 0);
    fmpz_set(alphapow->coeffs + 1, alpha);

    fmpz_mod_polyu3_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx);

    for (i = 0; i < r; i++)
    {
        fmpz_mod_polyu3_interp_reduce_2sm_bpoly(Bp + i, Bm + i,
                                                         B + i, alphapow, ctx);
    }

    success = fmpz_mod_bpoly_hlift(r, Ap, Bp, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = fmpz_mod_bpoly_hlift(r, Am, Bm, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    if (fmpz_mod_poly_degree(modulus, ctx) > 0)
    {
        fmpz_mod_poly_evaluate_fmpz(c, modulus, alpha, ctx);
    #if WANT_ASSERT
        {
            fmpz_t cc;
            fmpz_init(cc);
            fmpz_mod_neg(cc, alpha, ctx);
            fmpz_mod_poly_evaluate_fmpz(cc, modulus, cc, ctx);
            FLINT_ASSERT(fmpz_equal(cc, c));
            fmpz_clear(cc);
        }
    #endif

        fmpz_mod_mul(c, c, alpha, ctx);
        fmpz_mod_add(c, c, c, ctx);
        fmpz_mod_inv(c, c, ctx);
        fmpz_mod_poly_scalar_mul_fmpz(modulus, modulus, c, ctx);
        for (i = 0; i < r; i++)
        {
            fmpz_mod_polyu3n_interp_crt_2sm_bpoly(BBdegZ + i, BB + i, T,
                                  Bp + i, Bm + i, modulus, alphapow, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            fmpz_mod_polyu3n_interp_lift_2sm_bpoly(BBdegZ + i, BB + i,
                                               Bp + i, Bm + i, alpha, ctx);
        }
    }

    fmpz_mod_mul(c, alpha, alpha, ctx);
    fmpz_mod_neg(c, c, ctx);
    fmpz_mod_poly_shift_left_scalar_addmul_fmpz_mod(modulus, 2, c, ctx);

    j = BBdegZ[0];
    for (i = 1; i < r; i++)
        j += BBdegZ[i];

    if (j > AdegZ)
    {
        success = 0;
        goto cleanup;
    }

    if (fmpz_mod_poly_degree(modulus, ctx) <= AdegZ)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#ifdef FLINT_WANT_ASSERT
    if (success == 1)
    {
        fmpz_mod_polyu_t T1, T2, T3;
        fmpz_mod_polyu_init(T1);
        fmpz_mod_polyu_init(T2);
        fmpz_mod_polyu_init(T3);
        fmpz_mod_polyu_get_fmpz_mod_polyun(T2, BB + 0, ctx);
        fmpz_mod_polyu_get_fmpz_mod_polyun(T3, BB + 1, ctx);
        fmpz_mod_polyu_mul(T1, T2, T3, ctx);
        for (i = 2; i < r; i++)
        {
            fmpz_mod_polyu_get_fmpz_mod_polyun(T3, BB + i, ctx);
            fmpz_mod_polyu_mul(T2, T1, T3, ctx);
            fmpz_mod_polyu_swap(T2, T1);
        }
        FLINT_ASSERT(fmpz_mod_polyu_equal(A, T1));
        fmpz_mod_polyu_clear(T1);
        fmpz_mod_polyu_clear(T2);
        fmpz_mod_polyu_clear(T3);
    }
#endif

    fmpz_mod_polyun_clear(T, ctx);
    fmpz_mod_bpoly_clear(Ap, ctx);
    fmpz_mod_bpoly_clear(Am, ctx);

    for (i = 0; i < r; i++)
    {
        fmpz_mod_bpoly_clear(Bp + i, ctx);
        fmpz_mod_bpoly_clear(Bm + i, ctx);
    }

    flint_free(BBdegZ);
    flint_free(Bp);
    flint_free(Bm);

    fmpz_mod_poly_clear(modulus, ctx);
    fmpz_mod_poly_clear(alphapow, ctx);
    fmpz_mod_poly_clear(t1, ctx);
    fmpz_mod_poly_clear(t2, ctx);

    fmpz_mod_poly_stack_clear(St->poly_stack);
    fmpz_mod_bpoly_stack_clear(St->bpoly_stack);

    fmpz_clear(alpha);
    fmpz_clear(c);

    FLINT_ASSERT(success || fmpz_abs_fits_ui(fmpz_mod_ctx_modulus(ctx)));

    return success;
}
