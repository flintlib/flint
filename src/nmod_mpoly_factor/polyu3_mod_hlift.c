/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


#ifdef FLINT_WANT_ASSERT

/* take last variable out of dense storage */
static void n_polyu_set_n_polyun(n_polyu_t A, const n_polyun_t B)
{
    slong i, j;
    A->length = 0;
    for (i = 0; i < B->length; i++)
    {
        for (j = B->coeffs[i].length - 1; j >= 0; j--)
        {
            if (B->coeffs[i].coeffs[j] == 0)
                continue;
            n_polyu_fit_length(A, A->length + 1);
            A->coeffs[A->length] = B->coeffs[i].coeffs[j];
            A->exps[A->length] = B->exps[i] + j;
            A->length++;
        }
    }
}


static void n_polyu_sort_terms(n_polyu_t A)
{
    slong i, j;
    for (i = 1; i < A->length; i++)
    for (j = i; j > 0 && A->exps[j - 1] < A->exps[j]; j--)
    {
        FLINT_SWAP(mp_limb_t, A->coeffs[j - 1], A->coeffs[j]);
        FLINT_SWAP(ulong, A->exps[j - 1], A->exps[j]);
    }
    return;
}


int n_polyu_mod_is_canonical(const n_polyu_t A, nmod_t mod)
{
    slong i;
    if (A->length < 0)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i] >= mod.n || A->coeffs[i] <= 0)
            return 0;
        if (i > 0 && A->exps[i] >= A->exps[i - 1])
            return 0;
    }
    return 1;
}

static void n_polyu_mod_combine_like_terms(n_polyu_t A, nmod_t ctx)
{
    slong in, out;

    out = -1;

    for (in = 0; in < A->length; in++)
    {
        FLINT_ASSERT(in > out);

        if (out >= 0 && A->exps[out] == A->exps[in])
        {
            A->coeffs[out] = nmod_add(A->coeffs[out], A->coeffs[in], ctx);
        }
        else
        {
            if (out < 0 || A->coeffs[out] != 0)
                out++;

            if (out != in)
            {
                A->exps[out] = A->exps[in];
                A->coeffs[out] = A->coeffs[in];
            }
        }
    }

    if (out < 0 || A->coeffs[out] != 0)
        out++;

    A->length = out;
}


static int n_polyu_equal(n_polyu_t A, n_polyu_t B)
{
    slong i;
    if (A->length != B->length)
        return 0;
    for (i = 0; i < A->length; i++)
    {
        if (A->coeffs[i] != B->coeffs[i])
            return 0;
        if (A->exps[i] != B->exps[i])
            return 0;
    }
    return 1;
}


static void n_polyu_mod_mul(
    n_polyu_t A,
    const n_polyu_t B,
    const n_polyu_t C,
    nmod_t ctx)
{
    slong Ai, Bi, Ci;

    FLINT_ASSERT(A != B);
    FLINT_ASSERT(A != C);

    n_polyu_fit_length(A, B->length*C->length);
    Ai = 0;
    for (Bi = 0; Bi < B->length; Bi++)
    for (Ci = 0; Ci < C->length; Ci++)
    {
        A->exps[Ai]   = B->exps[Bi] + C->exps[Ci];
        A->coeffs[Ai] = nmod_mul(B->coeffs[Bi], C->coeffs[Ci], ctx);
        Ai++;
    }
    A->length = Ai;
    n_polyu_sort_terms(A);
    n_polyu_mod_combine_like_terms(A, ctx);
}

#endif


void n_poly_fill_powers(
    n_poly_t alphapow,
    slong target,
    nmod_t mod)
{
    if (target + 1 > alphapow->length)
    {
        slong k;
        slong oldlength = alphapow->length;
        FLINT_ASSERT(2 <= oldlength);
        n_poly_fit_length(alphapow, target + 1);
        for (k = oldlength; k <= target; k++)
        {
            alphapow->coeffs[k] = nmod_mul(alphapow->coeffs[k - 1],
                                                alphapow->coeffs[1], mod);
        }
        alphapow->length = target + 1;
    }
}


void n_polyu3_mod_interp_reduce_2sm_bpoly(
    n_bpoly_t Ap,
    n_bpoly_t Am,
    const n_polyu_t A,
    n_poly_t alpha,
    nmod_t mod)
{
    slong i;
    slong cur0, cur1, e0, e1, e2;
    ulong tp0, tp1, tp2, tm0, tm1, tm2, p1, p0;
    const mp_limb_t * Acoeffs = A->coeffs;
    const ulong * Aexps = A->exps;

    n_bpoly_zero(Ap);
    n_bpoly_zero(Am);

    FLINT_ASSERT(A->length > 0);

    i = 0;

    cur0 = extract_exp(Aexps[i], 2, 3);
    cur1 = extract_exp(Aexps[i], 1, 3);
    e2   = extract_exp(Aexps[i], 0, 3);

    n_poly_fill_powers(alpha, e2, mod);

    tp2 = tp1 = tp0 = tm2 = tm1 = tm0 = 0;
    if (e2 & 1)
        umul_ppmm(tm1, tm0, alpha->coeffs[e2], Acoeffs[i]);
    else
        umul_ppmm(tp1, tp0, alpha->coeffs[e2], Acoeffs[i]);

    for (i = 1; i < A->length; i++)
    {
        e0 = extract_exp(Aexps[i], 2, 3);
        e1 = extract_exp(Aexps[i], 1, 3);
        e2 = extract_exp(Aexps[i], 0, 3);

        FLINT_ASSERT(e0 <= cur0);
        if (e0 < cur0 || e1 < cur1)
        {
            NMOD_RED3(tp0, tp2, tp1, tp0, mod);
            NMOD_RED3(tm0, tm2, tm1, tm0, mod);

            n_bpoly_set_coeff(Ap, cur0, cur1, nmod_add(tp0, tm0, mod));
            n_bpoly_set_coeff(Am, cur0, cur1, nmod_sub(tp0, tm0, mod));

            tp2 = tp1 = tp0 = tm2 = tm1 = tm0 = 0;
        }
        else
        {
            FLINT_ASSERT(e0 == cur0);
            FLINT_ASSERT(e1 == cur1);
        }

        cur0 = e0;
        cur1 = e1;

        n_poly_fill_powers(alpha, e2, mod);

        FLINT_ASSERT(alpha->coeffs[e2] == nmod_pow_ui(alpha->coeffs[1], e2, mod));

        umul_ppmm(p1, p0, alpha->coeffs[e2], Acoeffs[i]);
        if (e2 & 1)
            add_sssaaaaaa(tm2, tm1, tm0, tm2, tm1, tm0, 0, p1, p0);
        else
            add_sssaaaaaa(tp2, tp1, tp0, tp2, tp1, tp0, 0, p1, p0);
    }

    NMOD_RED3(tp0, tp2, tp1, tp0, mod);
    NMOD_RED3(tm0, tm2, tm1, tm0, mod);
    n_bpoly_set_coeff(Ap, cur0, cur1, nmod_add(tp0, tm0, mod));
    n_bpoly_set_coeff(Am, cur0, cur1, nmod_sub(tp0, tm0, mod));
}


/*
    T(x0, x1, x2) is in F[x2][x0, x1]
    A(x0, x1) are B(x0, x1) are in F[x0, x1]
    set T so that
        T(x0, x1, x2) == A(x0, x1) mod (x2 - alpha)
        T(x0, x1, x2) == B(x0, x1) mod (x2 + alpha)
*/

void n_polyu3n_mod_interp_lift_2sm_bpoly(
    slong * lastdeg,
    n_polyun_t T,
    const n_bpoly_t A,
    const n_bpoly_t B,
    mp_limb_t alpha,
    nmod_t mod)
{
    slong lastlength = 0;
    n_poly_struct * Tcoeffs;
    ulong * Texps;
    slong Ti;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    n_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    mp_limb_t u, v, Avalue, Bvalue;
    mp_limb_t d0, d1;

    FLINT_ASSERT(2*alpha < mod.n);

    d0 = (1 + mod.n)/2;
    d1 = nmod_inv(nmod_add(alpha, alpha, mod), mod);

    n_polyun_fit_length(T, FLINT_MAX(A->length, B->length));
    Tcoeffs = T->coeffs;
    Texps = T->exps;

    Ti = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);
    bi = (Bi < 0) ? 0 : n_poly_degree(B->coeffs + Bi);

    while (Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            n_polyun_fit_length(T, Ti + FLINT_MAX(Ai, Bi) + 1);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        FLINT_ASSERT(Ai < 0 || Acoeffs[Ai].coeffs[ai] != 0);
        FLINT_ASSERT(Bi < 0 || Bcoeffs[Bi].coeffs[bi] != 0);

        Avalue = 0;
        if (Ai >= 0)
        {
            Avalue = Acoeffs[Ai].coeffs[ai];
            Texps[Ti] = pack_exp3(Ai, ai, 0);
        }

        Bvalue = 0;
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            if (Avalue == 0)
            {
                Bvalue = Bcoeffs[Bi].coeffs[bi];
                Texps[Ti] = Bexp;
            }
            else
            {
                if (Texps[Ti] <= Bexp)
                {
                    Bvalue = Bcoeffs[Bi].coeffs[bi];
                }
                if (Texps[Ti] < Bexp)
                {
                    Avalue = 0;
                    Texps[Ti] = Bexp;
                }
            }
        }

        u = nmod_sub(Avalue, Bvalue, mod);
        v = nmod_add(Avalue, Bvalue, mod);
        u = nmod_mul(u, d1, mod);
        v = nmod_mul(v, d0, mod);

        FLINT_ASSERT(u != 0 || v != 0);

        n_poly_fit_length(Tcoeffs + Ti, 2);
        Tcoeffs[Ti].coeffs[0] = v;
        Tcoeffs[Ti].coeffs[1] = u;
        Tcoeffs[Ti].length = 1 + (u != 0);
        lastlength = FLINT_MAX(lastlength, Tcoeffs[Ti].length);
        Ti++;

        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && Bcoeffs[Bi].coeffs[bi] == 0);
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = n_poly_degree(Bcoeffs + Bi);
            }
        }
    }
    T->length = Ti;

    FLINT_ASSERT(n_polyun_mod_is_canonical(T, mod));

    *lastdeg = lastlength - 1;

    return;
}


int n_polyu3n_mod_interp_crt_2sm_bpoly(
    slong * lastdeg,
    n_polyun_t F,
    n_polyun_t T,
    n_bpoly_t A,
    n_bpoly_t B,
    const n_poly_t modulus,
    n_poly_t alphapow,
    nmod_t mod)
{
    int changed = 0;
    slong lastlength = 0;
    n_poly_t zero;
    n_poly_struct * Tcoeffs;
    ulong * Texps;
    slong Ti;
    n_poly_struct * Fcoeffs = F->coeffs;
    ulong * Fexps = F->exps;
    slong Flen = F->length;
    slong Fi;
    n_poly_struct * Acoeffs = A->coeffs;
    slong Ai, ai;
    n_poly_struct * Bcoeffs = B->coeffs;
    slong Bi, bi;
    n_poly_struct * Fvalue;
    mp_limb_t u, v, Avalue, Bvalue, FvalueA, FvalueB;
    int texp_set, cmp;
    mp_limb_t alpha = alphapow->coeffs[1];

#ifdef FLINT_WANT_ASSERT
    u = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
    u = nmod_mul(u, alpha, mod);
    u = nmod_mul(u, 2, mod);
    FLINT_ASSERT(u == 1);
    u = n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod);
    u = nmod_mul(u, alpha, mod);
    u = nmod_mul(u, 2, mod);
    FLINT_ASSERT(u == 1);
#endif

    FLINT_ASSERT(n_bpoly_mod_is_canonical(A, mod));
    FLINT_ASSERT(n_bpoly_mod_is_canonical(B, mod));
    FLINT_ASSERT(n_polyun_mod_is_canonical(F, mod));

    n_polyun_fit_length(T, FLINT_MAX(Flen, A->length));
    Tcoeffs = T->coeffs;
    Texps = T->exps;

    zero->alloc = 0;
    zero->length = 0;
    zero->coeffs = NULL;

    Ti = Fi = 0;
    Ai = A->length - 1;
    Bi = B->length - 1;
    ai = (Ai < 0) ? 0 : n_poly_degree(A->coeffs + Ai);
    bi = (Bi < 0) ? 0 : n_poly_degree(B->coeffs + Bi);

    while (Fi < Flen || Ai >= 0 || Bi >= 0)
    {
        if (Ti >= T->alloc)
        {
            slong extra = Flen - Fi;
            extra = FLINT_MAX(extra, Ai);
            extra = FLINT_MAX(extra, Bi);
            n_polyun_fit_length(T, Ti + extra + 1);
            Tcoeffs = T->coeffs;
            Texps = T->exps;
        }

        FLINT_ASSERT(Fi >= Flen || Fcoeffs[Fi].length > 0);
        FLINT_ASSERT(Ai < 0 || Acoeffs[Ai].coeffs[ai] != 0);
        FLINT_ASSERT(Bi < 0 || Bcoeffs[Bi].coeffs[bi] != 0);

        Fvalue = zero;
        texp_set = 0;
        if (Fi < Flen)
        {
            Fvalue = Fcoeffs + Fi;
            texp_set = 1;
            Texps[Ti] = Fexps[Fi];
        }

        Avalue = 0;
        if (Ai >= 0)
        {
            ulong Aexp = pack_exp3(Ai, ai, 0);
            cmp = (!texp_set) ? -1 :
                  Texps[Ti] < Aexp ? -1 :
                  Texps[Ti] > Aexp ? 1 : 0;
            if (cmp <= 0)
            {
                Avalue = Acoeffs[Ai].coeffs[ai];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                texp_set = 1;
                Texps[Ti] = Aexp;
            }
        }

        Bvalue = 0;
        if (Bi >= 0)
        {
            ulong Bexp = pack_exp3(Bi, bi, 0);
            cmp = (!texp_set) ? -1 :
                  Texps[Ti] < Bexp ? -1 :
                  Texps[Ti] > Bexp ? 1 : 0;
            if (cmp <= 0)
            {
                Bvalue = Bcoeffs[Bi].coeffs[bi];
            }
            if (cmp < 0)
            {
                Fvalue = zero;
                Avalue = 0;
                texp_set = 1;
                Texps[Ti] = Bexp;
            }
        }

        FLINT_ASSERT(texp_set);

        n_poly_mod_eval2_pow(&FvalueA, &FvalueB, Fvalue, alphapow, mod);
        FvalueA = nmod_sub(FvalueA, Avalue, mod);
        FvalueB = nmod_sub(FvalueB, Bvalue, mod);
        u = nmod_sub(FvalueB, FvalueA, mod);
        v = nmod_mul(mod.n - alpha, nmod_add(FvalueB, FvalueA, mod), mod);
        if (u != 0 || v != 0)
        {
            changed = 1;
            n_poly_mod_addmul_linear(Tcoeffs + Ti, Fvalue, modulus, u, v, mod);
        }
        else
        {
            FLINT_ASSERT(!n_poly_is_zero(Fvalue));
            n_poly_set(Tcoeffs + Ti, Fvalue);
        }

        FLINT_ASSERT(Tcoeffs[Ti].length >= Fvalue->length);

        Fi += (Fvalue != zero);
        if (Avalue != 0)
        {
            do {
                ai--;
            } while (ai >= 0 && Acoeffs[Ai].coeffs[ai] == 0);
            if (ai < 0)
            {
                do {
                    Ai--;
                } while (Ai >= 0 && Acoeffs[Ai].length == 0);
                if (Ai >= 0)
                    ai = n_poly_degree(Acoeffs + Ai);
            }
        }
        if (Bvalue != 0)
        {
            do {
                bi--;
            } while (bi >= 0 && Bcoeffs[Bi].coeffs[bi] == 0);
            if (bi < 0)
            {
                do {
                    Bi--;
                } while (Bi >= 0 && Bcoeffs[Bi].length == 0);
                if (Bi >= 0)
                    bi = n_poly_degree(Bcoeffs + Bi);
            }
        }

        FLINT_ASSERT(!n_poly_is_zero(Tcoeffs + Ti));
        lastlength = FLINT_MAX(lastlength, Tcoeffs[Ti].length);
        Ti++;
    }
    T->length = Ti;

    if (changed)
        n_polyun_swap(T, F);

    FLINT_ASSERT(n_polyun_mod_is_canonical(F, mod));

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
int n_polyu3_mod_hlift2(
    n_polyun_t BB0,
    n_polyun_t BB1,
    n_polyu_t A,
    n_polyu_t B0,
    n_polyu_t B1,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    nmod_t ctx)
{
    int success, Eok;
    n_polyun_t T;
    n_bpoly_t Ap, Am, B0p, B0m, B1p, B1m;
    n_poly_t modulus, alphapow, t1, t2;
    mp_limb_t alpha, c;
    slong ldegBB0, ldegBB1;
    slong Adegy, Adegz, Adegx;
    slong bad_primes_left;
    n_poly_bpoly_stack_t St;
    nmod_eval_interp_t E;

    FLINT_ASSERT(n_polyu_mod_is_canonical(A, ctx));
    FLINT_ASSERT(n_polyu_mod_is_canonical(B0, ctx));
    FLINT_ASSERT(n_polyu_mod_is_canonical(B1, ctx));

    n_polyun_init(T);
    n_bpoly_init(Ap);
    n_bpoly_init(Am);
    n_bpoly_init(B0p);
    n_bpoly_init(B0m);
    n_bpoly_init(B1p);
    n_bpoly_init(B1m);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);
    nmod_eval_interp_init(E);

    Eok = nmod_eval_interp_set_degree_modulus(E, degree_inner, ctx);

    n_polyu3_degrees(&Adegy, &Adegx, &Adegz, A);

    if (Adegx != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), Adegz + 2));
    n_poly_one(modulus);

    FLINT_ASSERT((ctx.n & UWORD(1)) == UWORD(1));

    alpha = (ctx.n - UWORD(1))/UWORD(2);

    bad_primes_left = FLINT_MAX(5, Adegz);

    goto choose_prime;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_polyu3_mod_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx);
    n_polyu3_mod_interp_reduce_2sm_bpoly(B0p, B0m, B0, alphapow, ctx);
    n_polyu3_mod_interp_reduce_2sm_bpoly(B1p, B1m, B1, alphapow, ctx);

    success = Eok ? n_bpoly_mod_hlift2_cubic(Ap, B0p, B1p, beta, degree_inner, ctx, E, St) :
                    n_bpoly_mod_hlift2(Ap, B0p, B1p, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = Eok ? n_bpoly_mod_hlift2_cubic(Am, B0m, B1m, beta, degree_inner, ctx, E, St) :
                    n_bpoly_mod_hlift2(Am, B0m, B1m, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, ctx);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus,
                                                ctx.n - alpha, ctx));
        c = nmod_mul(c, alpha, ctx);
        c = nmod_add(c, c, ctx);
        c = n_invmod(c, ctx.n);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx);
        n_polyu3n_mod_interp_crt_2sm_bpoly(&ldegBB0, BB0, T, B0p, B0m,
                                                  modulus, alphapow, ctx);
        n_polyu3n_mod_interp_crt_2sm_bpoly(&ldegBB1, BB1, T, B1p, B1m,
                                                  modulus, alphapow, ctx);
    }
    else
    {
        n_polyu3n_mod_interp_lift_2sm_bpoly(&ldegBB0, BB0, B0p, B0m, alpha, ctx);
        n_polyu3n_mod_interp_lift_2sm_bpoly(&ldegBB1, BB1, B1p, B1m, alpha, ctx);
    }

    c = ctx.n - nmod_mul(alpha, alpha, ctx);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, ctx);

    if (ldegBB0 + ldegBB1 > Adegz)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) <= Adegz)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#ifdef FLINT_WANT_ASSERT
    if (success == 1)
    {
        n_polyu_t T1, T2, T3;
        n_polyu_init(T1);
        n_polyu_init(T2);
        n_polyu_init(T3);
        n_polyu_set_n_polyun(T2, BB0);
        n_polyu_set_n_polyun(T3, BB1);
        n_polyu_mod_mul(T1, T2, T3, ctx);
        FLINT_ASSERT(n_polyu_equal(A, T1));
        n_polyu_clear(T1);
        n_polyu_clear(T2);
        n_polyu_clear(T3);
    }
#endif

    n_polyun_clear(T);
    n_bpoly_clear(Ap);
    n_bpoly_clear(Am);
    n_bpoly_clear(B0p);
    n_bpoly_clear(B0m);
    n_bpoly_clear(B1p);
    n_bpoly_clear(B1m);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);
    nmod_eval_interp_clear(E);

    return success;
}


/* r factor version */
int n_polyu3_mod_hlift(
    slong r,
    n_polyun_struct * BB,
    n_polyu_t A,
    n_polyu_struct * B,
    mp_limb_t beta,
    slong degree_inner, /* required degree in x */
    nmod_t ctx)
{
    int success, Eok;
    slong i, j;
    n_polyun_t T;
    n_bpoly_struct * Bp, * Bm;
    n_bpoly_t Ap, Am;
    n_poly_t modulus, alphapow, t1, t2;
    mp_limb_t alpha, c;
    slong * BBdegZ;
    slong AdegY, AdegX, AdegZ;
    slong bad_primes_left;
    n_poly_bpoly_stack_t St;
    nmod_eval_interp_t E;

    if (r < 3)
        return n_polyu3_mod_hlift2(BB + 0, BB + 1, A, B + 0, B + 1,
                                                      beta, degree_inner, ctx);

    FLINT_ASSERT(n_polyu_mod_is_canonical(A, ctx));
    for (i = 0; i < r; i++)
        FLINT_ASSERT(n_polyu_mod_is_canonical(B + i, ctx));

    BBdegZ = (slong *) flint_malloc(r*sizeof(slong));
    Bp = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    Bm = (n_bpoly_struct *) flint_malloc(r*sizeof(n_bpoly_struct));
    for (i = 0; i < r; i++)
    {
        n_bpoly_init(Bp + i);
        n_bpoly_init(Bm + i);
    }
    n_polyun_init(T);
    n_bpoly_init(Ap);
    n_bpoly_init(Am);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);
    n_poly_stack_init(St->poly_stack);
    n_bpoly_stack_init(St->bpoly_stack);
    nmod_eval_interp_init(E);

    Eok = nmod_eval_interp_set_degree_modulus(E, degree_inner, ctx);

    n_polyu3_degrees(&AdegY, &AdegX, &AdegZ, A);
    if (AdegX != degree_inner)
    {
        success = -1;
        goto cleanup;
    }

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), AdegZ + 2));
    n_poly_one(modulus);

    FLINT_ASSERT((ctx.n & UWORD(1)) == UWORD(1));

    alpha = (ctx.n - UWORD(1))/UWORD(2);

    bad_primes_left = FLINT_MAX(5, AdegZ);

    goto choose_prime;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }
    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= ctx.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_polyu3_mod_interp_reduce_2sm_bpoly(Ap, Am, A, alphapow, ctx);

    for (i = 0; i < r; i++)
    {
        n_polyu3_mod_interp_reduce_2sm_bpoly(Bp + i, Bm + i,
                                                    B + i, alphapow, ctx);
    }

    success = Eok ? n_bpoly_mod_hlift_cubic(r, Ap, Bp, beta, degree_inner, ctx, E, St) :
                    n_bpoly_mod_hlift(r, Ap, Bp, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    success = Eok ? n_bpoly_mod_hlift_cubic(r, Am, Bm, beta, degree_inner, ctx, E, St) :
                    n_bpoly_mod_hlift(r, Am, Bm, beta, degree_inner, ctx, St);
    if (success < 1)
    {
        if (success == 0 || --bad_primes_left < 0)
            goto cleanup;
        goto choose_prime;
    }

    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, ctx);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus,
                                                ctx.n - alpha, ctx));
        c = nmod_mul(c, alpha, ctx);
        c = nmod_add(c, c, ctx);
        c = nmod_inv(c, ctx);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, ctx);
        for (i = 0; i < r; i++)
        {
            n_polyu3n_mod_interp_crt_2sm_bpoly(BBdegZ + i, BB + i, T,
                                  Bp + i, Bm + i, modulus, alphapow, ctx);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            n_polyu3n_mod_interp_lift_2sm_bpoly(BBdegZ + i, BB + i,
                                               Bp + i, Bm + i, alpha, ctx);
        }
    }

    c = ctx.n - nmod_mul(alpha, alpha, ctx);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, ctx);

    j = BBdegZ[0];
    for (i = 1; i < r; i++)
        j += BBdegZ[i];

    if (j > AdegZ)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) <= AdegZ)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#ifdef FLINT_WANT_ASSERT
    if (success == 1)
    {
        n_polyu_t T1, T2, T3;
        n_polyu_init(T1);
        n_polyu_init(T2);
        n_polyu_init(T3);
        n_polyu_set_n_polyun(T2, BB + 0);
        n_polyu_set_n_polyun(T3, BB + 1);
        n_polyu_mod_mul(T1, T2, T3, ctx);
        for (i = 2; i < r; i++)
        {
            n_polyu_set_n_polyun(T3, BB + i);
            n_polyu_mod_mul(T2, T1, T3, ctx);
            n_polyu_swap(T2, T1);
        }
        FLINT_ASSERT(n_polyu_equal(A, T1));
        n_polyu_clear(T1);
        n_polyu_clear(T2);
        n_polyu_clear(T3);
    }
#endif

    n_polyun_clear(T);
    n_bpoly_clear(Ap);
    n_bpoly_clear(Am);

    for (i = 0; i < r; i++)
    {
        n_bpoly_clear(Bp + i);
        n_bpoly_clear(Bm + i);
    }

    flint_free(BBdegZ);
    flint_free(Bp);
    flint_free(Bm);

    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);
    n_poly_stack_clear(St->poly_stack);
    n_bpoly_stack_clear(St->bpoly_stack);
    nmod_eval_interp_clear(E);

    return success;
}

