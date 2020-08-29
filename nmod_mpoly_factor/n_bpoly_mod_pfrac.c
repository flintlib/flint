/*
    Copyright (C) 2020 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "nmod_mpoly_factor.h"


int n_bpoly_mod_pfrac2(
    n_bpoly_t C1, n_bpoly_t C2,
    slong C1_deg1_bound, slong C2_deg1_bound,
    n_bpoly_t A,
    n_bpoly_t B1, n_bpoly_t B2,
    nmod_t mod)
{
    int success;
    slong A_deg1, B1_deg1, B2_deg1, C1_deg1, C2_deg1;
    slong bad_prime_count, bound;
    mp_limb_t alpha, c;
    n_poly_t Aevalp, B1evalp, B2evalp, C1evalp, C2evalp;
    n_poly_t Aevalm, B1evalm, B2evalm, C1evalm, C2evalm;
    n_poly_t modulus, alphapow, t1, t2;
    n_bpoly_t T;

    n_poly_init(Aevalp);
    n_poly_init(B1evalp);
    n_poly_init(B2evalp);
    n_poly_init(C1evalp);
    n_poly_init(C2evalp);
    n_poly_init(Aevalm);
    n_poly_init(B1evalm);
    n_poly_init(B2evalm);
    n_poly_init(C1evalm);
    n_poly_init(C2evalm);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);

    n_bpoly_init(T);

    A_deg1 = n_bpoly_degree1(A);
    B1_deg1 = n_bpoly_degree1(B1);
    B2_deg1 = n_bpoly_degree1(B2);

    bound = A_deg1;
    bound = FLINT_MAX(bound, C1_deg1_bound + B2_deg1);
    bound = FLINT_MAX(bound, B1_deg1 + C2_deg1_bound);
    bound += 1;

    bad_prime_count = 0;

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    if ((mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > bound)
    {
        success = n_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, mod);
    n_bpoly_mod_interp_reduce_2sm_poly(B1evalp, B1evalm, B1, alphapow, mod);
    n_bpoly_mod_interp_reduce_2sm_poly(B2evalp, B2evalm, B2, alphapow, mod);

    /* make sure evaluation point did not drop the degree of a Bi */
    if (n_poly_degree(B1evalp) < n_bpoly_degree0(B1) ||
        n_poly_degree(B1evalm) < n_bpoly_degree0(B1) ||
        n_poly_degree(B2evalp) < n_bpoly_degree0(B2) ||
        n_poly_degree(B2evalm) < n_bpoly_degree0(B2))
    {
        goto choose_prime;
    }

    /* image pfrac's */
    if (!n_poly_mod_invmod(t1, B2evalp, B1evalp, mod))
        goto bad_prime;
    _n_poly_mod_mul(t2, Aevalp, t1, mod);
    _n_poly_mod_rem(C1evalp, t2, B1evalp, mod);
    _n_poly_mod_mul(t2, B2evalp, C1evalp, mod);
    n_poly_mod_sub(Aevalp, Aevalp, t2, mod);
    _n_poly_mod_div(C2evalp, Aevalp, B1evalp, mod);

    if (!n_poly_mod_invmod(t1, B2evalm, B1evalm, mod))
        goto bad_prime;
    _n_poly_mod_mul(t2, Aevalm, t1, mod);
    _n_poly_mod_rem(C1evalm, t2, B1evalm, mod);
    _n_poly_mod_mul(t2, B2evalm, C1evalm, mod);
    n_poly_mod_sub(Aevalm, Aevalm, t2, mod);
    _n_poly_mod_div(C2evalm, Aevalm, B1evalm, mod);

    /* update interpolants */
    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod));
        c = nmod_mul(c, alpha, mod);
        c = nmod_add(c, c, mod);
        c = nmod_inv(c, mod);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

        n_bpoly_mod_interp_crt_2sm_poly(&C1_deg1, C1, T, C1evalp, C1evalm, modulus, alphapow, mod);
        n_bpoly_mod_interp_crt_2sm_poly(&C2_deg1, C2, T, C2evalp, C2evalm, modulus, alphapow, mod);
    }
    else
    {
        n_bpoly_mod_interp_lift_2sm_poly(&C1_deg1, C1, C1evalp, C1evalm, alpha, mod);
        n_bpoly_mod_interp_lift_2sm_poly(&C2_deg1, C2, C2evalp, C2evalm, alpha, mod);
    }

    c = mod.n - nmod_mul(alpha, alpha, mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, mod);

    if (C1_deg1 > C1_deg1_bound || C2_deg1 > C2_deg1_bound)
    {
        success = 0;
        goto cleanup;
    }

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success == 1)
    {
        n_bpoly_t T1, T2, T3;
        n_bpoly_init(T1);
        n_bpoly_init(T2);
        n_bpoly_init(T3);
        n_bpoly_set(T1, A);
        n_bpoly_mod_mul(T2, C1, B2, mod);
        n_bpoly_mod_sub(T3, T1, T2, mod);
        n_bpoly_swap(T1, T3);
        n_bpoly_mod_mul(T2, C2, B1, mod);
        n_bpoly_mod_sub(T3, T1, T2, mod);
        n_bpoly_swap(T1, T3);
        FLINT_ASSERT(n_bpoly_is_zero(T1));
        n_bpoly_clear(T1);
        n_bpoly_clear(T2);
        n_bpoly_clear(T3);
    }
#endif

    n_poly_clear(Aevalp);
    n_poly_clear(B1evalp);
    n_poly_clear(B2evalp);
    n_poly_clear(C1evalp);
    n_poly_clear(C2evalp);
    n_poly_clear(Aevalm);
    n_poly_clear(B1evalm);
    n_poly_clear(B2evalm);
    n_poly_clear(C1evalm);
    n_poly_clear(C2evalm);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);

    n_bpoly_clear(T);

    return success;
}


/*
    Try to solve A/(B[0]*...*B[r-1]) = C[0]/B[0] + ... + C[r-1]/B[r-1]
    for the C[i] in Fp[y][x].
    return:
        1: solution found with deg_y(Ci) <= ldegCibound
        0: no solution exists with deg_y(Ci) <= ldegCibound
       -1: could not find enough evaluation points where the Bi are pariwise prime
       -2: found no evaluation points where the Bi are pariwise prime
*/
int n_bpoly_mod_pfrac(
    slong r,
    n_bpoly_struct * C,
    slong * C_deg1_bound,
    n_bpoly_t A,
    n_bpoly_struct * B,
    nmod_t mod)
{
    int success;
    slong i, j, bad_prime_count, bound;
    mp_limb_t alpha, c;
    n_poly_struct Aevalp[1], * Bevalp, * Cevalp;
    n_poly_struct Aevalm[1], * Bevalm, * Cevalm;
    n_poly_t modulus, alphapow, t1, t2;
    n_bpoly_t T;
    slong * B_deg1, * C_deg1, B_deg1_total, A_deg1;
    TMP_INIT;

    if (r < 3)
    {
        return n_bpoly_mod_pfrac2(C + 0, C + 1, C_deg1_bound[0],
                                        C_deg1_bound[1], A, B + 0, B + 1, mod);
    }

    TMP_START;

    B_deg1 = (slong *) TMP_ALLOC(2*r*sizeof(slong));
    C_deg1 = B_deg1 + r;

    Bevalp = (n_poly_struct *) TMP_ALLOC(4*r*sizeof(n_poly_struct));
    Bevalm = Bevalp + r;
    Cevalp = Bevalm + r;
    Cevalm = Cevalp + r;

    n_poly_init(Aevalp);
    n_poly_init(Aevalm);
    n_poly_init(modulus);
    n_poly_init(alphapow);
    n_poly_init(t1);
    n_poly_init(t2);
    for (i = 0; i < r; i++)
    {
        n_poly_init(Bevalp + i);
        n_poly_init(Bevalm + i);
        n_poly_init(Cevalp + i);
        n_poly_init(Cevalm + i);
    }

    n_bpoly_init(T);

    A_deg1 = n_bpoly_degree1(A);
    B_deg1_total = 0;
    for (i = 0; i < r; i++)
    {
        B_deg1[i] = n_bpoly_degree1(B + i);
        B_deg1_total += B_deg1[i];
    }

    bound = A_deg1;
    for (i = 0; i < r; i++)
        bound = FLINT_MAX(bound, C_deg1_bound[i] + B_deg1_total - B_deg1[i]);
    bound += 1;

    bad_prime_count = 0;

    n_poly_fit_length(alphapow, FLINT_MAX(WORD(3), bound + 1));
    n_poly_one(modulus);

    if ((mod.n & UWORD(1)) == UWORD(0))
    {
        success = -1;
        goto cleanup;
    }

    alpha = (mod.n - UWORD(1))/UWORD(2);

    goto choose_prime;

bad_prime:

    if (bad_prime_count > 2*bound)
    {
        success = n_poly_degree(modulus) > 0 ? -1 : -2;
        goto cleanup;
    }

    bad_prime_count++;

choose_prime:

    if (alpha < 2)
    {
        success = -1;
        goto cleanup;
    }

    alpha--;

    FLINT_ASSERT(0 < alpha && alpha <= mod.n/2);
    FLINT_ASSERT(alphapow->alloc >= 2);
    alphapow->length = 2;
    alphapow->coeffs[0] = 1;
    alphapow->coeffs[1] = alpha;

    n_bpoly_mod_interp_reduce_2sm_poly(Aevalp, Aevalm, A, alphapow, mod);
    for (i = 0; i < r; i++)
    {
        n_bpoly_mod_interp_reduce_2sm_poly(Bevalp + i, Bevalm + i, B + i, alphapow, mod);
        if (n_poly_degree(Bevalp + i) < n_bpoly_degree0(B + i) ||
            n_poly_degree(Bevalm + i) < n_bpoly_degree0(B + i))
        {
            goto choose_prime;
        }
    }

    for (i = 0; i < r; i++)
    {
        n_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            _n_poly_mod_mul(t1, t2, Bevalp + j, mod);
            n_poly_swap(t1, t2);
        }
        if (!n_poly_mod_invmod(t1, t2, Bevalp + i, mod))
            goto bad_prime;
        _n_poly_mod_mul(t2, Aevalp, t1, mod);
        _n_poly_mod_rem(Cevalp + i, t2, Bevalp + i, mod);

        n_poly_one(t2);
        for (j = 0; j < r; j++)
        {
            if (j == i)
                continue;
            _n_poly_mod_mul(t1, t2, Bevalm + j, mod);
            n_poly_swap(t1, t2);
        }
        if (!n_poly_mod_invmod(t1, t2, Bevalm + i, mod))
            goto bad_prime;
        _n_poly_mod_mul(t2, Aevalm, t1, mod);
        _n_poly_mod_rem(Cevalm + i, t2, Bevalm + i, mod);
    }

    if (n_poly_degree(modulus) > 0)
    {
        c = n_poly_mod_evaluate_nmod(modulus, alpha, mod);
        FLINT_ASSERT(c == n_poly_mod_evaluate_nmod(modulus, mod.n - alpha, mod));
        c = nmod_mul(c, alpha, mod);
        c = nmod_add(c, c, mod);
        c = nmod_inv(c, mod);
        _n_poly_mod_scalar_mul_nmod(modulus, modulus, c, mod);

        for (i = 0; i < r; i++)
        {
            n_bpoly_mod_interp_crt_2sm_poly(C_deg1 + i, C + i, T, Cevalp + i,
                                           Cevalm + i, modulus, alphapow, mod);
        }
    }
    else
    {
        for (i = 0; i < r; i++)
        {
            n_bpoly_mod_interp_lift_2sm_poly(C_deg1 + i, C + i, Cevalp + i,
                                                       Cevalm + i, alpha, mod);
        }
    }

    c = mod.n - nmod_mul(alpha, alpha, mod);
    n_poly_mod_shift_left_scalar_addmul(modulus, 2, c, mod);

    for (i = 0; i < r; i++)
    {
        if (C_deg1[i] > C_deg1_bound[i])
        {
            success = 0;
            goto cleanup;
        }
    }

    if (n_poly_degree(modulus) < bound)
    {
        goto choose_prime;
    }

    success = 1;

cleanup:

#if FLINT_WANT_ASSERT
    if (success == 1)
    {
        n_bpoly_t T1, T2, T3;
        n_bpoly_init(T1);
        n_bpoly_init(T2);
        n_bpoly_init(T3);
        n_bpoly_set(T1, A);
        for (i = 0; i < j; i++)
        {
            n_bpoly_one(T2);
            for (j = 0; j < r; j++)
            {
                n_bpoly_mod_mul(T3, T2, j == i ? C + j : B + j, mod);
                n_bpoly_swap(T2, T3);
            }
            n_bpoly_mod_sub(T3, T1, T2, mod);
            n_bpoly_swap(T1, T3);
        }
        FLINT_ASSERT(n_bpoly_is_zero(T1));
        n_bpoly_clear(T1);
        n_bpoly_clear(T2);
        n_bpoly_clear(T3);
    }
#endif

    n_poly_clear(Aevalp);
    n_poly_clear(Aevalm);
    n_poly_clear(modulus);
    n_poly_clear(alphapow);
    n_poly_clear(t1);
    n_poly_clear(t2);
    for (i = 0; i < r; i++)
    {
        n_poly_clear(Bevalp + i);
        n_poly_clear(Bevalm + i);
        n_poly_clear(Cevalp + i);
        n_poly_clear(Cevalm + i);
    }

    n_bpoly_clear(T);

    TMP_END;

    return success;
}
