/*
    Copyright (C) 2024 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

#include "gr_poly.h"

#include "mpn_mod.h"

/* for wide add_ssss.... macros. todo; these ought to be provided
   everywhere */
#if FLINT_BITS == 64 && defined(__AVX2__)
#include "crt_helpers.h"
#endif

static void
mpn_mod_set_mpn2(mp_ptr res, mp_srcptr s, mp_size_t l, gr_ctx_t ctx)
{
    MPN_NORM(s, l);
    mpn_mod_set_mpn(res, s, l, ctx);
}

#define FLINT_MPN_MUL_3_2X2(R2, R1, R0, a1, a0, b1, b0) \
    do \
    { \
        mp_limb_t __tmp2, __tmp1; \
        umul_ppmm(R1, R0, a0, b0); \
        (R2) = (a1) * (b1); \
        umul_ppmm(__tmp2, __tmp1, a0, b1); \
        add_ssaaaa(R2, R1, R2, R1, __tmp2, __tmp1); \
        umul_ppmm(__tmp2, __tmp1, a1, b0); \
        add_ssaaaa(R2, R1, R2, R1, __tmp2, __tmp1); \
    } \
    while (0) \

static int _mpn_mod_poly_divrem_q0_preinv1(mp_ptr Q, mp_ptr R,
                          mp_srcptr A, mp_srcptr B, slong lenA, mp_srcptr invL, gr_ctx_t ctx)
{
    /* special case for lenA == 1 omitted since this is dealt with
       in the calling function */

    mp_size_t nlimbs = MPN_MOD_CTX_NLIMBS(ctx);

    int monic = mpn_mod_is_one(invL, ctx) == T_TRUE;

    if (monic)
        mpn_mod_set(Q, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(Q, A + (lenA - 1) * nlimbs, invL, ctx);

    if (R == A)
    {
        mp_limb_t t[MPN_MOD_MAX_LIMBS];

        mpn_mod_neg(t, Q, ctx);
        _mpn_mod_vec_addmul_scalar(R, B, lenA - 1, t, ctx);
    }
    else
    {
        _mpn_mod_vec_mul_scalar(R, B, lenA - 1, Q, ctx);
        _mpn_mod_vec_sub(R, A, R, lenA - 1, ctx);
    }

    return GR_SUCCESS;
}

static int _mpn_mod_poly_divrem_q1_preinv1(mp_ptr Q, mp_ptr R,
                          mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                          mp_srcptr invL, gr_ctx_t ctx)
{
    mp_limb_t q0[MPN_MOD_MAX_LIMBS];
    mp_limb_t q1[MPN_MOD_MAX_LIMBS];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS + 1];
    mp_limb_t u[2 * MPN_MOD_MAX_LIMBS];
    slong i;
    mp_size_t nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    int monic = mpn_mod_is_one(invL, ctx) == T_TRUE;

    /* special case for lenB == 1 omitted since this is dealt with
       in the calling function */

    if (monic)
        mpn_mod_set(q1, A + (lenA - 1) * nlimbs, ctx);
    else
        mpn_mod_mul(q1, A + (lenA - 1) * nlimbs, invL, ctx);

    mpn_mod_mul(t, q1, B + (lenB - 2) * nlimbs, ctx);
    mpn_mod_sub(t, t, A + (lenA - 2) * nlimbs, ctx);

    if (monic)
        mpn_mod_set(q0, t, ctx);
    else
        mpn_mod_mul(q0, t, invL, ctx);

    mpn_mod_mul(t, q0, B, ctx);
    mpn_mod_add(R, A, t, ctx);

    mpn_mod_neg(Q, q0, ctx);
    mpn_mod_set(Q + nlimbs, q1, ctx);
    mpn_mod_neg(q1, q1, ctx);

#if defined(add_sssssaaaaaaaaaa)
    if (nlimbs == 2)
    {
        slong bits = 2 * MPN_MOD_CTX_MODULUS_BITS(ctx) + 1;
        slong slimbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

        if (slimbs == 3)
        {
            for (i = 1; i < lenB - 1; i++)
            {
                mp_srcptr B1ptr = B + (i - 1) * nlimbs;
                mp_srcptr Bptr = B + i * nlimbs;
                mp_srcptr Aptr = A + i * nlimbs;

                FLINT_MPN_MUL_3_2X2(t[2], t[1], t[0], q1[1], q1[0], B1ptr[1], B1ptr[0]);
                add_sssaaaaaa(t[2], t[1], t[0], t[2], t[1], t[0], 0, Aptr[1], Aptr[0]);
                FLINT_MPN_MUL_3_2X2(u[2], u[1], u[0], q0[1], q0[0], Bptr[1], Bptr[0]);
                add_sssaaaaaa(t[2], t[1], t[0], t[2], t[1], t[0], u[2], u[1], u[0]);
                mpn_mod_set_mpn2(R + i * nlimbs, t, slimbs, ctx);
            }
        }
        else
        {
            for (i = 1; i < lenB - 1; i++)
            {
                mp_srcptr B1ptr = B + (i - 1) * nlimbs;
                mp_srcptr Bptr = B + i * nlimbs;
                mp_srcptr Aptr = A + i * nlimbs;

                FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], q1[1], q1[0], B1ptr[1], B1ptr[0]);
                add_ssssaaaaaaaa(t[3], t[2], t[1], t[0], t[3], t[2], t[1], t[0], 0, 0, Aptr[1], Aptr[0]);
                FLINT_MPN_MUL_2X2(u[3], u[2], u[1], u[0], q0[1], q0[0], Bptr[1], Bptr[0]);
                add_sssssaaaaaaaaaa(t[4], t[3], t[2], t[1], t[0], 0, t[3], t[2], t[1], t[0], 0, u[3], u[2], u[1], u[0]);
                mpn_mod_set_mpn2(R + i * nlimbs, t, slimbs, ctx);
            }
        }
    }
    else
#endif
    {
        for (i = 1; i < lenB - 1; i++)
        {
            flint_mpn_mul_n(t, q1, B + (i - 1) * nlimbs, nlimbs);
            flint_mpn_mul_n(u, q0, B + i * nlimbs, nlimbs);
            t[2 * nlimbs] = mpn_add_n(t, t, u, 2 * nlimbs);
            mp_limb_t cy = mpn_add_n(t, t, A + i * nlimbs, nlimbs);
            mpn_add_1(t + nlimbs, t + nlimbs, nlimbs + 1, cy);
            mpn_mod_set_mpn2(R + i * nlimbs, t, 2 * nlimbs + 1, ctx);
        }
    }

    return GR_SUCCESS;
}

int
_mpn_mod_poly_divrem_basecase_preinv1(mp_ptr Q, mp_ptr R,
                                     mp_srcptr A, slong lenA, mp_srcptr B, slong lenB,
                                     mp_srcptr invL,
                                     gr_ctx_t ctx)
{
    slong iR, i, j;
    slong bits;
    mp_limb_t r[MPN_MOD_MAX_LIMBS];
    mp_limb_t c[MPN_MOD_MAX_LIMBS];
    mp_limb_t t[2 * MPN_MOD_MAX_LIMBS];
    mp_size_t slimbs, nlimbs;
    mp_ptr W;
    int monic;
    TMP_INIT;

    if (lenB == 1)
    {
        monic = (mpn_mod_is_one(invL, ctx) == T_TRUE);

        if (monic)
            return _mpn_mod_vec_set(Q, A, lenA, ctx);
        else
            return _mpn_mod_vec_mul_scalar(Q, A, lenA, invL, ctx);
    }
    else if (lenA == lenB + 1)
    {
        return _mpn_mod_poly_divrem_q1_preinv1(Q, R, A, lenA, B, lenB, invL, ctx);
    }
    else if (lenA == lenB)
    {
        return _mpn_mod_poly_divrem_q0_preinv1(Q, R, A, B, lenB, invL, ctx);
    }

    monic = (mpn_mod_is_one(invL, ctx) == T_TRUE);

    /* todo: case distinctions */
    nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    bits = 2 * MPN_MOD_CTX_MODULUS_BITS(ctx) + FLINT_BIT_COUNT(lenA - lenB + 1);
    slimbs = (bits + FLINT_BITS - 1) / FLINT_BITS;

    FLINT_ASSERT((slimbs == 2 * nlimbs) || (slimbs == 2 * nlimbs + 1) || (slimbs == 2 * nlimbs - 1));

    TMP_START;
    W = TMP_ALLOC(lenA * slimbs * sizeof(mp_limb_t));

    mp_ptr R3 = W;
    mp_ptr ptrQ = Q - nlimbs * (lenB - 1);

    for (i = 0; i < lenA; i++)
    {
        flint_mpn_copyi(R3 + slimbs * i, A + nlimbs * i, nlimbs);
        flint_mpn_zero(R3 + slimbs * i + nlimbs, slimbs - nlimbs);
    }

    for (iR = lenA - 1; iR >= lenB - 1; )
    {
        mpn_mod_set_mpn2(r, R3 + slimbs * iR, slimbs, ctx);

        while ((iR + 1 >= lenB) && flint_mpn_zero_p(r, nlimbs))
        {
            flint_mpn_zero(ptrQ + iR * nlimbs, nlimbs);
            iR--;

            if (iR + 1 >= lenB)
                mpn_mod_set_mpn2(r, R3 + slimbs * iR, slimbs, ctx);
        }

        if (iR + 1 >= lenB)
        {
            /* note: unoptimized in nmod_poly */
            if (monic)
                flint_mpn_copyi(ptrQ + iR * nlimbs, r, nlimbs);
            else
                mpn_mod_mul(ptrQ + iR * nlimbs, r, invL, ctx);

            if (lenB > 1)
            {
                mpn_mod_neg(c, ptrQ + iR * nlimbs, ctx);

                /* todo: consider writing all products to a temporary
                   buffer and doing a single big mpn_add_n */

#if defined(add_sssssaaaaaaaaaa)
                if (nlimbs == 2)
                {
                    mp_limb_t t[4];

                    if (slimbs == 5)
                    {
                        for (j = 0; j < lenB - 1; j++)
                        {
                            mp_ptr a = R3 + slimbs * (iR - lenB + 1 + j);
                            mp_srcptr b = B + j * nlimbs;
                            FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], b[1], b[0], c[1], c[0]);
                            add_sssssaaaaaaaaaa(a[4], a[3], a[2], a[1], a[0], a[4], a[3], a[2], a[1], a[0], 0, t[3], t[2], t[1], t[0]);
                        }
                    }
                    else if (slimbs == 4)
                    {
                        for (j = 0; j < lenB - 1; j++)
                        {
                            mp_ptr a = R3 + slimbs * (iR - lenB + 1 + j);
                            mp_srcptr b = B + j * nlimbs;
                            FLINT_MPN_MUL_2X2(t[3], t[2], t[1], t[0], b[1], b[0], c[1], c[0]);
                            add_ssssaaaaaaaa(a[3], a[2], a[1], a[0], a[3], a[2], a[1], a[0], t[3], t[2], t[1], t[0]);
                        }
                    }
                    else
                    {
                        for (j = 0; j < lenB - 1; j++)
                        {
                            mp_ptr a = R3 + slimbs * (iR - lenB + 1 + j);
                            mp_srcptr b = B + j * nlimbs;
                            FLINT_MPN_MUL_3_2X2(t[2], t[1], t[0], b[1], b[0], c[1], c[0]);
                            add_sssaaaaaa(a[2], a[1], a[0], a[2], a[1], a[0], t[2], t[1], t[0]);
                        }
                    }
                }
                else
#endif
                {
                    if (slimbs == 2 * nlimbs + 1)
                    {
                        for (j = 0; j < lenB - 1; j++)
                        {
                            flint_mpn_mul_n(t, B + j * nlimbs, c, nlimbs);
                            R3[slimbs * (iR - lenB + 1 + j) + slimbs - 1] += mpn_add_n(R3 + slimbs * (iR - lenB + 1 + j), R3 + slimbs * (iR - lenB + 1 + j), t, 2 * nlimbs);
                        }
                    }
                    else  /* slimbs == 2 * nlimbs or 2 * nlimbs - 1 */
                    {
                        for (j = 0; j < lenB - 1; j++)
                        {
                            flint_mpn_mul_n(t, B + j * nlimbs, c, nlimbs);
                            mpn_add_n(R3 + slimbs * (iR - lenB + 1 + j), R3 + slimbs * (iR - lenB + 1 + j), t, slimbs);
                        }
                    }
                }
            }

            iR--;
        }
    }

    for (iR = 0; iR < lenB - 1; iR++)
        mpn_mod_set_mpn2(R + iR * nlimbs, R3 + slimbs * iR, slimbs, ctx);

    TMP_END;

    return GR_SUCCESS;
}

int
_mpn_mod_poly_divrem_basecase(mp_ptr Q, mp_ptr R, mp_srcptr A, slong lenA,
                                  mp_srcptr B, slong lenB, gr_ctx_t ctx)
{
    mp_limb_t invB[MPN_MOD_MAX_LIMBS];
    int status;

    status = mpn_mod_inv(invB, B + (lenB - 1) * MPN_MOD_CTX_NLIMBS(ctx), ctx);
    if (status != GR_SUCCESS)
        return status;

    return _mpn_mod_poly_divrem_basecase_preinv1(Q, R, A, lenA, B, lenB, invB, ctx);
}
