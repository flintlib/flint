/*
    Copyright (C) 2024 Fredrik Johansson
    Copyright (C) 2024 Éric Schost
    Copyright (C) 2024 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "mpn_mod.h"

/* compute c += (a1 + b1) * (a2 + b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_addadd(mp_ptr val0, mp_ptr val1, mp_ptr val2, mp_ptr c, mp_srcptr a1, mp_srcptr b1, mp_srcptr a2, mp_srcptr b2, mp_size_t nlimbs, int add_can_overflow_nlimbs)
{
    if (!add_can_overflow_nlimbs)
    {
        mpn_add_n(val1, a1, b1, nlimbs);
        mpn_add_n(val2, a2, b2, nlimbs);
        flint_mpn_mul_n(val0, val1, val2, nlimbs);
        c[2 * nlimbs] += mpn_add_n(c, c, val0, 2 * nlimbs);
    }
    else
    {
        val1[nlimbs] = mpn_add_n(val1, a1, b1, nlimbs);
        val2[nlimbs] = mpn_add_n(val2, a2, b2, nlimbs);
        flint_mpn_mul_n(val0, val1, val2, nlimbs + 1);
        /* we write 2 * nlimbs + 2 limbs to val0, but the top limb will be zero */
        FLINT_ASSERT(val0[2 * nlimbs + 1] == 0);
        mpn_add_n(c, c, val0, 2 * nlimbs + 1);
    }
}

/* compute c += (a1 - b1) * (a2 - b2) */
/* val0, val1, val2 are scratch space */
FLINT_FORCE_INLINE void
addmul_subsub(mp_ptr val0, mp_ptr val1, mp_ptr val2, mp_ptr c, mp_srcptr a1, mp_srcptr b1, mp_srcptr a2, mp_srcptr b2, mp_size_t nlimbs)
{
    int neg;
    neg = flint_mpn_signed_sub_n(val1, a1, b1, nlimbs);
    neg ^= flint_mpn_signed_sub_n(val2, a2, b2, nlimbs);

    flint_mpn_mul_n(val0, val1, val2, nlimbs);
    if (neg)
        c[2 * nlimbs] -= mpn_sub_n(c, c, val0, 2 * nlimbs);
    else
        c[2 * nlimbs] += mpn_add_n(c, c, val0, 2 * nlimbs);
}

int mpn_mod_mat_mul_waksman(gr_mat_t C, const gr_mat_t A, const gr_mat_t B, gr_ctx_t ctx)
{
    slong nlimbs = MPN_MOD_CTX_NLIMBS(ctx);
    /* Enough to hold any unreduced value, including one sign bit. TODO: tighten. */
    slong slimbs = 2 * nlimbs + 1;
    slong m = A->r;
    slong n = B->r;
    slong p = B->c;
    /* Normally the sum of two input entries fits in nlimbs, but we may need
       an extra limb for a carry bit. */
    int add_can_overflow_nlimbs = (MPN_MOD_CTX_NORM(ctx) == 0);

    if (m == 0 || n == 0 || p == 0)
        return gr_mat_zero(C, ctx);

    slong i, l, j, k;

    mp_ptr Ctmp = flint_calloc(slimbs * ((m * p) + (p + m) + 5), sizeof(mp_limb_t));
                                            /* Ctmp itself has m * p entries */
    mp_ptr Crow = Ctmp + slimbs * (m * p);  /* Crow has p entries */
    mp_ptr Ccol = Crow + slimbs * p;        /* Ccol has m entries */
    mp_ptr val0 = Ccol + slimbs * m;        /* val0 has room for 2 sums */
    mp_ptr val1 = val0 + 2 * slimbs;        /* val1 has room for 1 sum   */
    mp_ptr val2 = val1 + slimbs;            /* val2 has room for 1 sum   */
    mp_ptr crow = val2 + slimbs;            /* crow has room for 1 sum   */

#define A_ENTRY(ii, jj) (((mp_srcptr) A->rows[ii]) + (jj) * nlimbs)
#define B_ENTRY(ii, jj) (((mp_srcptr) B->rows[ii]) + (jj) * nlimbs)

#define C_ENTRY(ii, jj) (Ctmp + ((ii) * p + (jj)) * slimbs)
#define Crow_ENTRY(ii) (Crow + (ii) * slimbs)
#define Ccol_ENTRY(ii) (Ccol + (ii) * slimbs)

    slong np = n >> 1;

    for (j = 1; j <= np; j++)
    {
        slong j2 = (j << 1) - 1;

        for (k = 0; k < p; k++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(0, k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs, add_can_overflow_nlimbs);
            addmul_subsub(val0, val1, val2, Crow_ENTRY(k), A_ENTRY(0, j2-1), B_ENTRY(j2, k), A_ENTRY(0, j2), B_ENTRY(j2-1, k), nlimbs);
        }

        for (l = 1; l < m; l++)
        {
            addmul_addadd(val0, val1, val2, C_ENTRY(l, 0), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs, add_can_overflow_nlimbs);
            addmul_subsub(val0, val1, val2, Ccol_ENTRY(l), A_ENTRY(l, j2-1), B_ENTRY(j2, 0), A_ENTRY(l, j2), B_ENTRY(j2-1, 0), nlimbs);
        }

        for (k = 1; k < p; k++)
        {
            for (l = 1; l < m; l++)
            {
                addmul_addadd(val0, val1, val2, C_ENTRY(l, k), A_ENTRY(l, j2-1), B_ENTRY(j2, k), A_ENTRY(l, j2), B_ENTRY(j2-1, k), nlimbs, add_can_overflow_nlimbs);
            }
        }
    }

    for (l = 1; l < m; l++)
    {
        mpn_add_n(val1, Ccol_ENTRY(l), C_ENTRY(l, 0), slimbs);
        flint_mpn_signed_div2(Ccol_ENTRY(l), val1, slimbs);
        mpn_sub_n(C_ENTRY(l, 0), C_ENTRY(l, 0), Ccol_ENTRY(l), slimbs);
    }

    mpn_add_n(val1, Crow, C_ENTRY(0, 0), slimbs);
    flint_mpn_signed_div2(val0, val1, slimbs);
    mpn_sub_n(C_ENTRY(0, 0), C_ENTRY(0, 0), val0, slimbs);

    for (k = 1; k < p; k++)
    {
        mpn_add_n(crow, Crow_ENTRY(k), C_ENTRY(0, k), slimbs);
        flint_mpn_signed_div2(val1, crow, slimbs);
        mpn_sub_n(C_ENTRY(0, k), C_ENTRY(0, k), val1, slimbs);
        mpn_sub_n(crow, val1, val0, slimbs);

        for (l = 1; l < m; l++)
        {
            mpn_sub_n(val2, C_ENTRY(l, k), crow, slimbs);
            mpn_sub_n(C_ENTRY(l, k), val2, Ccol_ENTRY(l), slimbs);
        }
    }

    if ((n & 1) == 1)
    {
        for (l = 0; l < m; l++)
        {
            for (k = 0; k < p; k++)
            {
                flint_mpn_mul_n(val0, A_ENTRY(l, n-1), B_ENTRY(n-1, k), nlimbs);
                C_ENTRY(l, k)[2 * nlimbs] += mpn_add_n(C_ENTRY(l, k), C_ENTRY(l, k), val0, 2 * nlimbs);
            }
        }
    }

    /* Reduce and write output */
    for (i = 0; i < m; i++)
    {
        for (k = 0; k < p; k++)
        {
            mp_size_t d;
            mp_ptr Cptr = ((mp_ptr) C->rows[i]) + k * nlimbs;

            /* As currently implemented, there is no wraparound arithmetic.
               Were that the case, we would need something like

                d = slimbs;
                mpn_neg(C_ENTRY(i, k), C_ENTRY(i, k), d);
                MPN_NORM(C_ENTRY(i, k), d);
                mpn_mod_set_mpn(Cptr, C_ENTRY(i, k), d, ctx);
                mpn_mod_neg(Cptr, Cptr, ctx);

               when the leading limb indicates a negative value. */
            FLINT_ASSERT((slong) C_ENTRY(i, k)[slimbs - 1] >= 0);

            d = slimbs;
            MPN_NORM(C_ENTRY(i, k), d);
            mpn_mod_set_mpn(Cptr, C_ENTRY(i, k), d, ctx);
        }
    }

    flint_free(Ctmp);

    return GR_SUCCESS;
}
