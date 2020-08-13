/*
    Copyright (C) 2019 Daniel Schultz

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_mod.h"

void _fmpz_mod_mul1(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = nmod_mul(b0, c0, ctx->mod);
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

/*
    Multiplication modulo 2^FLINT_BITS is easy.
*/
void _fmpz_mod_mul2s(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a0, b0, c0;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    b0 = fmpz_get_ui(b);
    c0 = fmpz_get_ui(c);
    a0 = b0*c0;
    fmpz_set_ui(a, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}

/*
    Standard Barrett reduction: (set r = FLINT_BITS)

    We have n fits into 2 words and 2^r < n < 2^(2r). Therefore
    2^(3r) > 2^(4r) / n > 2^(2r) and the precomputed number
    ninv = floor(2^(4r) / n) fits into 3 words.
    The inputs b and c are < n and therefore fit into 2 words.

    The computation of a = b*c mod n is:

        x = b*c             x < n^2 and therefore fits into 4 words
        z = (x >> r)*ninv   z <= n*2^(3*r) and therefore fits into 5 words
        q = (z >> (3r))*n   q fits into 4 words
        x = x - q           x fits into 3 words after the subtraction
        at this point the canonical reduction in the range [0, n) is one of
            a = x, a = x - n, or a = x - 2n
*/
void _fmpz_mod_mul2(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    mp_limb_t a1, a0, b1, b0, c1, c0;
    mp_limb_t x3, x2, x1, x0;
    mp_limb_t q2, q1, q0;
    mp_limb_t z4, z3, z2, z1, z0;
    mp_limb_t t4, t3, t2, t1;
    mp_limb_t s3, s2, s1;
    mp_limb_t u4, u3, u2, u1;
    mp_limb_t v4, v3, v2, v1;

    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_get_uiui(&b1, &b0, b);
    fmpz_get_uiui(&c1, &c0, c);

    /* x[3:0] = b[1:0]*c[1:0] */
    umul_ppmm(t2, t1, b0, c1);
    umul_ppmm(s2, s1, b1, c0);
    umul_ppmm(x3, x2, b1, c1);
    umul_ppmm(x1, x0, b0, c0);
    t3 = 0;
    add_sssaaaaaa(t3, t2, t1, t3, t2, t1, 0, s2, s1);
    add_sssaaaaaa(x3, x2, x1, x3, x2, x1, t3, t2, t1);

    /* z[5:0] = x[3:1] * ninv[2:0], z[5] should end up zero */
    umul_ppmm(z1, z0, x1, ctx->ninv_limbs[0]);
    umul_ppmm(z3, z2, x2, ctx->ninv_limbs[1]);
    z4 = x3 * ctx->ninv_limbs[2];
    umul_ppmm(t3, t2, x3, ctx->ninv_limbs[0]);
    umul_ppmm(s3, s2, x1, ctx->ninv_limbs[2]);
    t4 = 0;
    add_sssaaaaaa(t4, t3, t2, t4, t3, t2, 0, s3, s2);
    umul_ppmm(u2, u1, x2, ctx->ninv_limbs[0]);
    umul_ppmm(u4, u3, x3, ctx->ninv_limbs[1]);
    add_sssaaaaaa(z4, z3, z2, z4, z3, z2, t4, t3, t2);
    umul_ppmm(v2, v1, x1, ctx->ninv_limbs[1]);
    umul_ppmm(v4, v3, x2, ctx->ninv_limbs[2]);
    add_ssssaaaaaaaa(z4, z3, z2, z1, z4, z3, z2, z1, u4, u3, u2, u1);
    add_ssssaaaaaaaa(z4, z3, z2, z1, z4, z3, z2, z1, v4, v3, v2, v1);

    /* q[3:0] = z[4:3] * n[1:0], q[3] is not needed */
    /* x[3:0] -= q[3:0], x[3] should end up zero */
    umul_ppmm(t2, t1, z3, ctx->n_limbs[1]);
    umul_ppmm(s2, s1, z4, ctx->n_limbs[0]);
    umul_ppmm(q1, q0, z3, ctx->n_limbs[0]);
    sub_ddmmss(x2, x1, x2, x1, t2, t1);
    q2 = z4 * ctx->n_limbs[1];
    sub_ddmmss(x2, x1, x2, x1, s2, s1);
    sub_dddmmmsss(x2, x1, x0, x2, x1, x0, q2, q1, q0);

    /* at most two subtractions of n, use q as temp space */
    sub_dddmmmsss(q2, q1, q0, x2, x1, x0, 0, ctx->n_limbs[1], ctx->n_limbs[0]);
    if ((slong)(q2) >= 0)
    {
        sub_dddmmmsss(x2, x1, x0, q2, q1, q0, 0, ctx->n_limbs[1], ctx->n_limbs[0]);
        if ((slong)(x2) >= 0)
        {
            a1 = x1;
            a0 = x0;
        }
        else
        {
            a1 = q1;
            a0 = q0;
        }
    }
    else
    {
        a1 = x1;
        a0 = x0;
    }

    fmpz_set_uiui(a, a1, a0);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));

}

void _fmpz_mod_mulN(fmpz_t a, const fmpz_t b, const fmpz_t c,
                                                     const fmpz_mod_ctx_t ctx)
{
    FLINT_ASSERT(fmpz_mod_is_canonical(b, ctx));
    FLINT_ASSERT(fmpz_mod_is_canonical(c, ctx));

    fmpz_mul(a, b, c);
    fmpz_mod(a, a, ctx->n);

    FLINT_ASSERT(fmpz_mod_is_canonical(a, ctx));
}
