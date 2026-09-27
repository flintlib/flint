/*
    Copyright (C) 2026 Maël Hostettler

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

/*
    Curves over Z/NZ.

    When the discriminant is a unit modulo N, E is an elliptic curve over
    the ring Z/NZ, and its points (in the projective plane over Z/NZ) form
    a group that splits along the factorisation of N,

        E(Z/NZ) = prod over p^k || N of E(Z/p^k Z)

    by the Chinese remainder theorem. Over Z/p^k Z, reduction modulo p maps
    onto E(F_p) -- every point there is smooth and lifts by Hensel's lemma
    -- and its kernel, the points congruent to the identity, has p^(k-1)
    elements (it is the formal group evaluated at pZ/p^kZ). So

        #E(Z/NZ) = prod over p^k || N of p^(k-1) #E(F_p),

    and each #E(F_p) is an ordinary point count over a prime field, done by
    gr_ec_ctx_cardinality with everything it knows.

    When the discriminant is not a unit, E has bad reduction at some p | N
    and is not an elliptic curve over Z/NZ at all; that is GR_DOMAIN.

    What this costs is factoring N. A prime power is always recognised (a
    perfect power test and a primality proof); anything else is factored
    in full up to GR_EC_CRT_FACTOR_BITS bits, and above that only when
    splitting off the small primes leaves a prime power.
*/

#include "fmpz.h"
#include "fmpz_factor.h"
#include "gr.h"
#include "gr_ec.h"
#include "impl.h"

#define GR_EC_CRT_FACTOR_BITS 128

/* is the base ring one of the Z/NZ types? */
static int
_is_zn(gr_ctx_t R)
{
    switch (R->which_ring)
    {
        case GR_CTX_NMOD:
        case GR_CTX_NMOD8:
        case GR_CTX_NMOD32:
        case GR_CTX_FMPZ_MOD:
        case GR_CTX_MPN_MOD:
            return 1;
        default:
            return 0;
    }
}

/* n = p^e with p prime, proved; *e = 0 if not */
static void
_prime_power(fmpz_t p, ulong * e, const fmpz_t n)
{
    fmpz_t r;
    int k;

    fmpz_init(r);
    fmpz_set(p, n);
    *e = 1;

    while (fmpz_cmp_ui(p, 3) > 0 && (k = fmpz_is_perfect_power(r, p)) > 1)
    {
        fmpz_swap(p, r);
        *e *= k;
    }

    if (fmpz_cmp_ui(p, 1) <= 0 || !fmpz_is_prime(p))
        *e = 0;

    fmpz_clear(r);
}

/*
    The factorisation of n into proved primes, or 0 when that is not
    affordable (see above) or a factor could not be proved prime.
*/
static int
_factor(fmpz_factor_t fac, const fmpz_t n)
{
    fmpz_t p;
    ulong e;
    slong i;
    int ok = 1;

    fmpz_init(p);

    /* the common case: n = p^k, with no search at all */
    _prime_power(p, &e, n);

    if (e != 0)
    {
        _fmpz_factor_append(fac, p, e);
        fmpz_clear(p);
        return 1;
    }

    if (fmpz_bits(n) <= GR_EC_CRT_FACTOR_BITS)
        fmpz_factor(fac, n);
    else if (!fmpz_factor_smooth(fac, n, 32, 1))
    {
        /* what is left must be a prime power */
        slong last = fac->num - 1;
        ulong ce = fac->exp[last];

        _prime_power(p, &e, fac->p + last);

        if (e == 0)
            ok = 0;
        else
        {
            fmpz_swap(fac->p + last, p);
            fac->exp[last] = ce * e;
        }
    }

    for (i = 0; i < fac->num && ok; i++)
        if (!fmpz_is_prime(fac->p + i))
            ok = 0;

    fmpz_clear(p);

    return ok;
}

/* #E(F_p) for E reduced modulo p */
static int
_count_mod_p(fmpz_t res, const fmpz_t p, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    gr_ctx_t F;
    gr_ec_ctx_t Ep;
    gr_ptr a;
    fmpz_t c;
    slong i;
    int status = GR_SUCCESS;

    if (fmpz_abs_fits_ui(p))
    {
        if (gr_ctx_init_nmod(F, fmpz_get_ui(p)) != GR_SUCCESS)
            return GR_UNABLE;
    }
    else
        gr_ctx_init_fmpz_mod(F, p);

    /* p is proved prime */
    status |= gr_ctx_set_is_field(F, T_TRUE);

    fmpz_init(c);
    a = gr_heap_init_vec(5, F);

    for (i = 0; i < 5; i++)
    {
        status |= gr_get_fmpz(c, GR_EC_COEFF(ctx, i), R);
        status |= gr_set_fmpz(GR_ENTRY(a, i, F->sizeof_elem), c, F);
    }

    if (status == GR_SUCCESS)
    {
        status = gr_ec_ctx_init(Ep, F, GR_ENTRY(a, 0, F->sizeof_elem),
                GR_ENTRY(a, 1, F->sizeof_elem), GR_ENTRY(a, 2, F->sizeof_elem),
                GR_ENTRY(a, 3, F->sizeof_elem), GR_ENTRY(a, 4, F->sizeof_elem));

        /* the discriminant is a unit modulo N, so the reduction is smooth */
        if (status == GR_SUCCESS)
        {
            status = gr_ec_ctx_cardinality(res, Ep);
            gr_ec_ctx_clear(Ep);
        }
        else
            status = GR_UNABLE;
    }

    gr_heap_clear_vec(a, 5, F);
    fmpz_clear(c);
    gr_ctx_clear(F);

    return status;
}

int
gr_ec_ctx_cardinality_crt(fmpz_t res, gr_ec_ctx_t ctx)
{
    gr_ctx_struct * R = GR_EC_ELEM_CTX(ctx);
    fmpz_factor_t fac;
    fmpz_t N, Np, t;
    slong i;
    int status = GR_SUCCESS;

    if (!_is_zn(R))
        return GR_DOMAIN;

    fmpz_init(N);
    fmpz_init(Np);
    fmpz_init(t);
    fmpz_factor_init(fac);

    if (gr_ctx_cardinality_fmpz(N, R) != GR_SUCCESS || fmpz_cmp_ui(N, 1) <= 0)
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    /*
        An elliptic curve over Z/NZ has a unit discriminant. Decided by a
        gcd rather than gr_is_invertible, which not every Z/NZ type can
        answer for a composite modulus.
    */
    if (gr_get_fmpz(t, GR_EC_DISC(ctx), R) != GR_SUCCESS)
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    fmpz_gcd(t, t, N);

    if (!fmpz_is_one(t))
    {
        status = GR_DOMAIN;
        goto cleanup;
    }

    if (!_factor(fac, N))
    {
        status = GR_UNABLE;
        goto cleanup;
    }

    /* prod p^(k-1) #E(F_p) */
    fmpz_one(res);

    for (i = 0; i < fac->num && status == GR_SUCCESS; i++)
    {
        status |= _count_mod_p(Np, fac->p + i, ctx);

        if (status == GR_SUCCESS)
        {
            fmpz_pow_ui(t, fac->p + i, fac->exp[i] - 1);
            fmpz_mul(t, t, Np);
            fmpz_mul(res, res, t);
        }
    }

cleanup:
    fmpz_factor_clear(fac);
    fmpz_clear(t);
    fmpz_clear(Np);
    fmpz_clear(N);

    return status;
}
