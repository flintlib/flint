/*
    Copyright (C) 2021 Fredrik Johansson
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "ulong_extras.h"
#include "mpn_extras.h"
#include "fmpz.h"
#include "fmpz_vec.h"

/*
    LCM of many integers, typically denominators to be cleared.

    The running lcm L only grows, and in the common cases (all
    denominators 1 or equal, or dividing a common denominator) each entry d
    divides it. So each entry is first tested for divisibility, and L is
    only updated, by L <- L (d / gcd(L, d)), when the test fails.

    The small entries are done first, in a single limb as long as L fits:
    an entry equal to the previous one is skipped, and d | L is tested
    with the Hensel quotient (L (d/2^k)^{-1} mod B, see gcd_vec.c) without
    a division. When L outgrows a limb, the small entries are packed into
    products M of a limb, so that one pass r = L mod M over L tests all of
    them (r mod d = 0). The large entries are then done in a separate
    multi-limb lcm, where the divisibility test goes in the direction
    suggested by the sizes (d | L for a shorter d, L | d for a longer d,
    which then replaces L, as for growing denominators like factorials),
    and the two lcms are combined at the end.
*/

/* multi-limb lcm {p, n}, n = 0 for 1, with room for alloc limbs */
typedef struct
{
    mp_ptr p;
    mp_size_t n;
    mp_size_t alloc;
    slong fails;    /* consecutive failed divisibility tests */
}
lcm_big_t;

/* After this many consecutive failed divisibility tests of multi-limb
   entries (as for random denominators, where every entry needs a gcd),
   only one entry in LCM_RETRY is tested. */
#define LCM_FAILS 2
#define LCM_RETRY 8

static void
_lcm_big_init(lcm_big_t * L)
{
    L->p = NULL;
    L->n = 0;
    L->alloc = 0;
    L->fails = 0;
}

static void
_lcm_big_clear(lcm_big_t * L)
{
    flint_free(L->p);
}

static void
_lcm_big_fit(lcm_big_t * L, mp_size_t n)
{
    if (n > L->alloc)
    {
        mp_size_t alloc = FLINT_MAX(n, 2 * L->alloc);
        L->p = flint_realloc(L->p, alloc * sizeof(mp_limb_t));
        L->alloc = alloc;
    }
}

/* L <- L m for a limb m */
static void
_lcm_big_mul_1(lcm_big_t * L, ulong m)
{
    mp_limb_t cy;

    _lcm_big_fit(L, L->n + 1);
    cy = mpn_mul_1(L->p, L->p, L->n, m);
    if (cy != 0)
        L->p[L->n++] = cy;
}

/* L <- lcm(L, d) for a single limb d >= 2, L >= 2 */
static void
_lcm_big_lcm_1(lcm_big_t * L, ulong d)
{
    ulong r = mpn_mod_1(L->p, L->n, d);

    if (r != 0)
        _lcm_big_mul_1(L, d / n_gcd(d, r));
}

/* L <- lcm(L, d) for a normalised d of dn >= 1 limbs, d > 1 */
static void
_lcm_big_lcm(lcm_big_t * L, mp_srcptr d, mp_size_t dn)
{
    mp_size_t Ln = L->n, gn, qn;
    mp_ptr g, q, t;
    TMP_INIT;

    if (Ln == 0)
    {
        _lcm_big_fit(L, dn);
        flint_mpn_copyi(L->p, d, dn);
        L->n = dn;
        return;
    }

    if (dn == 1)
    {
        if (Ln == 1)
        {
            /* both single limbs */
            ulong a = L->p[0], b = d[0], hi, lo;
            if (a % b == 0)
                return;
            b = b / n_gcd(a, b);
            umul_ppmm(hi, lo, a, b);
            _lcm_big_fit(L, 2);
            L->p[0] = lo;
            L->p[1] = hi;
            L->n = 1 + (hi != 0);
            return;
        }

        _lcm_big_lcm_1(L, d[0]);
        return;
    }

    if (dn == Ln && flint_mpn_equal_p(L->p, d, dn))
        return;

    if (L->fails >= LCM_FAILS && (L->fails++ - LCM_FAILS) % LCM_RETRY != 0)
        goto general;

    /* the divisibility test in the direction suggested by the sizes */
    if (dn > Ln)
    {
        if (flint_mpn_divisible(d, dn, L->p, Ln))
        {
            _lcm_big_fit(L, dn);
            flint_mpn_copyi(L->p, d, dn);
            L->n = dn;
            L->fails = 0;
            return;
        }
    }
    else if (dn == Ln)
    {
        if (mpn_cmp(L->p, d, dn) > 0)
        {
            if (flint_mpn_divisible(L->p, Ln, d, dn))
            {
                L->fails = 0;
                return;
            }
        }
        else if (flint_mpn_divisible(d, dn, L->p, Ln))
        {
            flint_mpn_copyi(L->p, d, dn);
            L->fails = 0;
            return;
        }
    }
    else
    {
        if (flint_mpn_divisible(L->p, Ln, d, dn))
        {
            L->fails = 0;
            return;
        }
    }

    if (L->fails < LCM_FAILS)
        L->fails++;

general:
    /* L <- L (d / gcd(L, d)) */
    TMP_START;
    g = TMP_ALLOC((FLINT_MIN(Ln, dn) + dn + Ln + dn) * sizeof(mp_limb_t));
    q = g + FLINT_MIN(Ln, dn);
    t = q + dn;

    gn = flint_mpn_gcd(g, L->p, Ln, d, dn);

    if (gn == 1 && g[0] == 1)
    {
        flint_mpn_copyi(q, d, dn);
        qn = dn;
    }
    else
    {
        qn = dn - gn + 1;
        flint_mpn_divexact(q, d, dn, g, gn);
        MPN_NORM(q, qn);
    }

    if (Ln >= qn)
        flint_mpn_mul(t, L->p, Ln, q, qn);
    else
        flint_mpn_mul(t, q, qn, L->p, Ln);
    Ln = Ln + qn;
    Ln -= (t[Ln - 1] == 0);

    _lcm_big_fit(L, Ln);
    flint_mpn_copyi(L->p, t, Ln);
    L->n = Ln;

    TMP_END;
}

/* the small entries packed into a product of a limb, tested against a
   multi-limb L by one remainder */
#define LCM_BATCH 64

typedef struct
{
    ulong d[LCM_BATCH];
    slong n;
    ulong M;
}
lcm_batch_t;

static void
_lcm_batch_flush(lcm_big_t * L, lcm_batch_t * B)
{
    ulong r, M = B->M, Minv, dj, rj, m, hi, lo;
    slong j;

    if (B->n == 0)
        return;

    r = mpn_mod_1(L->p, L->n, M);
    Minv = n_preinvert_limb(M);

    for (j = 0; j < B->n; j++)
    {
        dj = B->d[j];
        rj = r % dj;

        if (rj != 0)
        {
            /* L <- L m, and r = L mod M follows */
            m = dj / n_gcd(dj, rj);
            _lcm_big_mul_1(L, m);
            umul_ppmm(hi, lo, r, m);
            r = n_ll_mod_preinv(hi, lo, M, Minv);
        }
    }

    B->n = 0;
    B->M = 1;
}

/*
    res = lcm(inp, vec[0], vec[stride], ..., vec[(len - 1) stride]), with
    inp = NULL for 1. Returns 0 (with res set to 0) if some entry is zero.
*/
static void
_fmpz_vec_lcm_inp(fmpz_t res, const fmpz * vec, slong stride, slong len, const fmpz * inp)
{
    ulong L1, lastd, d, k, dd, hi, lo, r, m;
    unsigned int vL;
    slong i, nbig, fails = 0;
    lcm_big_t Ls, Lb;
    lcm_batch_t B;
    fmpz c;
    int multi;

    _lcm_big_init(&Ls);
    _lcm_big_init(&Lb);

    L1 = 1;
    nbig = 0;
    if (inp != NULL)
    {
        if (fmpz_is_zero(inp))
            goto zero;
        if (!COEFF_IS_MPZ(*inp))
            L1 = FLINT_ABS(*inp);
        else
        {
            mpz_ptr z = COEFF_TO_PTR(*inp);
            _lcm_big_lcm(&Lb, z->_mp_d, FLINT_ABS(z->_mp_size));
        }
    }

    /* the small entries, in a single limb while L fits */
    vL = flint_ctz(L1);
    lastd = 1;
    multi = 0;
    for (i = 0; i < len; i++)
    {
        c = vec[i * stride];

        if (COEFF_IS_MPZ(c))
        {
            nbig++;
            continue;
        }

        if (c == 0)
            goto zero;

        d = FLINT_ABS(c);
        if (d == lastd)
            continue;
        lastd = d;

        /* d = 2^k dd divides L iff k <= v(L) and the Hensel quotient
           L dd^{-1} mod B times dd gives L, i.e. its high part is 0; this
           avoids a division, which is only needed if the test fails, so
           after repeated failures the division is done directly */
        if (fails < LCM_FAILS)
        {
            k = flint_ctz(d);
            if (k <= vL)
            {
                dd = d >> k;
                if (n_mulhi(L1 * n_binvert(dd), dd) == 0)
                {
                    fails = 0;
                    continue;
                }
            }
            fails++;
            r = L1 % d;
        }
        else
        {
            r = L1 % d;
            if (r == 0)
            {
                fails = 0;
                continue;
            }
        }

        m = d / n_gcd(d, r);
        umul_ppmm(hi, lo, L1, m);
        if (hi == 0)
        {
            L1 = lo;
            vL = flint_ctz(L1);
        }
        else
        {
            _lcm_big_fit(&Ls, 2);
            Ls.p[0] = lo;
            Ls.p[1] = hi;
            Ls.n = 2;
            multi = 1;
            i++;
            break;
        }
    }

    /* L no longer fits in a limb: packed products of the entries */
    if (multi)
    {
        B.n = 0;
        B.M = 1;

        for ( ; i < len; i++)
        {
            c = vec[i * stride];

            if (COEFF_IS_MPZ(c))
            {
                nbig++;
                continue;
            }

            if (c == 0)
                goto zero;

            d = FLINT_ABS(c);
            if (d == lastd || d == 1)
                continue;
            lastd = d;

            umul_ppmm(hi, lo, B.M, d);
            if (hi != 0 || B.n == LCM_BATCH)
            {
                _lcm_batch_flush(&Ls, &B);
                lo = d;
            }

            B.d[B.n++] = d;
            B.M = lo;
        }

        _lcm_batch_flush(&Ls, &B);
    }

    /* the large entries */
    if (nbig != 0)
    {
        for (i = 0; i < len; i++)
        {
            c = vec[i * stride];
            if (COEFF_IS_MPZ(c))
            {
                mpz_ptr z = COEFF_TO_PTR(c);
                _lcm_big_lcm(&Lb, z->_mp_d, FLINT_ABS(z->_mp_size));
            }
        }
    }

    /* combine */
    if (Lb.n == 0)
    {
        if (!multi)
        {
            fmpz_set_ui(res, L1);
        }
        else
        {
            mpz_ptr z = _fmpz_promote(res);
            flint_mpn_copyi(FLINT_MPZ_REALLOC(z, Ls.n), Ls.p, Ls.n);
            z->_mp_size = Ls.n;
            _fmpz_demote_val(res);
        }
    }
    else
    {
        mpz_ptr z;

        if (!multi)
        {
            if (L1 != 1)
            {
                if (Lb.n == 1)
                    _lcm_big_lcm(&Lb, &L1, 1);
                else
                    _lcm_big_lcm_1(&Lb, L1);
            }
        }
        else
        {
            _lcm_big_lcm(&Lb, Ls.p, Ls.n);
        }

        z = _fmpz_promote(res);
        flint_mpn_copyi(FLINT_MPZ_REALLOC(z, Lb.n), Lb.p, Lb.n);
        z->_mp_size = Lb.n;
        _fmpz_demote_val(res);
    }

    _lcm_big_clear(&Ls);
    _lcm_big_clear(&Lb);
    return;

zero:
    _lcm_big_clear(&Ls);
    _lcm_big_clear(&Lb);
    fmpz_zero(res);
}

void
_fmpz_vec_lcm(fmpz_t res, const fmpz * vec, slong len)
{
    _fmpz_vec_lcm_inp(res, vec, 1, len, NULL);
}

void
_fmpz_vec_lcm_strided(fmpz_t res, const fmpz * vec, slong stride, slong len, const fmpz_t inp)
{
    _fmpz_vec_lcm_inp(res, vec, stride, len, inp);
}
