/*
    Copyright (C) 2026 Vincent Neiger

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "test_helpers.h"
#include "ulong_extras.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly/impl.h"

/* the number of levels of the recursion on blocks of 16 for len
   coefficients */
static slong
_kt_depth16(slong len)
{
    slong d = 0, m = (len + 15) / 16;

    while (m > 1)
    {
        m = (m + 1) / 2;
        d++;
    }

    return d;
}

/* a modulus, often at the limits of the one-word code for blocks of len
   coefficients and terms products (n below 2^(32 - d), 2^(31 - d) for
   middle products, and terms (n - 1)^2 < 2^64), or a power of 2, or any */
static ulong
_kt_modulus(flint_rand_t state, slong len, slong terms, int mid)
{
    ulong n, b1, b2;
    slong d = _kt_depth16(len);

    switch (n_randint(state, 4))
    {
        case 0:
            /* the largest n the one-word code takes, or a little less */
#if FLINT_BITS == 64
            b1 = (UWORD(1) << (mid ? 31 : 32)) >> d;
            b2 = n_sqrt(UWORD_MAX / (ulong) terms) + 1;
            n = FLINT_MIN(b1, b2);
            if (n > 2 && n_randint(state, 2))
                n -= n_randint(state, FLINT_MIN(n - 1, 100));
#else
            (void) b1;
            (void) b2;
            (void) d;
            n = n_randtest_not_zero(state);
#endif
            break;
        case 1:
            n = UWORD(1) << n_randint(state, 31);
            break;
        case 2:
            n = n_randtest_bits(state, 1 + n_randint(state, 30));
            break;
        default:
            n = n_randtest_bits(state, 1 + n_randint(state, FLINT_BITS));
    }

    return FLINT_MAX(n, 1);
}

/* a modulus with at most maxbits bits, often one where a rule changes:
   powers of 2 (including 2^52), 2^k + 1 and 2^k - 1, the bound 1515531528
   of the doubles, a prime of fft_small, or a prime p with a large power
   of 2 dividing p - 1 */
static ulong
_kt_rule_modulus(flint_rand_t state, slong maxbits)
{
    ulong n;
    slong k = 1 + n_randint(state, maxbits);

    switch (n_randint(state, 8))
    {
        case 0:
            n = UWORD(1) << (k - 1);
            break;
        case 1:
            n = (k >= 2) ? (UWORD(1) << (k - 1)) + 1 : 1;
            break;
        case 2:
            n = (k == FLINT_BITS) ? UWORD_MAX : (UWORD(1) << k) - 1;
            break;
#if FLINT_BITS == 64
        case 3:
            n = UWORD(1515531528) + n_randint(state, 2);
            break;
        case 4:
            n = (maxbits >= 50) ? UWORD(1108307720798209) : UWORD(167772161);
            break;
        case 5:
            n = (maxbits >= 30) ? UWORD(998244353) : UWORD(36175873);
            break;
#endif
        default:
            n = n_randtest_bits(state, k);
    }

    return FLINT_MAX(n, 1);
}

/* random residues, many n - 1 */
static void
_kt_randvec(nn_ptr v, slong len, flint_rand_t state, ulong n)
{
    slong i;
    int all_max = (n_randint(state, 4) == 0);

    for (i = 0; i < len; i++)
        v[i] = (all_max || n_randint(state, 4) == 0) ? n - 1 : n_randint(state, n);
}

/* a length, often a multiple of 16 or next to one */
static slong
_kt_len(flint_rand_t state, slong max)
{
    slong l = 1 + n_randint(state, max);

    if (n_randint(state, 2))
    {
        l = 16 * (l / 16) + (slong) n_randint(state, 3) - 1;
        l = FLINT_MAX(l, 1);
    }

    return FLINT_MIN(l, max);
}

/* the block length and the number of terms of the one-word code for a
   window of a product of la and lb coefficients, after the truncation of
   the operands to nhi and the trimming of a middle product, as in
   _nmod_poly_mulmid_karatsuba; returns whether it is a middle product */
static int
_kt_blocks(slong * s, slong * terms, slong la, slong lb, slong nlo, slong nhi)
{
    slong len1 = FLINT_MIN(la, nhi), len2 = FLINT_MIN(lb, nhi), t;

    if (nlo == 0)
    {
        *s = *terms = FLINT_MIN(len1, len2);
        return 0;
    }

    t = nlo - (len2 - 1);
    if (t > 0)
    {
        len1 -= t;
        nlo -= t;
        nhi -= t;
    }
    t = nlo - (len1 - 1);
    if (t > 0)
    {
        len2 -= t;
        nlo -= t;
        nhi -= t;
    }
    if (len1 < len2)
        FLINT_SWAP(slong, len1, len2);

    *s = FLINT_MIN(len2, nhi - nlo);
    *terms = len2;
    return 1;
}

TEST_FUNCTION_START(nmod_poly_mulmid_karatsuba, state)
{
    slong iter;

    /* any window, against the classical algorithm; the operands are
       followed by words >= n, which must not be read, and the coefficients
       after the window are not written */
    for (iter = 0; iter < 1000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nn_ptr a, b, c, d;
        slong la, lb, len, nlo, nhi, s, terms, i, max;
        int square, mid;

        max = n_randint(state, 4) ? 80 : n_randint(state, 4) ? 300 : 2100;
        la = _kt_len(state, max);
        lb = n_randint(state, 4) ? _kt_len(state, max) : la + _kt_len(state, 3 * max);
        if (n_randint(state, 2))
            FLINT_SWAP(slong, la, lb);
        len = la + lb - 1;

        switch (n_randint(state, 6))
        {
            case 0: nlo = 0; nhi = len; break;
            case 1: nlo = 0; nhi = 1 + n_randint(state, len); break;
            case 2: nlo = FLINT_MIN(la, lb) - 1; nhi = FLINT_MAX(la, lb); break;
            case 3: nlo = n_randint(state, len); nhi = len; break;
            case 4:
                /* a middle product of one or two blocks of outputs */
                nlo = n_randint(state, len);
                nhi = nlo + 1 + _kt_len(state, 48);
                nhi = FLINT_MIN(nhi, len);
                break;
            default:
                nlo = n_randint(state, len);
                nhi = nlo + 1 + n_randint(state, len - nlo);
        }

        /* at the limits of the one-word code for the blocks it uses */
        mid = _kt_blocks(&s, &terms, la, lb, nlo, nhi);
        nmod_init(&mod, _kt_modulus(state, s, terms, mid));

        square = (la == lb && nlo == 0 && n_randint(state, 4) == 0);

        a = _nmod_vec_init(la + 16);
        b = _nmod_vec_init(lb + 16);
        c = _nmod_vec_init(nhi - nlo + 1);
        d = _nmod_vec_init(nhi - nlo);

        _kt_randvec(a, la, state, mod.n);
        _kt_randvec(b, lb, state, mod.n);
        for (i = 0; i < 16; i++)
        {
            a[la + i] = UWORD_MAX - i;
            b[lb + i] = UWORD_MAX - i;
        }
        c[nhi - nlo] = 12345;

        if (square)
        {
            _nmod_poly_mulmid_karatsuba(c, a, la, a, la, nlo, nhi, mod);
            _nmod_poly_mulmid_classical(d, a, la, a, la, nlo, nhi, mod);
        }
        else
        {
            _nmod_poly_mulmid_karatsuba(c, a, la, b, lb, nlo, nhi, mod);
            _nmod_poly_mulmid_classical(d, a, la, b, lb, nlo, nhi, mod);
        }

        if (!_nmod_vec_equal(c, d, nhi - nlo) || c[nhi - nlo] != 12345)
        {
            for (i = 0; i < nhi - nlo && c[i] == d[i]; i++) ;
            flint_printf("FAIL (window):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, nlo = %wd, nhi = %wd, square = %d, first difference at %wd\n",
                         mod.n, la, lb, nlo, nhi, square, i);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

    /* the doubles (fp50), forced for any n < 2^50, any window */
    for (iter = 0; iter < 500 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nn_ptr a, b, c, d;
        slong la, lb, len, nlo, nhi, i, max;
        ulong n;

        max = n_randint(state, 4) ? 80 : n_randint(state, 4) ? 300 : 1100;
        la = _kt_len(state, max);
        lb = n_randint(state, 4) ? _kt_len(state, max) : la + _kt_len(state, 3 * max);
        if (n_randint(state, 2))
            FLINT_SWAP(slong, la, lb);
        len = la + lb - 1;

        switch (n_randint(state, 4))
        {
            case 0: nlo = 0; nhi = len; break;
            case 1: nlo = 0; nhi = 1 + n_randint(state, len); break;
            case 2: nlo = FLINT_MIN(la, lb) - 1; nhi = FLINT_MAX(la, lb); break;
            default:
                nlo = n_randint(state, len);
                nhi = nlo + 1 + n_randint(state, len - nlo);
        }

#if FLINT_BITS == 64
        switch (n_randint(state, 4))
        {
            case 0:
                n = (UWORD(1) << 50) - 1 - n_randint(state, 100);
                break;
            case 1:
                n = UWORD(1) << n_randint(state, 50);
                break;
            case 2:
                n = n_randtest_bits(state, 45 + n_randint(state, 6));
                break;
            default:
                n = n_randtest_bits(state, 1 + n_randint(state, 50));
        }
#else
        n = n_randtest_not_zero(state);
#endif
        nmod_init(&mod, FLINT_MAX(n, 1));

        a = _nmod_vec_init(la + 16);
        b = _nmod_vec_init(lb + 16);
        c = _nmod_vec_init(nhi - nlo + 1);
        d = _nmod_vec_init(nhi - nlo);

        _kt_randvec(a, la, state, mod.n);
        _kt_randvec(b, lb, state, mod.n);
        for (i = 0; i < 16; i++)
        {
            a[la + i] = UWORD_MAX - i;
            b[lb + i] = UWORD_MAX - i;
        }
        c[nhi - nlo] = 12345;

        _nmod_poly_mulmid_karatsuba_fp50(c, a, la, b, lb, nlo, nhi, mod);
        _nmod_poly_mulmid_classical(d, a, la, b, lb, nlo, nhi, mod);

        if (!_nmod_vec_equal(c, d, nhi - nlo) || c[nhi - nlo] != 12345)
        {
            for (i = 0; i < nhi - nlo && c[i] == d[i]; i++) ;
            flint_printf("FAIL (fp50):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, nlo = %wd, nhi = %wd, first difference at %wd\n",
                         mod.n, la, lb, nlo, nhi, i);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

    /* the dispatchers, against the classical algorithm, for every bit size
       (the IFMA kernel above 50 bits), often at the moduli where the rules
       change */
    for (iter = 0; iter < 300 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        nn_ptr a, b, c, d;
        slong la, lb, n, nlo, nhi, kind;

        if (n_randint(state, 2))
            nmod_init(&mod, _kt_rule_modulus(state, FLINT_BITS));
        else
            nmod_init(&mod, n_randtest_bits(state,
                                    2 + n_randint(state, FLINT_BITS - 1)));

        lb = 4 + n_randint(state, n_randint(state, 2) ? 60 : 600);
        la = lb + (n_randint(state, 2) ? n_randint(state, 4) : n_randint(state, 3 * lb));
        kind = n_randint(state, 4);

        a = _nmod_vec_init(la);
        b = _nmod_vec_init(lb);
        c = _nmod_vec_init(la + lb);
        d = _nmod_vec_init(la + lb);

        _kt_randvec(a, la, state, mod.n);
        _kt_randvec(b, lb, state, mod.n);

        if (kind == 0)
        {
            _nmod_poly_mul(c, a, la, b, lb, mod);
            _nmod_poly_mul_classical(d, a, la, b, lb, mod);
            n = la + lb - 1;
        }
        else if (kind == 1)
        {
            n = 1 + n_randint(state, la + lb - 1);
            _nmod_poly_mullow(c, a, la, b, lb, n, mod);
            _nmod_poly_mullow_classical(d, a, la, b, lb, n, mod);
        }
        else if (kind == 2)
        {
            _nmod_poly_mul(c, b, lb, b, lb, mod);
            _nmod_poly_mul_classical(d, b, lb, b, lb, mod);
            n = 2 * lb - 1;
        }
        else
        {
            nlo = n_randint(state, la + lb - 1);
            nhi = nlo + 1 + n_randint(state, la + lb - 1 - nlo);
            _nmod_poly_mulmid(c, a, la, b, lb, nlo, nhi, mod);
            _nmod_poly_mulmid_classical(d, a, la, b, lb, nlo, nhi, mod);
            n = nhi - nlo;
        }

        if (!_nmod_vec_equal(c, d, n))
        {
            flint_printf("FAIL (dispatch):\n");
            flint_printf("n = %wu, la = %wd, lb = %wd, kind = %wd\n", mod.n, la, lb, kind);
            fflush(stdout);
            flint_abort();
        }

        _nmod_vec_clear(a);
        _nmod_vec_clear(b);
        _nmod_vec_clear(c);
        _nmod_vec_clear(d);
    }

#if NMOD_POLY_HAVE_KARATSUBA
    /* _nmod_poly_karatsuba_in_range, checked before the rule, never declines
       what the rule takes */
    for (iter = 0; iter < 10000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        slong len1, len2, m;
        int kind = n_randint(state, 4);

        nmod_init(&mod, _kt_rule_modulus(state, 50));
        len2 = 8 + n_randint(state, n_randint(state, 2) ? 40 : 600);
        len1 = (kind == _NMOD_POLY_KIND_SQR) ? len2
             : len2 + (slong) n_randint(state, 4 * len2);
        if (kind == _NMOD_POLY_KIND_MUL || kind == _NMOD_POLY_KIND_SQR)
            m = len1 + len2 - 1;
        else if (kind == _NMOD_POLY_KIND_LOW)
            m = 8 + n_randint(state, len1 + len2 - 9);
        else
            m = 8 + n_randint(state, 4 * len2);

        if (_nmod_poly_want_karatsuba(len1, len2, m, kind, mod.n)
                && !_nmod_poly_karatsuba_in_range(len2, m, kind, NMOD_BITS(mod)))
        {
            flint_printf("FAIL (in range):\n");
            flint_printf("n = %wu, len1 = %wd, len2 = %wd, m = %wd, kind = %d\n",
                         mod.n, len1, len2, m, kind);
            fflush(stdout);
            flint_abort();
        }
    }
#endif

    /* _nmod_poly_dispatch, from the bit size of n, agrees with the choice
       from n of the functions out of line where it decides */
    for (iter = 0; iter < 20000 * flint_test_multiplier(); iter++)
    {
        nmod_t mod;
        slong len1, len2, m, bits;
        int kind = n_randint(state, 4), d, e;

        nmod_init(&mod, _kt_rule_modulus(state, FLINT_BITS));
        bits = NMOD_BITS(mod);
        len2 = 1 + n_randint(state, n_randint(state, 2) ? 64 : 1100);
        len1 = (kind == _NMOD_POLY_KIND_SQR) ? len2
             : len2 + (slong) (n_randint(state, 2) ? n_randint(state, 4)
                                           : n_randint(state, 4 * len2));
        if (kind == _NMOD_POLY_KIND_LOW && len2 < 2)
            kind = _NMOD_POLY_KIND_MUL;
        if (kind == _NMOD_POLY_KIND_MUL || kind == _NMOD_POLY_KIND_SQR)
            m = len1 + len2 - 1;
        else if (kind == _NMOD_POLY_KIND_LOW)
            m = len1 + n_randint(state, len2 - 1);
        else
            m = 6 + n_randint(state, FLINT_MAX(1, len1 + len2 - 7));

        d = _nmod_poly_dispatch(len1, len2, m, kind, bits);

        e = _nmod_poly_dispatch_n(len1, len2, m, kind, mod);

        if (d != _NMOD_POLY_DISPATCH_BY_N && d != e)
        {
            flint_printf("FAIL (dispatch from the bit size):\n");
            flint_printf("n = %wu, len1 = %wd, len2 = %wd, m = %wd, kind = %d, %d instead of %d\n",
                         mod.n, len1, len2, m, kind, d, e);
            fflush(stdout);
            flint_abort();
        }
    }

    /* the same, exhaustively for short lengths at the moduli where the
       rules change */
    {
        ulong ns[16];
        slong nn = 0, in, len1, len2, m, j;
        int kind, d, e;

        ns[nn++] = 2;
        ns[nn++] = UWORD(1) << 20;
        ns[nn++] = (UWORD(1) << 30) - 1;
        ns[nn++] = UWORD(1) << 30;
#if FLINT_BITS == 64
        ns[nn++] = UWORD(1515531528);
        ns[nn++] = UWORD(1515531529);
        ns[nn++] = (UWORD(1) << 32) + 1;
        ns[nn++] = UWORD(1108307720798209);
        ns[nn++] = (UWORD(1) << 50) - 1;
        ns[nn++] = UWORD(1) << 52;
        ns[nn++] = (UWORD(1) << 52) + 1;
        ns[nn++] = UWORD_MAX;
#endif

        for (in = 0; in < nn; in++)
        {
            nmod_t mod;

            nmod_init(&mod, ns[in]);

            for (len2 = 1; len2 <= 40; len2++)
            for (j = 0; j < 4; j++)
            for (kind = 0; kind < 4; kind++)
            {
                len1 = (kind == _NMOD_POLY_KIND_SQR) ? len2
                     : (j == 0) ? len2 : (j == 1) ? len2 + 3
                     : (j == 2) ? 2 * len2 + 1 : 4 * len2;

                for (m = 1; m <= len1 + len2 - 1; m++)
                {
                    if ((kind == _NMOD_POLY_KIND_MUL
                            || kind == _NMOD_POLY_KIND_SQR)
                            && m != len1 + len2 - 1)
                        continue;
                    if (kind == _NMOD_POLY_KIND_LOW
                            && (m < len1 || m == len1 + len2 - 1))
                        continue;

                    d = _nmod_poly_dispatch(len1, len2, m, kind, NMOD_BITS(mod));
                    e = _nmod_poly_dispatch_n(len1, len2, m, kind, mod);

                    if (d != _NMOD_POLY_DISPATCH_BY_N && d != e)
                    {
                        flint_printf("FAIL (dispatch from the bit size, boundary moduli):\n");
                        flint_printf("n = %wu, len1 = %wd, len2 = %wd, m = %wd, kind = %d, %d instead of %d\n",
                                     mod.n, len1, len2, m, kind, d, e);
                        fflush(stdout);
                        flint_abort();
                    }
                }
            }
        }
    }

    TEST_FUNCTION_END(state);
}
