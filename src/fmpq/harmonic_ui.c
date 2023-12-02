/*
    Copyright (C) 2011, 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpq.h"

#if FLINT_BITS == 64
#define FMPQ_HARMONIC_UI_TAB_SIZE 47
#else
#define FMPQ_HARMONIC_UI_TAB_SIZE 25
#endif

static const mp_limb_t fmpq_harmonic_ui_tab_num[] =
{
    0, 1, 3, 11, 25, 137, 49, 363, 761, 7129, 7381, 83711, 86021, 1145993,
    1171733, 1195757, 2436559, 42142223, 14274301, 275295799, 55835135,
    18858053, 19093197, 444316699, 1347822955,
#if FLINT64
    UWORD(34052522467), UWORD(34395742267), UWORD(312536252003),
    UWORD(315404588903), UWORD(9227046511387), UWORD(9304682830147),
    UWORD(290774257297357), UWORD(586061125622639), UWORD(53676090078349),
    UWORD(54062195834749), UWORD(54437269998109), UWORD(54801925434709),
    UWORD(2040798836801833), UWORD(2053580969474233), UWORD(2066035355155033),
    UWORD(2078178381193813), UWORD(85691034670497533), UWORD(12309312989335019),
    UWORD(532145396070491417), UWORD(5884182435213075787),
    UWORD(5914085889685464427), UWORD(5943339269060627227),
#endif
};

const mp_limb_t fmpq_harmonic_ui_tab_den[] =
{
    1, 1, 2, 6, 12, 60, 20, 140, 280, 2520, 2520, 27720, 27720, 360360,
    360360, 360360, 720720, 12252240, 4084080, 77597520, 15519504, 5173168,
    5173168, 118982864, 356948592,
#if FLINT64
    UWORD(8923714800), UWORD(8923714800), UWORD(80313433200),
    UWORD(80313433200), UWORD(2329089562800), UWORD(2329089562800),
    UWORD(72201776446800), UWORD(144403552893600), UWORD(13127595717600),
    UWORD(13127595717600), UWORD(13127595717600), UWORD(13127595717600),
    UWORD(485721041551200), UWORD(485721041551200), UWORD(485721041551200),
    UWORD(485721041551200), UWORD(19914562703599200), UWORD(2844937529085600),
    UWORD(122332313750680800), UWORD(1345655451257488800),
    UWORD(1345655451257488800), UWORD(1345655451257488800),
#endif
};

/*
The basic approach to compute H(n) quickly is to use a balanced sum.
To save some more time, we note that we only have to sum over the odd
terms since H(n) = H(floor(n/2))/2 + H_odd(n). Recursive application of this
formula results in a geometric series for the weight of each odd term 1/k:

    n/2   < k <= n         : weight 1
    n/4   < k <= n/2       : weight 3/2
    n/8   < k <= n/4       : weight 7/4
    n/16  < k <= n/8       : weight 15/8
         ...
    n/2^d < k <= n/2^(d-1) : weight (2^d-1)/2^(d-1)

Although not necessary, the implementation is simplified by always splitting
the interval exactly in half, since we then just have to increment d on every
subinterval that starts with a = 1. Below a threshold, we fall back to direct
summation of the odd fractions.

As a final optimization, we accumulate word-size partial sums in
single limbs in the basecase summation.

A basic Python implementation:

def harmonic_odd_direct(a, b, n, d):
    t, v = 0, 1
    if a == 1:
        for k in range(b-1-(b%2), 0, -2):
            while k <= (n >> d):
                d += 1
            r = 2**(d-1)*k
            t, v = ((2**d-1)*v + r*t), r*v
        return t, v
    else:
        a += (a % 2 == 0)
        for k in range(a, b, 2):
            t, v = (v+k*t), k*v
        return (2**d - 1) * t, 2**(d-1) * v

def harmonic_odd_balanced(a, b, n, d):
    if b - a < 50:
        return harmonic_odd_direct(a, b, n, d)
    m = (a+b) // 2
    t, v = harmonic_odd_balanced(a, m, n, d + (a==1))
    u, w = harmonic_odd_balanced(m, b, n, d)
    return (t*w + u*v), v*w

def harmonic(n):
    return harmonic_odd_balanced(1, n+1, n, 1)

*/

static void
harmonic_odd_direct(fmpz_t P, fmpz_t Q, ulong a, ulong b, ulong n, int d)
{
    mp_limb_t p, q, r, s, t, u, v, w = 0;
    slong k;

    fmpz_zero(P);
    fmpz_one(Q);

    p = 0;
    q = 1;

    if (a == 1)
    {
        for (k = b - 1 - (b % 2); k > 0; k -= 2)
        {
            while (k <= (n >> d))
                d++;

            r = (UWORD(1) << d) - UWORD(1);
            s = ((mp_limb_t) k) << (d-1);

            umul_ppmm(t, u, p, s);
            umul_ppmm(v, w, q, r);

            if (t == 0 && v == 0)
            {
                add_ssaaaa(t, u, t, u, v, w);

                if (t == 0)
                    umul_ppmm(v, w, q, s);
            }

            if (t == 0 && v == 0)
            {
                p = u;
                q = w;
            }
            else
            {
                fmpz_mul_ui(P, P, q);
                fmpz_addmul_ui(P, Q, p);
                fmpz_mul_ui(Q, Q, q);

                p = r;
                q = s;
            }
        }

        if (p != 0)
        {
            fmpz_mul_ui(P, P, q);
            fmpz_addmul_ui(P, Q, p);
            fmpz_mul_ui(Q, Q, q);
        }
    }
    else
    {
        a += (a % 2 == 0);

        for (k = a; k < b; k += 2)
        {
            umul_ppmm(t, u, p, k);
            v = 0;

            if (t == 0)
            {
                add_ssaaaa(t, u, t, u, 0, q);
                if (t == 0)
                    umul_ppmm(v, w, q, k);
            }

            if (t == 0 && v == 0)
            {
                p = u;
                q = w;
            }
            else
            {
                fmpz_mul_ui(P, P, q);
                fmpz_addmul_ui(P, Q, p);
                fmpz_mul_ui(Q, Q, q);

                p = 1;
                q = k;
            }
        }

        if (p != 0)
        {
            fmpz_mul_ui(P, P, q);
            fmpz_addmul_ui(P, Q, p);
            fmpz_mul_ui(Q, Q, q);
        }

        fmpz_mul_ui(P, P, (UWORD(1) << d) - UWORD(1));
        fmpz_mul_ui(Q, Q, UWORD(1) << (d - 1));
    }
}

static void
harmonic_odd_balanced(fmpz_t P, fmpz_t Q, ulong a, ulong b, ulong n, int d)
{
    if (b - a < 50)
    {
        harmonic_odd_direct(P, Q, a, b, n, d);
    }
    else
    {
        ulong m;
        fmpz_t R, S;

        fmpz_init(R);
        fmpz_init(S);

        m = a + (b - a) / 2;

        harmonic_odd_balanced(P, Q, a, m, n, d + (a==1));
        harmonic_odd_balanced(R, S, m, b, n, d);

        fmpz_mul(P, P, S);
        fmpz_addmul(P, Q, R);
        fmpz_mul(Q, Q, S);

        fmpz_clear(R);
        fmpz_clear(S);
    }
}

void
_fmpq_harmonic_ui(fmpz_t num, fmpz_t den, ulong n)
{
    if (n < FMPQ_HARMONIC_UI_TAB_SIZE)
    {
        fmpz_set_ui(num, fmpq_harmonic_ui_tab_num[n]);
        fmpz_set_ui(den, fmpq_harmonic_ui_tab_den[n]);
    }
    else
    {
        /* overflow */
        if ((slong) n < 0)
            flint_throw(FLINT_ERROR, "(%s)\n", __func__);

        harmonic_odd_balanced(num, den, 1, n + 1, n, 1);
        _fmpq_canonicalise(num, den);
    }
}

void
fmpq_harmonic_ui(fmpq_t x, ulong n)
{
    _fmpq_harmonic_ui(fmpq_numref(x), fmpq_denref(x), n);
}

