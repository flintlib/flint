/*
    Copyright (C) 2026 Fredrik Johansson
    Developed using Claude Fable 5.1

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <math.h>
#include "ulong_extras.h"
#include "fmpz.h"
#include "ecpp.h"

/*
    Discriminant pool for ECPP (in the spirit of fastECPP, Franke,
    Kleinjung, Morain, Wirth 2004).

    A fundamental discriminant D < 0 factors as D = q0 * prod p^*, where
    p^* = (-1)^{(p-1)/2} p runs over the odd primes dividing D and q0 is 1,
    -4, 8 or -8. The Kronecker symbols (p^* / n), (q0 / n) are the genus
    characters. n is represented by the principal form of discriminant D
    (which is what Cornacchia decides) only if all of them are 1, and in
    that case it is with probability 2^{g-1} / h, where g is the number of
    characters: one over the odd part o = h / 2^{g-1} of the class number.

    Consequently the square root of D modulo n is only ever needed when
    every p^* dividing D is a square modulo n, and it is the product of the
    square roots of the p^*, which are computed once per n. The pool
    therefore consists of the fundamental discriminants that are smooth
    over a fixed set of small primes: their square roots cost nothing
    beyond that fixed set, and a pool of thousands of D's is available
    for every n.

    The pool is sorted by increasing cost of realising a step with D: by
    number of units (D = -3, -4 first: more curves per solution), then
    odd part o (the degree of the factor of the class polynomial over the
    genus field whose root is needed, see genus_poly.c, and the inverse of
    the probability of a Cornacchia success), then class number, then |D|.
*/

/*
    Estimated cost of realising a step with D, in units of one scalar
    multiplication on the curve at a few thousand bits: the twists, the
    class polynomial numerics of the class field tower (cheap with the
    Weber invariant, i.e. for even D; about 30 units per hundred classes
    with j), and the descent through the prime factors p of h: a square
    root for p = 2, radicals or a small powering for p = 3, powerings of
    degree p otherwise. Class numbers with a prime factor above 13 are
    not used (cost above ECPP_DISC_COST_MAX).
*/
double
ecpp_disc_cost_v(const ecpp_disc_struct * d, int veven)
{
    return ecpp_disc_cost_n(d, veven, NULL, 4000);
}

/*
    With nmodp = {n mod 5, n mod 7, n mod 11, n mod 13} the actual descent
    costs are used: a level of prime degree p is a single p-th root (in
    F_n or F_{n^2}) when n = +-1 mod p, and a powering of a degree-p
    polynomial otherwise; without nmodp the average over n.
*/
double
ecpp_disc_cost_n(const ecpp_disc_struct * d, int veven, const ulong * nmodp, slong bits)
{
    /*
        The class polynomial numerics do not depend on n, while the unit
        (a scalar multiplication) shrinks with n: their relative cost at
        bits is that at 4000 bits times (4000 / bits)^2.6, so that for
        small n only small class numbers are worth it.
    */
    double numfac = pow(4000.0 / FLINT_MAX(bits, 200), 2.6);
    numfac = FLINT_MIN(numfac, 500.0);
    numfac = FLINT_MAX(numfac, 1.0);
    slong h = d->h, p;
    ulong absD = (ulong) (-d->D), m = absD >> 2;
    /*
        Weber's function: even D directly; odd D through the order of
        conductor 2 (discriminant 4D) when the v of 4n = t^2 + |D| v^2 is
        even, since then a curve with that endomorphism ring and the right
        order exists (any curve of the right order serves the proof). For
        D = 1 mod 8 the order has the same class number and v is always
        even (t^2 + |D| v^2 = 4n with |D| = 7 mod 8 forces t, v even); for
        D = 5 mod 8 the class number is three times that of D and v is
        even half of the time.
    */
    int weber = (((absD & 3) == 0) && ((m & 7) != 4 && (m & 7) != 0))
                || (((absD & 1) != 0) && veven);
    double cost = 1.5, genus;

    /* the twists are tried in turn: 6 for D = -3, 4 for D = -4, 2 else */
    if (d->D == -3)
        cost = 1.0 + 3.5;
    else if (d->D == -4)
        cost = 1.0 + 2.5;
    if (h == 1)
        return cost;
    if (weber && (absD & 1) != 0 && (absD & 7) == 3)
        h *= 3;
    /* the alternative: a root of the factor of H_D over the genus field,
       of degree o (by radicals up to 4) */
    genus = 1.5 + ((d->o <= 4) ? 0.3 * d->o : 1.4 * d->o + 4.0);
    /* class polynomial numerics: with the Weber invariant about 0.13 s for
       h = 256, with j about 3.5 s for h = 230 (dominated by the
       evaluations of j at high precision), so j is only worth it for
       small class numbers */
    if (weber)
        cost += 0.03 * h * numfac;
    else if (h <= 48)
        cost += 1.0 * h * numfac;
    else
        cost += 1e6;
    for (p = 2; h > 1; p++)
    {
        while (h % p == 0)
        {
            double powmod, cheap = 1.0, prob;
            slong idx = (p == 5) ? 0 : (p == 7) ? 1 : (p == 11) ? 2 : (p == 13) ? 3 : -1;
            h /= p;
            if (p == 2) { cost += 0.3; continue; }
            if (p == 3) { cost += 1.5; continue; }
            if (p > 13) { cost += 1e6; continue; }
            powmod = (p == 5) ? 6.0 : (p == 7) ? 10.0 : 3.0 * p;
            if (nmodp != NULL)
                cost += (bits >= ECPP_KUMMER_BITS && (nmodp[idx] == 1 || nmodp[idx] == (ulong) p - 1)) ? cheap : powmod;
            else
            {
                prob = 2.0 / (p - 1);       /* n = +-1 mod p */
                cost += prob * cheap + (1 - prob) * powmod;
            }
        }
    }
    return FLINT_MIN(cost, genus);
}

/* the optimistic estimate (v even), used to order and admit the pool */
double
ecpp_disc_cost(const ecpp_disc_struct * d)
{
    return ecpp_disc_cost_v(d, 1);
}

/* whether the class field tower (with the Weber invariant when possible)
   is the cheaper way to a curve for this discriminant; sets *Dt to the
   discriminant to build the tower on (D or 4D) */
int
ecpp_disc_use_tower(const ecpp_disc_struct * d, int veven, slong * Dt)
{
    ulong absD = (ulong) (-d->D), m = absD >> 2;
    int weber = (((absD & 3) == 0) && ((m & 7) != 4 && (m & 7) != 0))
                || (((absD & 1) != 0) && veven);
    *Dt = (((absD & 1) != 0) && veven) ? 4 * d->D : d->D;
    return d->h >= 2 && (weber || d->h <= 48);
}

static int
_disc_cmp(const void * x, const void * y)
{
    const ecpp_disc_struct * a = (const ecpp_disc_struct *) x;
    const ecpp_disc_struct * b = (const ecpp_disc_struct *) y;
    int wa = (a->D == -3) ? 6 : (a->D == -4) ? 4 : 2;
    int wb = (b->D == -3) ? 6 : (b->D == -4) ? 4 : 2;

    if (wa != wb)
        return (wa > wb) ? -1 : 1;
    if (a->o != b->o)
        return (a->o < b->o) ? -1 : 1;
    if (a->cost != b->cost)
        return (a->cost < b->cost) ? -1 : 1;
    return (a->D > b->D) ? -1 : (a->D < b->D);
}

/*
    Fills *table with the fundamental discriminants D, -Dmax <= D < 0,
    whose odd prime factors are among the nprimes odd primes in primes
    (increasing), with class number h(D) <= hmax and odd part o <= omax,
    sorted as described above. Returns the number of entries; the table is allocated with
    flint_malloc.

    Class numbers are obtained by counting reduced forms
    (a, b, c) of discriminant b^2 - 4ac = D, i.e. |b| <= a <= c with
    b >= 0 when |b| = a or a = c, in one pass over all (a, b, c).
*/
slong
ecpp_disc_table(ecpp_disc_struct ** table, const ulong * primes, slong nprimes,
                            slong Dmax, slong hmax, slong omax, double costmax)
{
    unsigned short * h;
    unsigned int * rem;
    slong a, b, D, amax, num, i, k;
    ecpp_disc_struct * t;

    h = flint_calloc(Dmax + 1, sizeof(unsigned short));
    rem = flint_malloc((Dmax + 1) * sizeof(unsigned int));

    /*
        Class numbers: reduced forms (a, b, c), |b| <= a <= c, with b >= 0
        when |b| = a or a = c; for fixed (a, b) the discriminants
        |D| = 4ac - b^2, c >= a, form an arithmetic progression with step
        4a. Forms with 0 < |b| < a < c count twice (b and -b). The range
        of D is processed in chunks that fit in the cache, as in Enge's CM
        library, which is what makes the sieve fast. Forms of fundamental
        discriminant are primitive, and only those are used, so no gcd
        check.
    */
    {
        const slong chunk = WORD(1) << 18;
        slong lo;

        amax = (slong) sqrt((double) Dmax / 3.0) + 1;
        for (lo = 0; lo <= Dmax; lo += chunk)
        {
            slong hi = FLINT_MIN(lo + chunk, Dmax + 1);
            unsigned short * hh = h;

            for (a = 1; a <= amax; a++)
            {
                slong step = 4 * a, D0;

                /* b = 0 and b = a: once each, c >= a */
                for (b = 0; b <= a; b += a)
                {
                    D0 = 4 * a * a - b * b;
                    if (D0 >= hi)
                        continue;
                    if (D0 < lo)
                        D0 += step * ((lo - D0 + step - 1) / step);
                    for (D = D0; D < hi; D += step)
                        hh[D]++;
                    if (a == 0)
                        break;
                }
                /* 0 < b < a: c = a once, c > a twice */
                for (b = 1; b < a; b++)
                {
                    D0 = 4 * a * a - b * b;
                    if (D0 >= hi)
                        continue;
                    if (D0 >= lo)
                    {
                        hh[D0]++;
                        D0 += step;
                    }
                    else
                        D0 += step * ((lo - D0 + step - 1) / step);
                    for (D = D0; D < hi; D += step)
                        hh[D] += 2;
                }
            }
        }
    }

    /* smooth part: rem[D] = D with the primes of the set divided out */
    for (D = 0; D <= Dmax; D++)
        rem[D] = (unsigned int) D;
    for (D = 2; D <= Dmax; D += 2)
        while (rem[D] % 2 == 0)
            rem[D] /= 2;
    for (i = 0; i < nprimes; i++)
    {
        ulong p = primes[i];
        for (D = p; D <= Dmax; D += p)
            while (rem[D] % p == 0)
                rem[D] /= p;
    }

    /*
        Fundamental: -D = 1 mod 4 squarefree, or -D = 4m with m = 2, 3
        mod 4 squarefree. Squarefreeness of the smooth ones is checked
        while factoring below.
    */
    num = 0;
    for (D = 3; D <= Dmax; D++)
        if (h[D] > 0 && h[D] <= hmax && rem[D] == 1
                && (D % 4 == 3 || (D % 4 == 0 && ((D / 4) % 4 == 1 || (D / 4) % 4 == 2))))
            num++;

    t = flint_malloc(FLINT_MAX(num, 1) * sizeof(ecpp_disc_struct));
    num = 0;
    for (D = 3; D <= Dmax; D++)
    {
        slong m, g, prodstar, o;

        if (!(h[D] > 0 && h[D] <= hmax && rem[D] == 1))
            continue;
        if (D % 4 == 3)
            m = D;
        else if (D % 4 == 0 && ((D / 4) % 4 == 1 || (D / 4) % 4 == 2))
            m = D / 4;
        else
            continue;

        /* factor m; the odd part of m must be squarefree */
        t[num].nfac = 0;
        g = 0;
        prodstar = 1;
        if (m % 2 == 0)
            m /= 2;
        for (i = 0; m > 1; i++)
        {
            slong p = primes[i];
            if (m % p == 0)
            {
                m /= p;
                if (m % p == 0)
                    break;
                if (t[num].nfac == ECPP_DISC_MAXFAC)
                    break;
                t[num].fac[t[num].nfac++] = (unsigned short) i;
                prodstar *= (p % 4 == 1) ? p : -p;
                g++;
            }
        }
        if (m != 1)
            continue;

        t[num].D = -D;
        t[num].h = h[D];
        t[num].q0 = (-D) / prodstar;    /* 1, -4, 8 or -8 */
        if (t[num].q0 != 1)
            g++;
        /* o = h / 2^{g-1} (h is divisible by 2^{g-1} for fundamental D) */
        o = h[D];
        for (k = 1; k < g; k++)
            o = (o + 1) / 2;
        t[num].g = g;
        t[num].o = FLINT_MAX(o, 1);
        t[num].cost = ecpp_disc_cost(t + num);
        /* admitted if the genus factor is small, or if the class field
           tower is cheap: class number with small prime factors only */
        if (t[num].o > omax && !(t[num].cost < costmax))
            continue;
        num++;
    }

    qsort(t, num, sizeof(ecpp_disc_struct), _disc_cmp);

    flint_free(h);
    flint_free(rem);
    *table = t;
    return num;
}
