/*=============================================================================

    This file is part of FLINT.

    FLINT is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    FLINT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with FLINT; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301 USA

=============================================================================*/
/******************************************************************************

    Copyright (C) 2010 Fredrik Johansson

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "mpn_extras.h"

/*
The basic approach to compute H(n) quickly is to use a balanced sum
S(a,b) = S(a,m) + S(m,b) where m = floor((a+b)/2), analogous to using a
balanced product for n factorial. To reduce overhead, we can work with
unnormalized numerators and denominators and perform a final GCD reduction
only at the end. In fact, this way the computation of the denominator
(before the final GCD reduction) exactly amounts to a balanced product
for n factorial.

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

Note on memory allocation:

The maximum size in bits required for each partial numerator/denominator
between a and b can be bounded by bits(n)*(b-a+2), since each
denominator product is bounded by bits(n)*(b-a+1) and H_n is bounded
in magnitude by bits(n). This assumes that we don't introduce any spurious
factors.

Some more optimization is possible:
+ At the bottom level, we could add two terms in one step most of the
  time since 1/k + 1/(k+2) will usually have single-limb
  numerator and denominator.
+ More generally, we could remove multiples of 3, 5, ... from the summation,
  although this would result in a considerably more complicated algorithm, and
  the returns would probably diminish quickly.
+ Some more trailing zeros could be eliminated from numerators/denominators.
  However, these account for less than 1% of the operand sizes for reasonably
  large inputs.
+ It may be faster to normalize some of the partial sums.
+ The memory management could be improved.

*/

void
flint_mpn_harmonic_odd_direct(mp_ptr t, mp_size_t * tsize,
                        mp_ptr v, mp_size_t * vsize,
                        len_t a, len_t b, len_t n, int d)
{
    mp_size_t ts, vs;

    *t = 0UL;
    *v = 1UL;
    ts = 1;
    vs = 1;

    if (a == 1)
    {
        mp_limb_t r, s;
        len_t k;
        for (k = b - 1 - (b % 2); k > 0; k -= 2)
        {
            while (k <= (n >> d))
                d += 1;
            r = ((mp_limb_t) k) << (d-1);
            s = (1UL << d) - 1UL;
            MPN_MUL_1(t, ts, t, ts, r);
            MPN_ADDMUL_1(t, ts, v, vs, s);
            MPN_MUL_1(v, vs, v, vs, r);
        }
    }
    else
    {
        a += (a % 2 == 0);
        for ( ; a < b; a += 2)
        {
            MPN_MUL_1(t, ts, t, ts, a);
            MPN_ADD(t, ts, t, ts, v, vs);
            MPN_MUL_1(v, vs, v, vs, a);
        }
        MPN_MUL_1(t, ts, t, ts, (1UL<<d) - 1UL);
        MPN_MUL_1(v, vs, v, vs, 1UL << (d-1));
    }

    *tsize = ts;
    *vsize = vs;
}

void
flint_mpn_harmonic_odd_balanced(mp_ptr t, mp_size_t * tsize,
                          mp_ptr v, mp_size_t * vsize,
                          len_t a, len_t b, len_t n, int d)
{
    mp_ptr p, q, r, s;
    mp_size_t ps, qs, rs, ss, ts, vs;
    len_t m, tmpsize;

    if (b - a < 50)
    {
        flint_mpn_harmonic_odd_direct(t, tsize, v, vsize, a, b, n, d);
        return;
    }

    m = (a + b) / 2;

    tmpsize = sizeof(mp_limb_t) * FLINT_BIT_COUNT(n) * (b-a+2) / FLINT_BITS;
    p = flint_malloc(tmpsize);  /* Twice because we re-use it for q*r */
    q = flint_malloc(tmpsize / 2 + 1);
    r = flint_malloc(tmpsize / 2 + 1);
    s = flint_malloc(tmpsize / 2 + 1);

    flint_mpn_harmonic_odd_balanced(p, &ps, q, &qs, a, m, n, d + (a == 1));
    flint_mpn_harmonic_odd_balanced(r, &rs, s, &ss, m, b, n, d);

    MPN_MUL(t, ts, p, ps, s, ss);
    MPN_MUL(p, ps, q, qs, r, rs);
    MPN_ADD(v, vs, t, ts, p, ps);
    MPN_SET(t, ts, v, vs);
    MPN_MUL(v, vs, q, qs, s, ss);

    *tsize = ts;
    *vsize = vs;

    flint_free(p);
    flint_free(q);
    flint_free(r);
    flint_free(s);
}
