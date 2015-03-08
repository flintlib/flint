#include <stdio.h>
#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>


mp_limb_t n_rootrem(mp_limb_t * r, mp_limb_t x, unsigned int k)       /* taken from the link, you gave */
{
    mp_limb_t a, b, m, p, hi, lo;
    unsigned int bc;
    int i;

    if (k <= 1 || x <= 1)
    {
        *r = 0;
        return x;
    }

    if (k == 2)
    {
        return n_sqrtrem(r, x);
    }

    if (k >= FLINT_BITS || (UWORD(1) << k) > x)
    {
        *r = x - 1;
        return 1;
    }

    a = 2;
    b = UWORD(1) << ((FLINT_BIT_COUNT(x) + k - 1) / k);

    while (a < b)
    {
        m = a + (b - a) / 2;

        p = m + 1;

        for (i = 1; i < k; i++)
        {
            umul_ppmm(hi, lo, p, m + 1);

            if (hi != 0)
                goto too_large;
            else
                p = lo;
        }

        if (p == x)
        {
            *r = 0;
            return m + 1;
        }
        else if (p > x)
        {
            too_large:
            b = m;
        }
        else
        {
            a = m + 1;
        }
    }

    p = a;
    for (i = 1; i < k; i++)
        p *= a;

    *r = x - p;
    return a;
}

mp_limb_t n_is_perfect_power(mp_limb_t * r, mp_limb_t x)
{
    int p[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};                       /* prime less than 64 */
    int q[] = {3, 7, 11, 29, 23, 53, 103, 191, 277, 349, 373, 593, 739, 947, 1129, 1697, 1889, 1831};     /* prime such that, q[i] % p[i] = 1 */
    int i, b = FLINT_BIT_COUNT(x);                                                                                                /* if x = y ^ p[i] for some y => (y ^ (q[i] - 1)) % q[i] = 1 or 0 */
    double qinv;
    mp_limb_t root;

    for (i = 0; i < 18 && p[i] <= b; i++)
    {
        qinv = n_precompute_inverse(q[i]);
        if (n_powmod_precomp(x, (q[i] - 1) / p[i], q[i], qinv) <= 1)
        {
            root = n_rootrem(r, x, p[i]);
            if (*r == 0)
            {
                return root;
            }
        }
    }
    return 0;
}

