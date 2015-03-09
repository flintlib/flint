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
    int i, b = FLINT_BIT_COUNT(x);
    mp_limb_t root;

    for (i = b; i > 0; i--)
    {
        root = n_rootrem(r, x, i);
        if (*r == 0)
        {
            *r = root;                                              /* set r to root */
            return i;                                               /* return power */
        }
    }
    return 0;
}
