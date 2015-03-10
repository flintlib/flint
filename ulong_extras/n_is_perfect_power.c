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

int n_is_perfect_power(mp_limb_t * r, mp_limb_t x)
{
    int p[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61};                       /* prime less than 64 */
    int q[] = {3, 7, 11, 29, 23, 53, 103, 191, 47, 59, 311, 149, 83, 173, 283, 107, 709, 367};     /* prime such that, q[i] % p[i] = 1 */
    int z[] = {1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1, 1, 1, 4, 1, 1, 1, 1};                            /* if x = y ^ p[i] for some y => (y ^ (q[i] - 1)) % q[i] = 1 or 0 */
    int v[] = {0, 1, 2, 3, 2, 5, 6, 7, 2, 3, 10, 11, 12, 13, 14, 15, 2, 17, 18, 19, 20};                  /* v[i] is root of i and z[i] is it's power */
    mp_limb_t qmodulo[18][25] = { { 3 }, { 67 }, { 1027 }, { 268570627 }, { 4194307 }, { 1082130434 , 1048577 },
                            { 2 , 50380800 , 0 , 65 }, { 130 , 131200 , 262144 , 8192 , 16793600 , 1090519041 },
                            { 2 , 16385 }, { 2 , 67108865 }, { 66 , 1048592 , 2147483648 , 0 , 0 , 0 , 16777216 , 0 , 524296 , 4325377 },
                            { 2 , 4096 , 0 , 512 , 1048577 }, { 2 , 0 , 262145 }, { 2 , 0 , 536936448 , 0 , 0 , 4097 },
                            { 2 , 12288 , 0 , 0 , 0 , 0 , 0 , 49152 , 67108865 }, { 2 , 0 , 0 , 1025 },
                            { 2 , 0 , 134217729 , 0 , 0 , 134217728 , 0 , 24 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 6 , 1024 , 0 , 0 , 1056 , 0 , 0 , 17 },
                            { 2 , 0 , 1572864 , 0 , 0 , 0 , 0 , 0 , 402653184 , 0 , 0 , 16385 }

                          };             /* kth bit of qmodulo[i][j] is 1 if  y = 32 * (j - 1) + k, y ^ ((p[i] - 1) / q[i])  % q[i] is <= 1,  else 0 */
                                        /* i.e q[i][0]th elements bits from (LSB to MSB )represent 0.1.2....31, q[i][1]th represent 32.33....64 residue class modulo i and so on */

    int i, power = 1,  j = 1, pos, k;
    mp_limb_t root, prev;

    for (i = 0; i < 18; i++)
    {
        if (x <= 20)
        {
            power *= z[x];
            x = v[x];
            break;
        }

        k = x % q[i];                                     /* class of x congruent modulo q[i] */
        pos = k / 32 + (k % 32 == 0) ? 0 : 1;            /* position of correponding class of q[i] in qmodulo table */

        if (pos == 0) pos = 1;                           /* in case q[i] divides x */

        if ((qmodulo[i][pos - 1] & (1 << (k % 32))))     /* check if bit corresponding to k is set*/
        {
            prev = x;
            x = n_rootrem(r, x, p[i]);
            if (*r == 0)
            {
                j = 1;                                                       /* hold subsequent power of p[i] */
                while (*r == 0)
                {                                                            /* store product of all power that satisfy test */
                    j = j * p[i];
                    prev = x;
                    x =  n_rootrem(r, x, p[i]);                              /* test for subsequent power of prime p[i] */
                }
                power *= j;
                x = prev;
            }
            else x = prev;
        }
    }

    *r = x;                                                             /* store root in  r */
    return power;                                                       /* return power */
}

