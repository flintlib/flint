#include <stdio.h>
#include <stdlib.h>
#include <flint.h>
#include <ulong_extras.h>


/* program computes primitive root for the 2, 4, p,  p ^ k, 2 * (p ^ k), where p is an odd-prime & k >= 1 */

int n_is_perfect_power(mp_limb_t * root, mp_limb_t n)
{
    int i, k;
    mp_limb_t a;

    if (n_is_prime(n))
    {
        *root = n;                                  /* check if n is a prime*/
        return 1;
    }

    for (i = 3; i <= sqrt(n); i+= 2)               /*iterating through possible roots*/
    {
        a = i;
        for (k = 2; a < n; k++) a *= i;
        if (a == n)
        {
            *root = i;
            return k;
        }
    }
    return 0;
}

mp_limb_t n_primitive_root_prefactor(mp_limb_t n, n_factor_t * factors, mp_limb_t * phi)
{
    slong i;
    int found;
    mp_limb_t result, a, pm1;
    double pinv;

    pm1 = *phi;
    pinv = n_precompute_inverse(n);

    for (a = 2; a < n; a++)
    {
        if (n_gcd(n, a) == 1)
        {
            found = 1;
            for (i = 0; i < factors->num; i++)
            {
                result = n_powmod_precomp(a, pm1 / factors->p[i], n, pinv);
                if (result == 1)
                {
                    found = 0;
                    break;
                }
            }
            if (found)
            {
                return a;
            }
        }
    }
    flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
    abort();
}

mp_limb_t n_primitive_root(mp_limb_t n)
{
    if (n_is_prime(n)) return n_primitive_root_prime(n);
    else
    {
        if (n == 4) return 3;
        else if (n % 4 == 0)
        {
            flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
            return 0;
        }
        else
        {
            mp_limb_t a, m, root, phi;
            int k;
            n_factor_t factors;
            n_factor_init(&factors);
            if (n % 2 == 0)
            {
                m = 2;
                n /= 2;
            }
            else
            {
                m = 1;
            }
            k = n_is_perfect_power(&root, n);
            if (k == 0 || n_is_prime(root) != 1)
            {
                flint_printf("Exception (n_primitive_root_prime_prefactor).  root not found.\n");
                return 0;
            }
            else
            {
                phi = n;
                phi = (phi / root) * (root - 1);                      /* calculating  phi(n), using formula */
                n_factor(&factors, root - 1, 1);                      /* factoring (root - 1), which is a factor of phi(n) */
                n *= m;
                if (k > 1)                                             /* adding extra factor if exist */
                {
                    factors.num++;
                    factors.p[factors.num] = root;
                    factors.exp[factors.num] = 1;
                }
                a = n_primitive_root_prefactor(n, &factors, &phi);
                return a;
            }
        }
    }
}
