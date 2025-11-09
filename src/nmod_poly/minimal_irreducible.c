/*
    Copyright (C) 2025 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include <stdint.h>
#include "ulong_extras.h"
#include "nmod.h"
#include "nmod_vec.h"
#include "nmod_poly.h"
#include "nmod_poly/impl.h"
#include "nmod_poly_factor.h"

static ulong
n_multiplicative_order(ulong x, ulong p, ulong pinv, n_factor_t * p1fac)
{
    ulong m, q, mm;
    slong i;
    m = p - 1;

    for (i = 0; i < p1fac->num; i++)
    {
        q = p1fac->p[i];

        while (m % q == 0)
        {
            mm = m / q;
            if (n_powmod2_preinv(x, mm, p, pinv) != 1)
                break;
            m = mm;
        }
    }

    return m;
}

int nmod_poly_irreducible_binomial(nmod_poly_t res, ulong n)
{
    ulong q = res->mod.n;
    ulong e, a, qinv;
    n_factor_t nfac, q1fac;
    slong i;

    if (n == 0)
        return 0;

    if (n == 1)
    {
        nmod_poly_one(res);
        nmod_poly_set_coeff_ui(res, 1, 1);
        return 1;
    }

    /* There is an irreducible binomial x^n + a mod q iff
       p | q - 1 for each prime p | n and
       if 4 | n then 4 | q - 1. */
    if (n % 4 == 0)
        if (q % 4 != 1)
            return 0;
    n_factor_init(&nfac);
    n_factor(&nfac, n, 1);
    for (i = 0; i < nfac.num; i++)
        if ((q - 1) % nfac.p[i] != 0)
            return 0;

    nmod_poly_zero(res);
    nmod_poly_set_coeff_ui(res, n, 1);

    /* Let e = ord_q(-a). Then we have irreducibility iff
       (p divides e  but  p does not divide (q-1)/e)
       for each prime factor p of n. */
    n_factor_init(&q1fac);
    n_factor(&q1fac, q - 1, 1);
    qinv = n_preinvert_limb(q);

    for (a = 1; a < q; a++)
    {
        e = n_multiplicative_order(q - a, q, qinv, &q1fac);

        for (i = 0; i < nfac.num; i++)
            if (!(e % nfac.p[i] == 0 && ((q - 1) / e) % nfac.p[i] != 0))
                goto next;

        nmod_poly_set_coeff_ui(res, 0, a);
        return 1;

        next:
            continue;
    }

    flint_throw(FLINT_ERROR, "Failed to construct irreducible binomial (n = %wu, p = %wu)", n, q);
}



/* Todo: refine and extend this table */
static const unsigned int sieve_limit_tab[][10] = {
     /* p     1    2     3     4     5      6      7     8     */
    {   2,    1,   1,    1,    1,    24,    24,    256,  256,  0 },
    {   3,    1,   1,    1,    1,    128,   512,   512,  3072, 0 },
    {   5,    1,   1,    1,    48,   512,   1536,  3072, 0 },
    {   7,    1,   1,    1,    64,   1536,  3072,  0 },
    {  11,    1,   1,    64,   768,  2048,  0 },
    {  13,    1,   1,    48,   768,  0 },
    {  17,    1,   1,    384,  768,  0 },
    {  19,    1,   16,   96,   1536, 0 },
    {  23,    1,   16,   256,  2048, 0 },
    {  29,    1,   12,   512,  0 },
    {  31,    1,   12,   384,  0 },
    {  37,    1,   12,   384,  0 },
    {  41,    1,   24,   1024, 0 },
    {  43,    1,   32,   512,  0 },
    {  47,    1,   32,   768,  0 },
    {  53,    1,   48,   1024, 0 },
    {  59,    1,   32,   1024, 0 },
    {  61,    1,   64,   1024, 0 },
    {  67,    1,   64,   1536, 0 },
    {  71,    1,   64,   1536, 0 },
    {  73,    1,   64,   1536, 0 },
    {  79,    1,   64,   1536, 0 },
    {  83,    1,   64,   2048, 0 },
    {  89,    1,   64,   2048, 0 },
    {  97,    4,   128,  2048, 0 },
    { 101,    8,   128,  3072, 0 },
    { 103,    4,   96,   3072, 0 },
    { 107,    4,   96,   0 },
    { 127,    4,   128,  0 },
    { 179,    4,   192,  0 },
    { 269,    8,   256,  0 },
    { 283,    8,   384,  0 },
    { 431,    8,   512,  0 },
    { 487,    8,   768,  0 },
    { 653,    12,  768,  0 },
    { 727,    16,  1024, 0 },
    { 953,    16,  1536, 0 },
    { 1009,   24,  2048, 0 },
    { 1549,   32, 0 },
    { 2333,   48, 0 },
    { 3877,   64, 0 },
    { 7789,   96, 0 },
    { 14107,  128, 0 },
    { 28069,  192, 0 },
    { 35831,  256, 0 },
    { 82457,  384, 0 },
    { 121949, 768, 0 },
    { 180247, 1024, 0 },
    { 0 },
};

/* Return b such that we want to do trial division by all irreducible
   factors up to degree b. */
static ulong trinomial_sieve_limit(ulong p, ulong n)
{
    slong i;
    ulong ptab, b, limit;

    for (i = 0; ; i++)
    {
        ptab = sieve_limit_tab[i][0];
        if (ptab == 0)
            return 0;

        if (ptab >= p)
        {
            limit = 0;
            for (b = 1; sieve_limit_tab[i][b] != 0; b++)
            {
                if (n >= sieve_limit_tab[i][b])
                    limit = b;
            }

            return limit;
        }
    }
}


/* Return the product of all irreducible polynomial of degree r. */
/* Well-known fact: x^(p^r) - x contains all factors of degree dividing r,
   so we just need to divide out the smaller products. */
static void
nmod_poly_product_all_irreducibles_deg(nmod_poly_t res, ulong r)
{
    ulong k, d, step, p = res->mod.n;

    FLINT_ASSERT(r >= 1 && r <= 8);

    nmod_poly_zero(res);

    if (r == 1)
    {
        nmod_poly_set_coeff_ui(res, 1, 1);
        nmod_poly_set_coeff_ui(res, p, 1);
        return;
    }

    if (r == 4 || r == 6)
    {
        d = n_pow(p, r) - n_pow(p, 2);
        step = p * p - 1;
    }
    else if (r == 8)
    {
        d = n_pow(p, r) - n_pow(p, 4);
        step = n_pow(p, 4) - 1;
    }
    else
    {
        d = n_pow(p, r) - p;
        step = p - 1;
    }

    for (k = 0; k <= d; k += step)
        nmod_poly_set_coeff_ui(res, k, 1);

    if (r == 6)
    {
        nmod_poly_t t;
        nmod_poly_init(t, p);
        nmod_poly_product_all_irreducibles_deg(t, 3);
        nmod_poly_divexact(res, res, t);
        nmod_poly_clear(t);
    }
}

/* Hack: nmod_poly_gcd is not optimised for unbalanced GCD with the smaller
   operand sparse, so we deal with this here. */
static void
_nmod_poly_inplace_rem_sparse_monic(nn_ptr R,
        slong lenA, nn_srcptr Bcoeffs, const slong * Bexps, slong nzB, slong lenB, nmod_t mod)
{
    slong i, j, k;
    slong n = lenB - 1;
    ulong c;

    for (i = lenA - 1; i >= n; i--)
    {
        c = R[i];

        for (k = nzB - 2; k >= 0; k--)
        {
            j = Bexps[k];
            R[j + i - n] = nmod_sub(R[j + i - n], nmod_mul(c, Bcoeffs[k], mod), mod);
        }
    }
}

static void
nmod_poly_gcd_with_sparse(nmod_poly_t res, const nmod_poly_t A, const nmod_poly_t B,
        nn_srcptr Bcoeffs, const slong * Bexps, slong nzB)
{
    if (A->length <= B->length)
    {
        nmod_poly_gcd(res, A, B);
    }
    else
    {
        nn_ptr R;
        slong Rlen = B->length - 1;
        nmod_poly_t T;
        TMP_INIT;
        TMP_START;
        R = TMP_ALLOC(sizeof(ulong) * A->length);
        _nmod_vec_set(R, A->coeffs, A->length);
        _nmod_poly_inplace_rem_sparse_monic(R, A->length, Bcoeffs, Bexps, nzB, B->length, res->mod);
        NMOD_VEC_NORM(R, Rlen);
        T->coeffs = R;
        T->length = Rlen;
        T->alloc = Rlen;
        nmod_poly_gcd(res, T, B);
        TMP_END;
    }
}

/* Discriminant of x^n + ax^k + b  mod p > 2 (Swan's formula). */
static ulong
_nmod_poly_trinomial_discriminant(ulong n, ulong k, ulong a, ulong b, nmod_t mod)
{
    ulong d, n1, k1;
    ulong D1, D2, D3;

    d = n_gcd(n, k);
    n1 = n / d;
    k1 = k / d;

    D1 = nmod_pow_ui(b, k - 1, mod);
    D2 = nmod_mul(nmod_pow_ui(n, n1, mod), nmod_pow_ui(b, n1 - k1, mod), mod);
    D3 = nmod_mul(nmod_pow_ui(n - k, n1 - k1, mod), nmod_pow_ui(k, k1, mod), mod);
    D3 = nmod_mul(D3, nmod_pow_ui(a, n1, mod), mod);
    if (n1 % 2)
        D3 = nmod_neg(D3, mod);
    D2 = nmod_sub(D2, D3, mod);
    D2 = nmod_pow_ui(D2, d, mod);
    D1 = nmod_mul(D1, D2, mod);
    if ((n * (n - 1) / 2) % 2)
        D1 = nmod_neg(D1, mod);

    return D1;
}

/* Exceptional n for which no irreducible trinomial exists. These tables are
   not needed for correctness, but speed up the search. */

/* All exceptions for n <= 10000 */
static const uint16_t no_trinomial_mod3[] = {
    49, 57, 65, 68, 75, 98, 105, 123, 129, 130, 132, 149, 161, 175, 189, 197,
    207, 212, 213, 221, 223, 231, 233, 264, 267, 276, 281, 292, 297, 298, 303, 309,
    311, 319, 332, 343, 391, 394, 397, 401, 404, 405, 411, 415, 426, 435, 436, 437,
    439, 441, 442, 453, 459, 463, 489, 497, 513, 521, 524, 528, 534, 540, 543, 545,
    548, 551, 552, 556, 559, 561, 571, 572, 591, 597, 612, 613, 619, 622, 627, 628,
    633, 636, 641, 653, 657, 664, 671, 683, 686, 692, 699, 703, 708, 715, 723, 725,
    732, 737, 747, 748, 753, 761, 783, 794, 795, 799, 802, 804, 809, 822, 823, 825,
    836, 847, 851, 852, 855, 874, 878, 879, 881, 891, 893, 903, 905, 908, 909, 918,
    933, 951, 953, 957, 961, 965, 975, 979, 981, 989, 994, 996, 1003, 1005, 1011, 1029,
    1039, 1041, 1047, 1051, 1056, 1059, 1065, 1068, 1073, 1075, 1080, 1086, 1097, 1099, 1102, 1104,
    1109, 1111, 1125, 1127, 1129, 1135, 1140, 1142, 1147, 1153, 1164, 1193, 1195, 1203, 1221, 1224,
    1225, 1226, 1227, 1233, 1238, 1239, 1241, 1252, 1253, 1261, 1263, 1267, 1269, 1272, 1277, 1281,
    1297, 1308, 1316, 1342, 1345, 1356, 1357, 1363, 1365, 1366, 1385, 1389, 1393, 1411, 1416, 1419,
    1431, 1435, 1445, 1450, 1452, 1464, 1467, 1471, 1481, 1484, 1491, 1500, 1505, 1507, 1523, 1524,
    1533, 1555, 1565, 1566, 1572, 1581, 1587, 1589, 1603, 1605, 1608, 1611, 1617, 1623, 1628, 1639,
    1644, 1646, 1647, 1649, 1673, 1677, 1687, 1689, 1694, 1695, 1697, 1704, 1709, 1713, 1723, 1725,
    1731, 1733, 1747, 1748, 1749, 1762, 1785, 1786, 1788, 1805, 1810, 1829, 1831, 1836, 1851, 1869,
    1873, 1876, 1877, 1883, 1887, 1889, 1892, 1902, 1905, 1906, 1914, 1917, 1922, 1925, 1930, 1937,
    1957, 1958, 1967, 1971, 1975, 1978, 1992, 2004, 2005, 2007, 2010, 2019, 2023, 2031, 2047, 2057,
    2061, 2071, 2078, 2103, 2108, 2112, 2116, 2118, 2124, 2130, 2136, 2139, 2143, 2146, 2150, 2153,
    2155, 2160, 2172, 2181, 2189, 2198, 2203, 2204, 2205, 2208, 2213, 2217, 2218, 2220, 2222, 2229,
    2241, 2254, 2257, 2258, 2259, 2265, 2270, 2273, 2280, 2283, 2284, 2297, 2299, 2308, 2311, 2316,
    2323, 2328, 2331, 2333, 2339, 2348, 2357, 2367, 2369, 2375, 2386, 2390, 2395, 2405, 2413, 2417,
    2428, 2429, 2433, 2441, 2448, 2451, 2454, 2455, 2457, 2460, 2469, 2476, 2478, 2479, 2487, 2489,
    2501, 2511, 2517, 2526, 2527, 2534, 2538, 2541, 2544, 2554, 2555, 2557, 2559, 2561, 2571, 2591,
    2595, 2611, 2616, 2633, 2636, 2643, 2644, 2649, 2657, 2668, 2671, 2681, 2684, 2690, 2700, 2708,
    2711, 2712, 2715, 2733, 2745, 2755, 2769, 2770, 2773, 2779, 2781, 2786, 2787, 2791, 2801, 2803,
    2809, 2815, 2825, 2827, 2831, 2832, 2838, 2843, 2844, 2862, 2863, 2870, 2890, 2904, 2908, 2909,
    2921, 2923, 2924, 2928, 2935, 2942, 2955, 2959, 2961, 2962, 2965, 2993, 3000, 3003, 3010, 3014,
    3017, 3027, 3029, 3033, 3043, 3048, 3049, 3051, 3053, 3060, 3063, 3067, 3069, 3076, 3077, 3087,
    3091, 3093, 3116, 3117, 3119, 3123, 3124, 3125, 3127, 3132, 3140, 3144, 3151, 3165, 3167, 3171,
    3178, 3183, 3185, 3204, 3206, 3207, 3212, 3216, 3219, 3222, 3228, 3244, 3246, 3249, 3252, 3261,
    3269, 3281, 3288, 3291, 3292, 3297, 3308, 3329, 3332, 3346, 3348, 3357, 3367, 3369, 3387, 3396,
    3401, 3408, 3413, 3418, 3423, 3428, 3429, 3435, 3439, 3446, 3449, 3459, 3466, 3471, 3481, 3487,
    3507, 3508, 3540, 3556, 3567, 3572, 3576, 3585, 3593, 3596, 3597, 3603, 3605, 3609, 3610, 3615,
    3619, 3621, 3643, 3662, 3663, 3668, 3672, 3684, 3687, 3689, 3691, 3703, 3711, 3716, 3717, 3724,
    3737, 3746, 3749, 3751, 3754, 3756, 3764, 3765, 3766, 3774, 3777, 3778, 3780, 3787, 3791, 3797,
    3804, 3810, 3825, 3828, 3837, 3844, 3853, 3855, 3859, 3860, 3861, 3871, 3875, 3876, 3885, 3891,
    3909, 3914, 3919, 3929, 3931, 3932, 3934, 3941, 3945, 3950, 3953, 3955, 3963, 3972, 3977, 3981,
    3984, 3999, 4004, 4008, 4010, 4017, 4020, 4023, 4025, 4027, 4029, 4039, 4041, 4062, 4065, 4073,
    4119, 4122, 4123, 4131, 4133, 4135, 4137, 4140, 4142, 4143, 4156, 4157, 4161, 4171, 4197, 4216,
    4217, 4219, 4224, 4227, 4232, 4236, 4244, 4245, 4248, 4257, 4260, 4261, 4272, 4278, 4286, 4289,
    4293, 4305, 4308, 4313, 4315, 4320, 4327, 4344, 4362, 4363, 4380, 4385, 4389, 4396, 4416, 4426,
    4433, 4436, 4440, 4444, 4445, 4447, 4452, 4457, 4459, 4469, 4471, 4473, 4476, 4479, 4487, 4493,
    4495, 4499, 4503, 4508, 4519, 4523, 4530, 4531, 4540, 4545, 4546, 4547, 4548, 4560, 4565, 4569,
    4571, 4581, 4591, 4594, 4597, 4601, 4613, 4622, 4627, 4632, 4635, 4646, 4647, 4651, 4652, 4653,
    4656, 4663, 4665, 4666, 4668, 4678, 4684, 4697, 4707, 4711, 4714, 4721, 4738, 4747, 4750, 4753,
    4771, 4772, 4773, 4781, 4791, 4796, 4817, 4833, 4834, 4836, 4841, 4843, 4855, 4858, 4860, 4876,
    4882, 4896, 4899, 4901, 4902, 4905, 4907, 4908, 4910, 4920, 4923, 4935, 4938, 4949, 4956, 4971,
    4975, 4985, 4987, 4989, 4993, 4997, 5002, 5009, 5011, 5015, 5017, 5023, 5041, 5045, 5052, 5054,
    5055, 5059, 5061, 5068, 5069, 5076, 5079, 5083, 5084, 5087, 5088, 5093, 5097, 5108, 5122, 5132,
    5151, 5155, 5165, 5171, 5182, 5193, 5199, 5211, 5212, 5219, 5222, 5228, 5232, 5257, 5261, 5268,
    5286, 5297, 5301, 5303, 5313, 5314, 5321, 5333, 5342, 5345, 5349, 5357, 5361, 5362, 5369, 5372,
    5373, 5391, 5393, 5400, 5416, 5417, 5422, 5424, 5430, 5435, 5453, 5457, 5465, 5466, 5481, 5487,
    5491, 5492, 5509, 5510, 5513, 5515, 5524, 5525, 5538, 5539, 5541, 5546, 5553, 5561, 5565, 5567,
    5571, 5575, 5582, 5588, 5601, 5602, 5606, 5612, 5618, 5621, 5630, 5635, 5652, 5655, 5659, 5662,
    5664, 5673, 5676, 5679, 5683, 5688, 5697, 5703, 5707, 5713, 5721, 5724, 5726, 5733, 5739, 5755,
    5763, 5764, 5765, 5780, 5787, 5791, 5799, 5801, 5804, 5808, 5818, 5823, 5825, 5828, 5837, 5839,
    5842, 5844, 5846, 5852, 5856, 5870, 5876, 5892, 5895, 5905, 5910, 5918, 5919, 5923, 5925, 5929,
    5930, 5931, 5933, 5935, 5941, 5945, 5955, 5957, 5967, 5969, 5989, 5995, 6000, 6003, 6005, 6029,
    6031, 6034, 6036, 6039, 6052, 6054, 6060, 6066, 6079, 6096, 6099, 6101, 6102, 6106, 6111, 6116,
    6120, 6134, 6137, 6141, 6161, 6163, 6165, 6177, 6180, 6186, 6191, 6207, 6211, 6212, 6223, 6243,
    6250, 6251, 6254, 6255, 6259, 6264, 6265, 6283, 6287, 6288, 6301, 6317, 6329, 6330, 6334, 6356,
    6361, 6363, 6365, 6370, 6372, 6377, 6379, 6393, 6408, 6412, 6414, 6425, 6427, 6432, 6444, 6456,
    6463, 6465, 6473, 6476, 6477, 6483, 6485, 6489, 6492, 6493, 6495, 6504, 6508, 6524, 6531, 6547,
    6548, 6551, 6559, 6567, 6576, 6579, 6585, 6587, 6603, 6609, 6619, 6621, 6629, 6652, 6675, 6677,
    6687, 6689, 6696, 6714, 6716, 6734, 6738, 6739, 6743, 6772, 6775, 6783, 6792, 6795, 6802, 6811,
    6816, 6826, 6847, 6856, 6858, 6867, 6868, 6893, 6898, 6905, 6909, 6919, 6937, 6948, 6949, 6955,
    6957, 6967, 6974, 6999, 7005, 7011, 7025, 7037, 7041, 7067, 7075, 7079, 7080, 7087, 7093, 7099,
    7111, 7113, 7121, 7123, 7141, 7149, 7152, 7161, 7167, 7181, 7183, 7186, 7197, 7203, 7207, 7219,
    7221, 7229, 7230, 7231, 7235, 7238, 7241, 7252, 7263, 7267, 7275, 7286, 7293, 7301, 7305, 7308,
    7313, 7339, 7344, 7345, 7348, 7349, 7359, 7363, 7364, 7368, 7373, 7380, 7382, 7383, 7385, 7395,
    7406, 7407, 7412, 7429, 7433, 7435, 7452, 7467, 7471, 7474, 7483, 7491, 7493, 7495, 7500, 7512,
    7530, 7532, 7533, 7548, 7553, 7556, 7559, 7560, 7567, 7575, 7594, 7597, 7604, 7608, 7613, 7617,
    7620, 7625, 7628, 7629, 7635, 7637, 7639, 7641, 7649, 7653, 7656, 7663, 7669, 7673, 7683, 7687,
    7688, 7692, 7697, 7706, 7707, 7710, 7713, 7717, 7727, 7731, 7732, 7735, 7737, 7742, 7745, 7750,
    7752, 7755, 7757, 7761, 7769, 7775, 7779, 7797, 7803, 7818, 7819, 7827, 7833, 7838, 7853, 7855,
    7857, 7858, 7862, 7867, 7876, 7881, 7882, 7889, 7890, 7895, 7903, 7924, 7926, 7944, 7945, 7948,
    7951, 7959, 7965, 7968, 7971, 7972, 7973, 7980, 7981, 7993, 7997, 7998, 8004, 8009, 8016, 8021,
    8023, 8027, 8033, 8037, 8040, 8050, 8052, 8057, 8059, 8068, 8075, 8078, 8079, 8081, 8084, 8091,
    8095, 8097, 8100, 8103, 8109, 8119, 8124, 8129, 8132, 8143, 8146, 8155, 8156, 8167, 8179, 8199,
    8201, 8203, 8205, 8209, 8215, 8220, 8238, 8239, 8244, 8246, 8251, 8253, 8259, 8261, 8262, 8265,
    8266, 8268, 8270, 8271, 8273, 8276, 8280, 8284, 8289, 8321, 8323, 8333, 8342, 8343, 8356, 8371,
    8388, 8396, 8405, 8417, 8429, 8434, 8436, 8438, 8441, 8448, 8451, 8453, 8468, 8469, 8472, 8475,
    8496, 8513, 8514, 8520, 8522, 8527, 8532, 8533, 8544, 8549, 8553, 8556, 8557, 8571, 8573, 8575,
    8586, 8589, 8597, 8609, 8610, 8616, 8626, 8630, 8640, 8643, 8647, 8649, 8657, 8668, 8673, 8681,
    8685, 8688, 8691, 8693, 8695, 8700, 8709, 8717, 8724, 8726, 8733, 8743, 8757, 8760, 8764, 8767,
    8769, 8772, 8779, 8781, 8789, 8799, 8801, 8803, 8805, 8807, 8823, 8825, 8827, 8832, 8835, 8844,
    8845, 8847, 8849, 8852, 8857, 8865, 8880, 8889, 8890, 8894, 8904, 8911, 8913, 8914, 8918, 8923,
    8925, 8927, 8935, 8937, 8940, 8942, 8943, 8949, 8952, 8958, 8961, 8971, 8974, 8985, 8986, 8990,
    8993, 8997, 8998, 9004, 9006, 9012, 9019, 9021, 9033, 9039, 9041, 9046, 9057, 9060, 9065, 9069,
    9092, 9094, 9096, 9099, 9103, 9105, 9116, 9117, 9120, 9124, 9127, 9137, 9138, 9142, 9145, 9147,
    9148, 9153, 9156, 9163, 9180, 9182, 9187, 9189, 9194, 9213, 9223, 9225, 9235, 9239, 9254, 9261,
    9264, 9270, 9285, 9294, 9302, 9305, 9309, 9312, 9317, 9329, 9330, 9335, 9336, 9340, 9341, 9351,
    9353, 9367, 9372, 9377, 9387, 9389, 9394, 9396, 9399, 9405, 9407, 9411, 9413, 9414, 9420, 9422,
    9433, 9435, 9436, 9444, 9447, 9452, 9471, 9476, 9492, 9497, 9499, 9506, 9507, 9509, 9513, 9524,
    9529, 9540, 9542, 9544, 9545, 9546, 9549, 9559, 9561, 9562, 9569, 9571, 9579, 9582, 9587, 9595,
    9596, 9597, 9605, 9612, 9615, 9617, 9619, 9621, 9637, 9639, 9644, 9653, 9655, 9665, 9666, 9669,
    9672, 9677, 9682, 9684, 9686, 9697, 9701, 9705, 9708, 9710, 9711, 9716, 9720, 9737, 9749, 9759,
    9773, 9779, 9780, 9783, 9792, 9796, 9798, 9802, 9804, 9807, 9809, 9813, 9816, 9819, 9831, 9840,
    9847, 9849, 9855, 9875, 9876, 9877, 9879, 9881, 9893, 9901, 9907, 9911, 9912, 9921, 9924, 9931,
    9933, 9943, 9950, 9951, 9953, 9957, 9970, 9975, 9977, 9978, 9979, 9980, 0
};

/* All exceptions for n <= 5000 */
static const uint16_t no_trinomial_mod5[] = {
    35, 70, 123, 125, 140, 181, 191, 209, 213, 219, 237, 249, 250, 253, 265, 273,
    280, 285, 307, 319, 345, 351, 362, 375, 382, 407, 413, 415, 418, 421, 429, 433,
    441, 445, 457, 461, 498, 500, 519, 530, 539, 547, 560, 570, 581, 638, 663, 683,
    690, 699, 723, 724, 725, 735, 750, 753, 765, 779, 787, 795, 803, 815, 817, 830,
    833, 835, 859, 866, 885, 890, 915, 917, 927, 957, 987, 993, 1000, 1015, 1037, 1038,
    1043, 1055, 1060, 1087, 1095, 1113, 1120, 1121, 1137, 1140, 1159, 1161, 1173, 1183, 1185, 1203,
    1249, 1257, 1261, 1267, 1275, 1276, 1293, 1299, 1335, 1351, 1357, 1359, 1365, 1380, 1381, 1387,
    1407, 1427, 1443, 1449, 1450, 1470, 1473, 1479, 1500, 1501, 1503, 1515, 1517, 1530, 1537, 1557,
    1558, 1575, 1590, 1606, 1619, 1621, 1630, 1631, 1653, 1660, 1665, 1670, 1671, 1672, 1683, 1687,
    1691, 1718, 1719, 1727, 1732, 1749, 1769, 1770, 1780, 1787, 1827, 1830, 1831, 1833, 1834, 1837,
    1855, 1869, 1873, 1877, 1887, 1899, 1901, 1911, 1914, 1929, 1959, 1967, 1977, 1986, 2000, 2025,
    2030, 2037, 2041, 2076, 2087, 2091, 2097, 2110, 2120, 2141, 2151, 2157, 2163, 2173, 2187, 2190,
    2240, 2242, 2251, 2257, 2261, 2280, 2319, 2327, 2339, 2341, 2345, 2353, 2359, 2367, 2370, 2379,
    2383, 2407, 2421, 2433, 2447, 2459, 2463, 2485, 2491, 2517, 2519, 2522, 2535, 2541, 2550, 2552,
    2575, 2583, 2595, 2601, 2607, 2621, 2631, 2645, 2653, 2670, 2673, 2683, 2702, 2705, 2717, 2730,
    2745, 2755, 2757, 2760, 2762, 2763, 2775, 2809, 2843, 2847, 2867, 2871, 2877, 2879, 2885, 2891,
    2895, 2899, 2900, 2917, 2940, 2967, 2991, 3000, 3002, 3003, 3027, 3030, 3035, 3049, 3051, 3053,
    3060, 3074, 3095, 3103, 3125, 3131, 3150, 3165, 3173, 3180, 3203, 3212, 3213, 3242, 3260, 3275,
    3295, 3320, 3327, 3330, 3340, 3344, 3361, 3366, 3377, 3387, 3395, 3403, 3435, 3439, 3447, 3455,
    3463, 3464, 3469, 3498, 3525, 3527, 3540, 3553, 3555, 3560, 3574, 3645, 3660, 3675, 3687, 3710,
    3715, 3723, 3727, 3746, 3759, 3795, 3802, 3803, 3807, 3823, 3831, 3877, 3885, 3889, 3895, 3901,
    3907, 3911, 3925, 3943, 3945, 3951, 3953, 3969, 3972, 3984, 4000, 4013, 4029, 4041, 4050, 4060,
    4077, 4079, 4095, 4111, 4137, 4147, 4151, 4163, 4165, 4174, 4211, 4220, 4227, 4240, 4263, 4282,
    4293, 4305, 4325, 4346, 4367, 4380, 4391, 4409, 4425, 4443, 4465, 4479, 4480, 4484, 4487, 4489,
    4497, 4503, 4560, 4563, 4573, 4605, 4635, 4647, 4655, 4663, 4682, 4687, 4690, 4691, 4697, 4706,
    4740, 4755, 4767, 4803, 4807, 4825, 4839, 4851, 4879, 4894, 4901, 4905, 4911, 4915, 4918, 4947,
    4951, 4953, 4970, 4981, 4995, 4997, 0
};

/* All exceptions for n <= 5000 */
static const uint16_t no_trinomial_mod7[] = {
    124, 163, 268, 301, 359, 364, 372, 385, 389, 407, 431, 476, 489, 527, 563, 601,
    628, 649, 728, 829, 835, 853, 862, 952, 995, 1031, 1063, 1077, 1079, 1115, 1259, 1285,
    1343, 1355, 1372, 1384, 1444, 1456, 1467, 1477, 1669, 1703, 1705, 1706, 1735, 1777, 1781, 1803,
    1825, 1828, 1884, 1889, 1904, 1975, 2003, 2044, 2057, 2083, 2093, 2105, 2158, 2189, 2281, 2339,
    2395, 2505, 2507, 2513, 2564, 2668, 2744, 2755, 2767, 2773, 2807, 2812, 2833, 2863, 2867, 2912,
    2985, 3031, 3091, 3093, 3157, 3203, 3263, 3271, 3289, 3332, 3345, 3406, 3463, 3487, 3515, 3563,
    3650, 3656, 3716, 3719, 3763, 3797, 3808, 3836, 3875, 3892, 3907, 3916, 3971, 4045, 4055, 4088,
    4116, 4166, 4204, 4210, 4249, 4321, 4340, 4391, 4435, 4607, 4673, 4681, 4687, 4708, 4748, 4763,
    4780, 4781, 4961, 0
};

/* We skip the mod 2 table. Although the exceptions are quite dense,
   the performance impact is barely noticeable (most time will be spent
   searching pentanomials). It would be more worthwhile to just include
   a table of the irreducible polynomials themselves over GF(2). */

/* Exceptions are also known mod 11, 13 and 17, but they are very sparse
   and hence not really worth it. */

/*
static const uint16_t no_trinomial_mod2[] = {
    8, 13, 16, 19, 24, 26, 27, 32, 37, 38, 40, 43, 45, 48, 50, 51,
    53, 56, 59, 61, 64, 67, 69, 70, 72, 75, 77, 78, 80, 82, 83, 85,
    // ...
};

// All exceptions for n <= 5000
static const uint16_t no_trinomial_mod11[] = {
    219, 291, 467, 599, 761, 909, 1137, 1243, 2859, 2873, 2883, 3047, 3531, 3651, 3932, 4139,
    4207, 4247, 4383, 4419, 4508, 4545, 4551, 4564, 4741, 0
};

// All exceptions for n <= 2000
static const uint16_t no_trinomial_mod13[] = {
    833, 0
};

// All exceptions for n <= 2000
static const uint16_t no_trinomial_mod17[] = {
    231, 375, 675, 717, 1941, 0
}; */

int nmod_poly_irreducible_trinomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong a, b, c, k, l;
    ulong a2, b2, e2n, cn, ck, cn2, ck2, D;
    slong i;
    nn_ptr invctab = NULL;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[3];
    slong Bexps[3];
    ulong sieve_limit;

    if (n <= 1)
        return 0;

    if (p == 2 && (n % 8 == 0))
        return 0;

    /* Ciet, Quisquater and Sica, A Short Note on Irreducible Trinomials in Binary Fields */
    if (p == 2 && ((n % 24 == 13) || (n % 24 == 19)) && n_is_prime(n))
        return 0;

    if (p == 3 || p == 5 || p == 7)
    {
        const uint16_t * tab;

        if (p == 3)      tab = no_trinomial_mod3;
        else if (p == 5) tab = no_trinomial_mod5;
        else             tab = no_trinomial_mod7;

        for (i = 0; tab[i] != 0 && tab[i] <= n; i++)
            if (tab[i] == n)
                return 0;
    }

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        invctab = _nmod_vec_init(p);
        for (c = 1; c < p; c++)
            invctab[c] = nmod_inv(c, mod);

        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

#define ALREADY_CHECKED(a,b,a2,b2) ((a2) < (a) || ((a2) == (a) && (b2) < (b)))

    /* Iterate over x^n + a x^k + b */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            /* Mirror symmetry: f(x) is irreducible iff x^n f(1/x) is irreducible */
            for (k = 1; k * 2 <= n; k++)
            {
                if (p > 2)
                {
                    /* See if we have already checked f(-x). */
                    /* n even, k even   x^n + ax^k + b ->  x^n + ax^k + b */
                    /* n even, k odd    x^n + ax^k + b ->  x^n - ax^k + b */
                    /* n odd,  k even   x^n + ax^k + b ->  x^n - ax^k - b */
                    /* n odd,  k odd    x^n + ax^k + b ->  x^n + ax^k - b */
                    b2 = (n % 2) ? p - b : b;
                    a2 = ((n ^ k) % 2) ? p - a : a;
                    if (ALREADY_CHECKED(a, b, a2, b2))
                        goto next_trinomial;
                }

                /* Swan-Stickelberger */
                if (p == 2)
                {
                    /* Perfect square */
                    if (n % 2 == 0 && k % 2 == 0)
                        goto next_trinomial;

                    if (n % 2 == k % 2)
                        l = n - k;
                    else
                        l = k;

                    if (n % 2 == 0 && l % 2 == 1 && n != 2 * l && ((n * l) / 2) % 4 <= 1)
                        goto next_trinomial;
                    if (n % 2 == 1 && l % 2 == 0 && ((2*n) % l != 0) && (n % 8 == 3 || n % 8 == 5))
                        goto next_trinomial;
                    if (n % 2 == 1 && l % 2 == 0 && ((2*n) % l == 0) && (n % 8 == 1 || n % 8 == 7))
                        goto next_trinomial;
                }
                else
                {
                    D = _nmod_poly_trinomial_discriminant(n, k, a, b, mod);

                    if (nmod_pow_ui(D, (p - 1) / 2, mod) == ((n % 2 == 1) ? p - 1 : 1))
                        goto next_trinomial;
                }

                if (want_pruning && p > 2)
                {
                    for (c = 1; 2 * c < p; c++)
                    {
                        ck = nmod_pow_ui(c, k, mod);
                        cn = nmod_mul(ck, nmod_pow_ui(c, n - k, mod), mod);

                        if (c >= 2)
                        {
                            /* See if we have already checked f(c*x). */
                            e2n = invctab[cn];
                            b2 = nmod_mul(b, e2n, mod);
                            a2 = nmod_mul(nmod_mul(a, ck, mod), e2n, mod);
                            if (ALREADY_CHECKED(a, b, a2, b2))
                                goto next_trinomial;

                            /* See if we have already checked f(-c*x). */
                            b2 = (n % 2) ? nmod_neg(b2, mod) : b2;
                            a2 = ((n ^ k) % 2) ? nmod_neg(a2, mod) : a2;
                            if (ALREADY_CHECKED(a, b, a2, b2))
                                goto next_trinomial;
                        }

                        /* Trial division: check if c^n + a c^k + b = 0. */
                        if (nmod_add(cn, nmod_mul(a, ck, mod), mod) == p - b)
                            goto next_trinomial;

                        /* Trial division: check if (-c)^n + a (-c)^k + b = 0. */
                        cn2 = (n % 2) ? nmod_neg(cn, mod) : cn;
                        ck2 = (k % 2) ? nmod_neg(ck, mod) : ck;
                        if (nmod_add(cn2, nmod_mul(a, ck2, mod), mod) == p - b)
                            goto next_trinomial;
                    }
                }

                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, 0, b);
                nmod_poly_set_coeff_ui(res, k, a);

                for (i = 0; i < num_sieve; i++)
                {
                    Bcoeffs[0] = b;
                    Bcoeffs[1] = a;
                    Bcoeffs[2] = 1;
                    Bexps[0] = 0;
                    Bexps[1] = k;
                    Bexps[2] = n;

                    nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 3);
                    if (!nmod_poly_is_one(h))
                        goto next_trinomial;
                }

                if (nmod_poly_is_irreducible_ddf(res))
                {
                    found = 1;
                    goto cleanup;
                }

                next_trinomial:
                    continue;
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        _nmod_vec_clear(invctab);
        for (i = 0; i < num_sieve; i++)
            nmod_poly_clear(sieve[i]);
        nmod_poly_clear(h);
    }

    return found;
}

int nmod_poly_irreducible_tetranomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong a, b, c, d, k, l;
    ulong dk, dl, dn, dk2, dl2, dn2, s;
    slong i;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[4];
    slong Bexps[4];
    ulong sieve_limit;

    if (n <= 2 || p == 2)
        return 0;

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

    /* Iterate over all x^n + a x^k + b x^l + c */
    for (a = 1; a < p; a++)
    {
        for (b = 1; b < p; b++)
        {
            for (c = 1; c < p; c++)
            {
                for (k = 2; k < n; k++)
                {
                    for (l = 1; l < k; l++)
                    {
                        if (want_pruning)
                        {
                            for (d = 1; 2 * d < p; d++)
                            {
                                dk = nmod_pow_ui(d, k, mod);
                                dl = nmod_pow_ui(d, l, mod);
                                dn = nmod_pow_ui(d, n, mod);

                                /* Trial division: check if d^n + a d^k + b d^l + c = 0. */
                                s = nmod_add(dn, nmod_mul(a, dk, mod), mod);
                                s = nmod_add(s, nmod_mul(b, dl, mod), mod);
                                if (s == p - c)
                                    goto next_tetranomial;

                                dn2 = (n % 2) ? nmod_neg(dn, mod) : dn;
                                dk2 = (k % 2) ? nmod_neg(dk, mod) : dk;
                                dl2 = (l % 2) ? nmod_neg(dl, mod) : dl;
                                /* Trial division: check if (-d)^n + a (-d)^k + b (-d)^l + c = 0. */
                                s = nmod_add(dn2, nmod_mul(a, dk2, mod), mod);
                                s = nmod_add(s, nmod_mul(b, dl2, mod), mod);
                                if (s == p - c)
                                    goto next_tetranomial;
                            }
                        }

                        nmod_poly_zero(res);
                        nmod_poly_set_coeff_ui(res, n, 1);
                        nmod_poly_set_coeff_ui(res, k, a);
                        nmod_poly_set_coeff_ui(res, l, b);
                        nmod_poly_set_coeff_ui(res, 0, c);

                        for (i = 0; i < num_sieve; i++)
                        {
                            Bcoeffs[0] = c;
                            Bcoeffs[1] = b;
                            Bcoeffs[2] = a;
                            Bcoeffs[3] = 1;
                            Bexps[0] = 0;
                            Bexps[1] = l;
                            Bexps[2] = k;
                            Bexps[3] = n;

                            nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 4);
                            if (!nmod_poly_is_one(h))
                                goto next_tetranomial;
                        }

                        {
                            ulong D, D1, D2;
                            D = nmod_poly_discriminant(res);
                            D1 = nmod_pow_ui(D, (p - 1) / 2, mod);
                            D2 = (n % 2 == 1) ? p - 1 : 1; 
                            if (D1 == D2)
                                goto next_tetranomial;
                        }

                        if (nmod_poly_is_irreducible_ddf(res))
                        {
                            found = 1;
                            goto cleanup;
                        }

                        next_tetranomial:
                            continue;
                    }
                }
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        for (i = 0; i < num_sieve; i++)
            nmod_poly_clear(sieve[i]);
        nmod_poly_clear(h);
    }

    return found;
}

/* Only intended for GF(2), so we only bother with coefficients 1 */
int nmod_poly_irreducible_pentanomial(nmod_poly_t res, ulong n)
{
    ulong p = res->mod.n;
    nmod_t mod = res->mod;
    ulong k, l, m, a;
    slong i;
    int want_pruning;
    int found = 0;
    nmod_poly_t sieve[10], h;
    slong num_sieve = 0;
    ulong Bcoeffs[5];
    slong Bexps[5];
    ulong sieve_limit;

    if (n <= 3)
        return 0;

    sieve_limit = trinomial_sieve_limit(p, n);
    want_pruning = (sieve_limit > 0);

    if (want_pruning)
    {
        nmod_poly_init_mod(h, mod);

        for (a = 2; a <= sieve_limit && a < n; a++)
        {
            nmod_poly_init_mod(sieve[num_sieve], mod);
            nmod_poly_product_all_irreducibles_deg(sieve[num_sieve], a);
            num_sieve++;
        }
    }

    Bcoeffs[0] = 1;
    Bcoeffs[1] = 1;
    Bcoeffs[2] = 1;
    Bcoeffs[3] = 1;
    Bcoeffs[4] = 1;

    /* Iterate over all x^n + x^k + x^l + x^m + 1 */
    for (k = 3; k < n; k++)
    {
        for (l = 2; l < k; l++)
        {
            for (m = 1; m < l; m++)
            {
                /* Perfect square */
                if (p == 2 && n % 2 == 0 && k % 2 == 0 && l % 2 == 0 && m % 2 == 0)
                    goto next_pentanomial;

                nmod_poly_zero(res);
                nmod_poly_set_coeff_ui(res, n, 1);
                nmod_poly_set_coeff_ui(res, k, 1);
                nmod_poly_set_coeff_ui(res, l, 1);
                nmod_poly_set_coeff_ui(res, m, 1);
                nmod_poly_set_coeff_ui(res, 0, 1);

                for (i = 0; i < num_sieve; i++)
                {
                    Bexps[0] = 0;
                    Bexps[1] = m;
                    Bexps[2] = l;
                    Bexps[3] = k;
                    Bexps[4] = n;

                    nmod_poly_gcd_with_sparse(h, sieve[i], res, Bcoeffs, Bexps, 5);
                    if (!nmod_poly_is_one(h))
                        goto next_pentanomial;
                }

                if (nmod_poly_is_irreducible_ddf(res))
                {
                    found = 1;
                    goto cleanup;
                }

                next_pentanomial:
                    continue;
            }
        }
    }

cleanup:
    if (want_pruning)
    {
        for (i = 0; i < num_sieve; i++)
            nmod_poly_clear(sieve[i]);
        nmod_poly_clear(h);
    }

    return found;
}

void
nmod_poly_minimal_irreducible(nmod_poly_t res, ulong n)
{
    FLINT_ASSERT(n != 0);

    if (n == 1)
    {
        nmod_poly_zero(res);
        nmod_poly_set_coeff_ui(res, 1, 1);
        return;
    }

    if (nmod_poly_irreducible_binomial(res, n))
        return;

    if (nmod_poly_irreducible_trinomial(res, n))
        return;

    if (nmod_poly_irreducible_tetranomial(res, n))
        return;

    if (nmod_poly_irreducible_pentanomial(res, n))
        return;

    /* The conjecture is that a tetranomial (p != 2) or pentanomial (p = 2)
       always works. */
    flint_throw(FLINT_ERROR, "Failed to construct minimal irreducible polynomial (n = %wu, p = %wu)", n, res->mod.n);
}

