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

    Copyright (C) 2011 Fredrik Johansson

******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "flint.h"
#include "ulong_extras.h"
#include "fmpz.h"
#include "arith.h"
#include "mpfr.h"
#include "math.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

#define MAX_32BIT 58

static const int lookup_table[MAX_32BIT][28] =
{
 {-1, 1}, {1, 1}, {1, 2}, {0, 1}, {-1, 2, 4}, {-1, 2}, 
 {-1, -4, 4, 8}, {-1, 0, 2}, {1, -6, 0, 8}, {-1, -2, 4}, 
 {1, 6, -12, -32, 16, 32}, {-3, 0, 4}, {-1, 6, 24, -32, -80, 32, 64}, 
 {1, -4, -4, 8}, {1, 8, -16, -8, 16}, {1, 0, -8, 0, 8}, 
 {1, -8, -40, 80, 240, -192, -448, 128, 256}, {-1, -6, 0, 8}, 
 {1, 10, -40, -160, 240, 672, -448, -1024, 256, 512}, {5, 0, -20, 0, 16}, 
 {1, -16, 32, 48, -96, -32, 64}, {-1, 6, 12, -32, -16, 32}, 
 {-1, -12, 60, 280, -560, -1792, 1792, 4608, -2304, -5120, 1024, 2048}, 
 {1, 0, -16, 0, 16}, {-1, 10, 100, -40, -800, 32, 2240, 0, -2560, 0, 
  1024}, {-1, -6, 24, 32, -80, -32, 64}, 
 {1, 18, 0, -240, 0, 864, 0, -1152, 0, 512}, {-7, 0, 56, 0, -112, 0, 64}, 
 {-1, 14, 112, -448, -2016, 4032, 13440, -15360, -42240, 28160, 67584, 
  -24576, -53248, 8192, 16384}, {1, -8, -16, 8, 16}, 
 {-1, -16, 112, 672, -2016, -8064, 13440, 42240, -42240, -112640, 67584, 
  159744, -53248, -114688, 16384, 32768}, 
 {1, 0, -32, 0, 160, 0, -256, 0, 128}, 
 {1, -24, 48, 344, -688, -1088, 2176, 1280, -2560, -512, 1024}, 
 {1, 8, -40, -80, 240, 192, -448, -128, 256}, 
 {1, 16, -160, -368, 1760, 2272, -7232, -5504, 13824, 5632, -12288, 
  -2048, 4096}, {-3, 0, 36, 0, -96, 0, 64}, 
 {-1, 18, 180, -960, -5280, 14784, 59136, -101376, -329472, 366080, 
  1025024, -745472, -1863680, 860160, 1966080, -524288, -1114112, 131072, 
  262144}, {-1, 10, 40, -160, -240, 672, 448, -1024, -256, 512}, 
 {1, 24, -48, -632, 1264, 3296, -6592, -6784, 13568, 6144, -12288, -2048, 
  4096}, {1, 0, -48, 0, 304, 0, -512, 0, 256}, 
 {1, -20, -220, 1320, 7920, -25344, -109824, 219648, 768768, -1025024, 
  -3075072, 2795520, 7454720, -4587520, -11141120, 4456448, 10027008, 
  -2359296, -4980736, 524288, 1048576}, {1, 16, 32, -48, -96, 32, 64}, 
 {1, 22, -220, -1760, 7920, 41184, -109824, -439296, 768768, 2562560, 
  -3075072, -8945664, 7454720, 19496960, -11141120, -26738688, 10027008, 
  22413312, -4980736, -10485760, 1048576, 2097152}, 
 {-11, 0, 220, 0, -1232, 0, 2816, 0, -2816, 0, 1024}, 
 {1, -24, -144, 248, 1680, -864, -7168, 1152, 13824, -512, -12288, 0, 
  4096}, {1, -12, -60, 280, 560, -1792, -1792, 4608, 2304, -5120, -1024, 
  2048}, {-1, -24, 264, 2288, -11440, -64064, 192192, 823680, -1647360, 
  -5857280, 8200192, 25346048, -25346048, -70189056, 50135040, 127008768, 
  -63504384, -149422080, 49807360, 110100480, -22020096, -46137344, 
  4194304, 8388608}, {1, 0, -64, 0, 320, 0, -512, 0, 256}, 
 {-1, 28, 196, -2968, -3136, 66304, 18816, -658816, -53760, 3587584, 
  78848, -11741184, -57344, 24084480, 16384, -31195136, 0, 24772608, 0, 
  -11010048, 0, 2097152}, {-1, -10, 100, 40, -800, -32, 2240, 0, -2560, 
  0, 1024}, {1, 32, -64, -1504, 3008, 16832, -33664, -76288, 152576, 
  173568, -347136, -210944, 421888, 131072, -262144, -32768, 65536}, 
 {13, 0, -364, 0, 2912, 0, -9984, 0, 16640, 0, -13312, 0, 4096}, 
 {-1, 26, 364, -2912, -21840, 96096, 512512, -1464320, -6223360, 
  12446720, 44808192, -65175552, -206389248, 222265344, 635043840, 
  -508035072, -1333592064, 784465920, 1917583360, -807403520, 
  -1857028096, 530579456, 1157627904, -201326592, -419430400, 33554432, 
  67108864}, {-1, 18, 0, -240, 0, 864, 0, -1152, 0, 512}, 
 {1, 24, -432, -1208, 15216, 28064, -185024, -263424, 1149184, 1250304, 
  -4177920, -3356672, 9375744, 5324800, -13123584, -4947968, 11141120, 
  2490368, -5242880, -524288, 1048576}, 
 {1, 0, -96, 0, 1376, 0, -6656, 0, 13568, 0, -12288, 0, 4096}, 
 {1, -40, 80, 2120, -4240, -31648, 63296, 194432, -388864, -613376, 
  1226752, 1087488, -2174976, -1097728, 2195456, 589824, -1179648, 
  -131072, 262144}, {-1, -14, 112, 448, -2016, -4032, 13440, 15360, 
  -42240, -28160, 67584, 24576, -53248, -8192, 16384}
};

/* The coefficients in 2^d * \prod_{i=1}^d (x - cos(a_i)) are
   easily bounded using the binomial theorem. */
static len_t
magnitude_bound(len_t d)
{
    len_t res;
    fmpz_t t;
    fmpz_init(t);
    fmpz_bin_uiui(t, d, d / 2);
    res = fmpz_bits(t);
    fmpz_clear(t);
    return FLINT_ABS(res) + d;
}

static void
fmpz_mul_or_div_2exp(fmpz_t x, fmpz_t y, len_t s)
{
    if (s >= 0)
        fmpz_mul_2exp(x, y, s);
    else
        fmpz_fdiv_q_2exp(x, y, -s);
}


/* Balanced product of linear factors (x+alpha_i) using
   fixed-point arithmetic with prec bits */
static void
balanced_product(fmpz * c, fmpz * alpha, len_t len, len_t prec)
{
    if (len == 1)
    {
        fmpz_one(c + 1);
        fmpz_mul_2exp(c + 1, c + 1, prec);
        fmpz_set(c, alpha);
    }
    else if (len == 2)
    {
        fmpz_mul(c, alpha, alpha + 1);
        fmpz_fdiv_q_2exp(c, c, prec);
        fmpz_add(c + 1, alpha, alpha + 1);
        fmpz_one(c + 2);
        fmpz_mul_2exp(c + 2, c + 2, prec);
    }
    else
    {
        fmpz *L, *R;
        len_t i, m;

        m = len / 2;
        L = _fmpz_vec_init(len + 2);
        R = L + m + 1;

        balanced_product(L, alpha, m, prec);
        balanced_product(R, alpha + m, len - m, prec);
        _fmpz_poly_mul(c, R, len - m + 1, L, m + 1);

        for (i = 0; i < len + 1; i++)
            fmpz_fdiv_q_2exp(c + i, c + i, prec);

        _fmpz_vec_clear(L, len + 2);
    }
}

void
_arith_cos_minpoly(fmpz * coeffs, len_t d, ulong n)
{
    len_t i, j, prec, exp;
    fmpz * alpha;
    fmpz_t half;
    mpfr_t t, u;

    if (n <= MAX_32BIT)
    {
        for (i = 0; i <= d; i++)
            fmpz_set_si(coeffs + i, lookup_table[n - 1][i]);
        return;
    }

    /* Direct formula for odd primes > 3 */
    if (n_is_prime(n))
    {
        len_t s = (n - 1) / 2;

        switch (s % 4)
        {
            case 0:
                fmpz_set_si(coeffs, 1L);
                fmpz_set_si(coeffs + 1, -s);
                break;
            case 1:
                fmpz_set_si(coeffs, 1L);
                fmpz_set_si(coeffs + 1, s + 1);
                break;
            case 2:
                fmpz_set_si(coeffs, -1L);
                fmpz_set_si(coeffs + 1, s);
                break;
            case 3:
                fmpz_set_si(coeffs, -1L);
                fmpz_set_si(coeffs + 1, -s - 1);
                break;
        }

        for (i = 2; i <= s; i++)
        {
            len_t b = (s - i) % 2;
            fmpz_mul2_uiui(coeffs + i, coeffs + i - 2, s+i-b, s+2-b-i);
            fmpz_divexact2_uiui(coeffs + i, coeffs + i, i, i-1);
            fmpz_neg(coeffs + i, coeffs + i);
        }

        return;
    }

    prec = magnitude_bound(d) + 5 + FLINT_BIT_COUNT(d);

    alpha = _fmpz_vec_init(d);
    fmpz_init(half);
    mpfr_init2(t, prec);
    mpfr_init2(u, prec);

    fmpz_one(half);
    fmpz_mul_2exp(half, half, prec - 1);
    mpfr_const_pi(t, prec);
    mpfr_div_ui(t, t, n, MPFR_RNDN);

    for (i = j = 0; j < d; i++)
    {
        if (n_gcd(n, i) == 1)
        {
            mpfr_mul_ui(u, t, 2 * i, MPFR_RNDN);
            mpfr_cos(u, u, MPFR_RNDN);
            mpfr_neg(u, u, MPFR_RNDN);
            exp = mpfr_get_z_2exp(_fmpz_promote(alpha + j), u);
            _fmpz_demote_val(alpha + j);
            fmpz_mul_or_div_2exp(alpha + j, alpha + j, exp + prec);
            j++;
        }
    }

    balanced_product(coeffs, alpha, d, prec);

    /* Scale and round */
    for (i = 0; i < d + 1; i++)
    {
        len_t r = d;
        if ((n & (n - 1)) == 0)
            r--;
        fmpz_mul_2exp(coeffs + i, coeffs + i, r);
        fmpz_add(coeffs + i, coeffs + i, half);
        fmpz_fdiv_q_2exp(coeffs + i, coeffs + i, prec);
    }

    fmpz_clear(half);
    mpfr_clear(t);
    mpfr_clear(u);
    _fmpz_vec_clear(alpha, d);
}

void
arith_cos_minpoly(fmpz_poly_t poly, ulong n)
{
    if (n == 0)
    {
        fmpz_poly_set_ui(poly, 1UL);
    }
    else
    {
        len_t d = (n <= 2) ? 1 : n_euler_phi(n) / 2;

        fmpz_poly_fit_length(poly, d + 1);
        _arith_cos_minpoly(poly->coeffs, d, n);
        _fmpz_poly_set_length(poly, d + 1);
    }
}
