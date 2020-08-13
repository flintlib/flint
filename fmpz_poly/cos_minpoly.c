/*
    Copyright (C) 2014 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "fmpz_poly.h"
#include "ulong_extras.h"

/* Use a lookup table for small n. We skip 53, 59 and 61, as the
   coefficients do not fit in 16 bits. */
#define MINPOLY_TAB_NUM 65
#define MINPOLY_TAB_MAX_LEN 24

static const char
minpoly_len_tab[MINPOLY_TAB_NUM] = {
    1, 2, 2, 2, 2, 3, 2, 4, 3, 4, 3, 6, 3, 7, 4, 5, 5, 9, 4, 10, 5, 7, 6,
    12, 5, 11, 7, 10, 7, 15, 5, 16, 9, 11, 9, 13, 7, 19, 10, 13, 9, 21, 7,
    22, 11, 13, 12, 24, 9, 22, 11, 17, 13, 27, 10, 21, 13, 19, 15, 30, 9,
    31, 16, 19, 17
};

static const short
minpoly_tab[MINPOLY_TAB_NUM][MINPOLY_TAB_MAX_LEN] = {
    {1},
    {-2, 1},
    {2, 1},
    {1, 1},
    {0, 1},
    {-1, 1, 1},
    {-1, 1},
    {-1, -2, 1, 1},
    {-2, 0, 1},
    {1, -3, 0, 1},
    {-1, -1, 1},
    {1, 3, -3, -4, 1, 1},
    {-3, 0, 1},
    {-1, 3, 6, -4, -5, 1, 1},
    {1, -2, -1, 1},
    {1, 4, -4, -1, 1},
    {2, 0, -4, 0, 1},
    {1, -4, -10, 10, 15, -6, -7, 1, 1},
    {-1, -3, 0, 1},
    {1, 5, -10, -20, 15, 21, -7, -8, 1, 1},
    {5, 0, -5, 0, 1},
    {1, -8, 8, 6, -6, -1, 1},
    {-1, 3, 3, -4, -1, 1},
    {-1, -6, 15, 35, -35, -56, 28, 36, -9, -10, 1, 1},
    {1, 0, -4, 0, 1},
    {-1, 5, 25, -5, -50, 1, 35, 0, -10, 0, 1},
    {-1, -3, 6, 4, -5, -1, 1},
    {1, 9, 0, -30, 0, 27, 0, -9, 0, 1},
    {-7, 0, 14, 0, -7, 0, 1},
    {-1, 7, 28, -56, -126, 126, 210, -120, -165, 55, 66, -12, -13, 1, 1},
    {1, -4, -4, 1, 1},
    {-1, -8, 28, 84, -126, -252, 210, 330, -165, -220, 66, 78, -13, -14, 1, 1},
    {2, 0, -16, 0, 20, 0, -8, 0, 1},
    {1, -12, 12, 43, -43, -34, 34, 10, -10, -1, 1},
    {1, 4, -10, -10, 15, 6, -7, -1, 1},
    {1, 8, -40, -46, 110, 71, -113, -43, 54, 11, -12, -1, 1},
    {-3, 0, 9, 0, -6, 0, 1},
    {-1, 9, 45, -120, -330, 462, 924, -792, -1287, 715, 1001, -364, -455, 105,
        120, -16, -17, 1, 1},
    {-1, 5, 10, -20, -15, 21, 7, -8, -1, 1},
    {1, 12, -12, -79, 79, 103, -103, -53, 53, 12, -12, -1, 1},
    {1, 0, -12, 0, 19, 0, -8, 0, 1},
    {1, -10, -55, 165, 495, -792, -1716, 1716, 3003, -2002, -3003, 1365, 1820,
        -560, -680, 136, 153, -18, -19, 1, 1},
    {1, 8, 8, -6, -6, 1, 1},
    {1, 11, -55, -220, 495, 1287, -1716, -3432, 3003, 5005, -3003, -4368,
        1820, 2380, -680, -816, 153, 171, -19, -20, 1, 1},
    {-11, 0, 55, 0, -77, 0, 44, 0, -11, 0, 1},
    {1, -12, -36, 31, 105, -27, -112, 9, 54, -1, -12, 0, 1},
    {1, -6, -15, 35, 35, -56, -28, 36, 9, -10, -1, 1},
    {-1, -12, 66, 286, -715, -2002, 3003, 6435, -6435, -11440, 8008, 12376,
        -6188, -8568, 3060, 3876, -969, -1140, 190, 210, -21, -22, 1, 1},
    {1, 0, -16, 0, 20, 0, -8, 0, 1},
    {-1, 14, 49, -371, -196, 2072, 294, -5147, -210, 7007, 77, -5733, -14,
        2940, 1, -952, 0, 189, 0, -21, 0, 1},
    {-1, -5, 25, 5, -50, -1, 35, 0, -10, 0, 1},
    {1, 16, -16, -188, 188, 526, -526, -596, 596, 339, -339, -103, 103, 16,
        -16, -1, 1},
    {13, 0, -91, 0, 182, 0, -156, 0, 65, 0, -13, 0, 1},
    {0},
    {-1, 9, 0, -30, 0, 27, 0, -9, 0, 1},
    {1, 12, -108, -151, 951, 877, -2891, -2058, 4489, 2442, -4080, -1639,
        2289, 650, -801, -151, 170, 19, -20, -1, 1},
    {1, 0, -24, 0, 86, 0, -104, 0, 53, 0, -12, 0, 1},
    {1, -20, 20, 265, -265, -989, 989, 1519, -1519, -1198, 1198, 531, -531,
        -134, 134, 18, -18, -1, 1},
    {-1, -7, 28, 56, -126, -126, 210, 120, -165, -55, 66, 12, -13, -1, 1},
    {0},
    {1, 0, -8, 0, 14, 0, -7, 0, 1},
    {0},
    {1, -8, -28, 84, 126, -252, -210, 330, 165, -220, -66, 78, 13, -14, -1, 1},
    {1, 24, 72, -170, -534, 405, 1385, -459, -1782, 276, 1287, -90, -546, 15,
        135, -1, -18, 0, 1},
    {2, 0, -64, 0, 336, 0, -672, 0, 660, 0, -352, 0, 104, 0, -16, 0, 1},
};

/* Recurrence for coefficients in rescaled Chebyshev polynomials */
#define CHEB_NEXT(y, x, m, k) \
    fmpz_mul2_uiui(y, x, m - 2*k + 1, m - 2*k + 2); \
    fmpz_divexact2_uiui(y, y, k, m - k); \
    fmpz_neg(y, y);  \

/* Computes the monic integer polynomial
    n odd: 2 (T(s+1,x/2) - T(s,x/2)), s = (n - 1) / 2
    n even: 2 (T(s+1,x/2) - T(s-1,x/2)), s = n / 2 */
static void
chebyshev_sum(fmpz * a, ulong n)
{
    ulong s, k, m;

    if (n == 1)
    {
        fmpz_set_si(a, -2);
        fmpz_one(a + 1);
        return;
    }

    if (n == 2)
    {
        fmpz_set_si(a, -4);
        fmpz_zero(a + 1);
        fmpz_one(a + 2);
        return;
    }

    s = n / 2;
    m = s + 1;

    fmpz_one(a + m);
    for (k = 1; k <= m / 2; k++)
    {
        CHEB_NEXT(a + m - 2 * k, a + m - 2 * k + 2, m, k);
    }

    if (n % 2 == 1)
    {
        m = s;

        fmpz_set_si(a + m, -1);
        for (k = 1; k <= m / 2; k++)
        {
            CHEB_NEXT(a + m - 2 * k, a + m - 2 * k + 2, m, k);
        }
    }
    else
    {
        m = s - 1;

        /* Use the top coefficient as scratch space. */
        for (k = 1; k <= m / 2; k++)
        {
            CHEB_NEXT(a + m + 2, a + m + 2, m, k);
            fmpz_sub(a + m - 2*k, a + m - 2*k, a + m + 2);
        }

        for (k = 1 - (m % 2); k < m + 2; k += 2)
            fmpz_zero(a + k);

        fmpz_sub_ui(a + m, a + m, 1);
        /* Set the top coefficient again. */
        fmpz_one(a + m + 2);
    }
}

#define MUL_TMP(P, Plen, T, Tlen) \
    fmpz * swap; \
    if (Plen >= Tlen) \
        _fmpz_poly_mul(U, P, Plen, T, Tlen); \
    else \
        _fmpz_poly_mul(U, T, Tlen, P, Plen); \
    Plen = Plen + Tlen - 1; \
    swap = P; P = U; U = swap; \

void
_fmpz_poly_cos_minpoly(fmpz * f, ulong n)
{
    fmpz *P, *Q, *T, *U;
    int *mu;
    ulong Pdeg, Qdeg;
    ulong Plen, Qlen, Tlen;
    ulong d;

    if (n < MINPOLY_TAB_NUM && minpoly_len_tab[n] <= MINPOLY_TAB_MAX_LEN)
    {
        for (d = 0; d < minpoly_len_tab[n]; d++)
            fmpz_set_si(f + d, minpoly_tab[n][d]);
        return;
    }

    /* Compute values of the Moebius function. We do this as a precomputation
       as it allows us to bound in advance the degrees of the numerator and
       denominator. */
    mu = flint_calloc(n + 1, sizeof(int));
    Pdeg = Qdeg = 0;

    for (d = 1; d <= n; d++)
    {
        if (n % d == 0)
        {
            mu[d] = n_moebius_mu(n / d);
            if (mu[d] == 1)
                Pdeg += (d / 2 + 1);
            else if (mu[d] == -1)
                Qdeg += (d / 2 + 1);
        }
    }

    /* We use two extra arrays as scratch space (note that Qdeg < Pdeg). */
    P = _fmpz_vec_init(Pdeg + 1);
    Q = _fmpz_vec_init(Pdeg + 1);
    T = _fmpz_vec_init(Pdeg + 1);
    U = _fmpz_vec_init(Pdeg + 1);

    Plen = Qlen = 1;
    fmpz_one(P);
    fmpz_one(Q);

    for (d = 1; d <= n; d++)
    {
        if (n % d == 0 && mu[d] != 0)
        {
            chebyshev_sum(T, d);
            Tlen = d / 2 + 2;

            if (mu[d] > 0)
            {
                MUL_TMP(P, Plen, T, Tlen);
            }
            else
            {
                MUL_TMP(Q, Qlen, T, Tlen);
            }
        }
    }

    _fmpz_poly_div(f, P, Plen, Q, Qlen, 0);

    _fmpz_vec_clear(P, Pdeg + 1);
    _fmpz_vec_clear(Q, Pdeg + 1);
    _fmpz_vec_clear(T, Pdeg + 1);
    _fmpz_vec_clear(U, Pdeg + 1);
    flint_free(mu);
}

void
fmpz_poly_cos_minpoly(fmpz_poly_t f, ulong n)
{
    slong len = (n < MINPOLY_TAB_NUM) ?
        minpoly_len_tab[n] : n_euler_phi(n) / 2 + 1;

    fmpz_poly_fit_length(f, len);
    _fmpz_poly_cos_minpoly(f->coeffs, n);
    _fmpz_poly_set_length(f, len);
}
