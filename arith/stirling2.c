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
#include <mpir.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "ulong_extras.h"


void fmpz_addmul_alt(fmpz_t s, fmpz_t t, fmpz_t u, int parity)
{
    if (parity % 2)
        fmpz_submul(s, t, u);
    else
        fmpz_addmul(s, t, u);
}


void _fmpz_stirling2_powsum(fmpz_t s, long n, long k)
{
    fmpz_t t, u;
    long j;

    fmpz_init(t);
    fmpz_init(u);
    fmpz_set_ui(t, 1UL);
    fmpz_set_ui(s, 0UL);

    for (j = 1; j < k/2 + 1; j++)
    {
        fmpz_mul_ui(t, t, k+1-j);
        fmpz_divexact_ui(t, t, j);
        fmpz_set_ui(u, j);
        fmpz_pow_ui(u, u, n);
        fmpz_addmul_alt(s, t, u, k+j);
        if (2*j != k)
        {
            /* C(k,j) = C(k,k-j) */
            fmpz_set_ui(u, k-j);
            fmpz_pow_ui(u, u, n);
            fmpz_addmul_alt(s, t, u, j);
        }
    }

    /* Last term not included because loop starts from 1 */
    fmpz_set_ui(u, k);
    fmpz_pow_ui(u, u, n);
    fmpz_add(s, s, u);
    fmpz_fac_ui(t, k);
    fmpz_divexact(s, s, t);
    fmpz_clear(t);
    fmpz_clear(u);
}

void _fmpz_stirling2_powsum_odd(fmpz_t s, long n, long k)
{
    fmpz_t t, u;
    fmpz * bc;
    long j, m, max_bc;

    fmpz_init(t);
    fmpz_init(u);
    max_bc = (k+1) / 2;

    bc = _fmpz_vec_init(max_bc + 1);
    fmpz_set_ui(bc, 1UL);
    for (j = 1; j < max_bc + 1; j++)
    {
        fmpz_set(bc+j, bc+j-1);
        fmpz_mul_ui(bc+j, bc+j, k+1-j);
        fmpz_divexact_ui(bc+j, bc+j, j);
    }

    fmpz_zero(s);
    for (j = 1; j < k + 1; j += 2)
    {
        fmpz_set_ui(u, j);
        fmpz_pow_ui(u, u, n);
        m = j;
        /* Process each m = 2^p * j */
        while (1)
        {
            if (m > max_bc)
                fmpz_addmul_alt(s, bc+k-m, u, k + m);
            else
                fmpz_addmul_alt(s, bc+m, u, k + m);
            m *= 2;
            if (m > k)
                break;
            fmpz_mul_2exp(u, u, n);
        }
    }

    _fmpz_vec_clear(bc, max_bc + 1);
    fmpz_fac_ui(t, k);
    fmpz_divexact(s, s, t);
    fmpz_clear(t);
    fmpz_clear(u);
}



void fmpz_stirling2(fmpz_t s, long n, long k)
{
    if (n < 0 || k < 0 || k > n)
    {
        fmpz_zero(s);
        return;
    }

    /* Topmost diagonals */
    if (k >= n - 1)
    {
        if (k == n)
            fmpz_set_ui(s, 1UL);
        else /* k == n - 1 */
        {
            /* S(n,n-1) = binomial(n,2) */
            fmpz_set_ui(s, n);
            fmpz_mul_ui(s, s, n-1);
            fmpz_divexact_ui(s, s, 2UL);
        }
        return;
    }

    /* Leftmost columns */
    if (k <= 2)
    {
        if (k < 2)
            fmpz_set_ui(s, k);
        else
        {
            /* S(n,2) = 2^(n-1)-1 */
            fmpz_set_ui(s, 1UL);
            fmpz_mul_2exp(s, s, n-1);
            fmpz_sub_ui(s, s, 1UL);
        }
        return;
    }

    if (n < 200)
        _fmpz_stirling2_powsum(s, n, k);
    else
        _fmpz_stirling2_powsum_odd(s, n, k);
}


void
fmpz_stirling2_vec(fmpz * row, long n, long klen)
{
    long k;
    for (k = 0; k < klen; k++)
        fmpz_stirling2(row + k, n, k);
}
