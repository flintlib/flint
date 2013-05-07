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
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_vec.h"
#include "arith.h"
#include "ulong_extras.h"

void _arith_stirling_next_row(fmpz * new, fmpz * prev, long n, 
                                               long klen, int kind);

static __inline__ void
_fmpz_addmul_alt(fmpz_t s, fmpz_t t, fmpz_t u, int parity)
{
    if (parity % 2)
        fmpz_submul(s, t, u);
    else
        fmpz_addmul(s, t, u);
}

static void
_fmpz_stirling2_powsum(fmpz_t s, long n, long k)
{
    fmpz_t t, u;
    fmpz * bc;
    long j, m, max_bc;

    fmpz_init(t);
    fmpz_init(u);
    max_bc = (k+1) / 2;

    bc = _fmpz_vec_init(max_bc + 1);
    fmpz_one(bc);
    for (j = 1; j <= max_bc; j++)
    {
        fmpz_set(bc+j, bc+j-1);
        fmpz_mul_ui(bc+j, bc+j, k+1-j);
        fmpz_divexact_ui(bc+j, bc+j, j);
    }

    fmpz_zero(s);
    for (j = 1; j <= k; j += 2)
    {
        fmpz_set_ui(u, j);
        fmpz_pow_ui(u, u, n);
        m = j;
        /* Process each m = 2^p * j */
        while (1)
        {
            if (m > max_bc)
                _fmpz_addmul_alt(s, bc+k-m, u, k + m);
            else
                _fmpz_addmul_alt(s, bc+m, u, k + m);
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

void
arith_stirling_number_2(fmpz_t s, long n, long k)
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
            fmpz_one(s);
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
            fmpz_one(s);
            fmpz_mul_2exp(s, s, n-1);
            fmpz_sub_ui(s, s, 1UL);
        }
        return;
    }

    _fmpz_stirling2_powsum(s, n, k);
}

void
arith_stirling_number_2_vec(fmpz * row, long n, long klen)
{
    long m;

    klen = FLINT_MIN(klen, n+1);

    if (klen < 1)
        return;

    fmpz_one(row);
    for (m = 1; m <= n; m++)
        _arith_stirling_next_row(row, row, m, klen, 2);
}
