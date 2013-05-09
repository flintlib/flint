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

    Copyright (C) 2010 Sebastian Pancratz

******************************************************************************/

#include <stdlib.h>
#include <gmp.h>
#include "flint.h"
#include "fmpz.h"
#include "fmpz_vec.h"
#include "fmpz_poly.h"

/*
    Assumptions.

    Suppose that $len1 \geq 3$ and $len2 \geq 2$.


    Definitions.

    Define a sequence $(n_i)$ by $n_1 = \ceil{len1 / 2}$, 
    $n_2 = \ceil{n_1 / 2}$, etc. all the way to 
    $n_K = \ceil{n_{K-1} / 2} = 2$.  Thus, $K = \ceil{\log_2 len1} - 1$. 
    Note that we can write $n_i = \ceil{len1 / 2^i}$.


    Rough description (of the allocation process, or the algorithm).

    Step 1.
    For $0 \leq i < n_1$, set h[i] to something of length at most len2.
    Set pow to $poly2^2$.
    
    Step n.
    For $0 \leq i < n_n$, set h[i] to something of length at most the length 
    of $poly2^{2^n - 1}$.
    Set pow to $poly^{2^n}$.
    
    Step K.
    For $0 \leq i < n_K$, set h[i] to something of length at most the length 
    of $poly2^{2^K - 1}$.
    Set pow to $poly^{2^K}$.


    Analysis of the space requirements.

    Let $S$ be the over all space we need, measured in number of coefficients.
    Then 
    \begin{align*}
    S & = 2 \times \bigl[ (2^K - 1) (len2 - 1) + 1 \bigr] 
        + \sum_{i=1}^{K-1} (n_i - n_{i+1}) \bigl[(2^i - 1) (len2 - 1) + 1\bigr] \\
      & = 2 \times \bigl[ (2^K - 1) (len2 - 1) + 1 \bigr] 
        + (len2 - 1) \sum_{i=1}^{K-1} (n_i - n_{i+1}) (2^i - 1) + n_1 - n_K.
    \end{align*}

    If $K = 1$, or equivalently $len1$ is 3 or 4, then $S = 2 \times len2$.
    Otherwise, we can bound $n_i - n_{i+1}$ from above as follows.  For any 
    non-negative integer $x$, 
    \begin{equation*}
    \ceil{x / 2^i} - \ceil{x / 2^{i+1}} \leq x/2^i - x/2^{i+1} = x / 2^{i+1}.
    \end{equation*}

    Thus, 
    \begin{align*}
    S & \leq 2 \times \bigl[ (2^K - 1) (len2 - 1) + 1 \bigr] 
           + (len2 - 1) \times len1 \times \sum_{i=1}^{K-1} (1/2 - 1/2^{i+1}) \\
      & \leq 2 \times \bigl[ (2^K - 1) (len2 - 1) + 1 \bigr] 
           + (len2 - 1) \times len1 \times (K/2 + 1).
    \end{align*}
 */

void
_fmpz_poly_compose_divconquer(fmpz * res, const fmpz * poly1, len_t len1, 
                                          const fmpz * poly2, len_t len2)
{
    len_t i, j, k, n;
    len_t *hlen, alloc, powlen;
    fmpz *v, **h, *pow, *temp;

    if (len1 <= 2 || len2 <= 1)
    {
        if (len1 == 1)
            fmpz_set(res, poly1);
        else if (len2 == 1)
            _fmpz_poly_evaluate_fmpz(res, poly1, len1, poly2);
        else  /* len1 == 2 */
            _fmpz_poly_compose_horner(res, poly1, len1, poly2, len2);
        return;
    }

    /* Initialisation */
    
    hlen = (len_t *) flint_malloc(((len1 + 1) / 2) * sizeof(len_t));
    
    k = FLINT_CLOG2(len1) - 1;
    
    hlen[0] = hlen[1] = ((1 << k) - 1) * (len2 - 1) + 1;
    for (i = k - 1; i > 0; i--)
    {
        len_t hi = (len1 + (1 << i) - 1) / (1 << i);
        for (n = (hi + 1) / 2; n < hi; n++)
            hlen[n] = ((1 << i) - 1) * (len2 - 1) + 1;
    }
    powlen = (1 << k) * (len2 - 1) + 1;
    
    alloc = 0;
    for (i = 0; i < (len1 + 1) / 2; i++)
        alloc += hlen[i];

    v = _fmpz_vec_init(alloc +  2 * powlen);
    h = (fmpz **) flint_malloc(((len1 + 1) / 2) * sizeof(fmpz *));
    h[0] = v;
    for (i = 0; i < (len1 - 1) / 2; i++)
    {
        h[i + 1] = h[i] + hlen[i];
        hlen[i]  = 0;
    }
    hlen[(len1 - 1) / 2] = 0;
    pow  = v + alloc;
    temp = pow + powlen;
    
    /* Let's start the actual work */
    
    for (i = 0, j = 0; i < len1 / 2; i++, j += 2)
    {
        if (poly1[j + 1] != 0L)
        {
            _fmpz_vec_scalar_mul_fmpz(h[i], poly2, len2, poly1 + j + 1);
            fmpz_add(h[i], h[i], poly1 + j);
            hlen[i] = len2;
        }
        else if (poly1[j] != 0L)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }
    if ((len1 & 1L))
    {
        if (poly1[j] != 0L)
        {
            fmpz_set(h[i], poly1 + j);
            hlen[i] = 1;
        }
    }
    
    _fmpz_poly_sqr(pow, poly2, len2);
    powlen = 2 * len2 - 1;
    
    for (n = (len1 + 1) / 2; n > 2; n = (n + 1) / 2)
    {
        if (hlen[1] > 0)
        {
            len_t templen = powlen + hlen[1] - 1;
            _fmpz_poly_mul(temp, pow, powlen, h[1], hlen[1]);
            _fmpz_poly_add(h[0], temp, templen, h[0], hlen[0]);
            hlen[0] = FLINT_MAX(hlen[0], templen);
        }
        
        for (i = 1; i < n / 2; i++)
        {
            if (hlen[2*i + 1] > 0)
            {
                _fmpz_poly_mul(h[i], pow, powlen, h[2*i + 1], hlen[2*i + 1]);
                hlen[i] = hlen[2*i + 1] + powlen - 1;
            } else
                hlen[i] = 0;
            _fmpz_poly_add(h[i], h[i], hlen[i], h[2*i], hlen[2*i]);
            hlen[i] = FLINT_MAX(hlen[i], hlen[2*i]);
        }
        if ((n & 1L))
        {
            _fmpz_vec_set(h[i], h[2*i], hlen[2*i]);
            hlen[i] = hlen[2*i];
        }
        
        _fmpz_poly_sqr(temp, pow, powlen);
        powlen += powlen - 1;
        {
            fmpz * t = pow;
            pow      = temp;
            temp     = t;
        }
    }

    _fmpz_poly_mul(res, pow, powlen, h[1], hlen[1]);
    _fmpz_vec_add(res, res, h[0], hlen[0]);
    
    _fmpz_vec_clear(v, alloc + 2 * powlen);
    flint_free(h);
    flint_free(hlen);
}

void
fmpz_poly_compose_divconquer(fmpz_poly_t res, 
                             const fmpz_poly_t poly1, const fmpz_poly_t poly2)
{
    const len_t len1 = poly1->length;
    const len_t len2 = poly2->length;
    len_t lenr;
    
    if (len1 == 0)
    {
        fmpz_poly_zero(res);
        return;
    }
    if (len1 == 1 || len2 == 0)
    {
        fmpz_poly_set_fmpz(res, poly1->coeffs);
        return;
    }
    
    lenr = (len1 - 1) * (len2 - 1) + 1;
    
    if (res != poly1 && res != poly2)
    {
        fmpz_poly_fit_length(res, lenr);
        _fmpz_poly_compose_divconquer(res->coeffs, poly1->coeffs, len1, 
                                                   poly2->coeffs, len2);
        _fmpz_poly_set_length(res, lenr);
        _fmpz_poly_normalise(res);
    }
    else
    {
        fmpz_poly_t t;
        fmpz_poly_init2(t, lenr);
        _fmpz_poly_compose_divconquer(t->coeffs, poly1->coeffs, len1,
                                                 poly2->coeffs, len2);
        _fmpz_poly_set_length(t, lenr);
        _fmpz_poly_normalise(t);
        fmpz_poly_swap(res, t);
        fmpz_poly_clear(t);
    }
}
