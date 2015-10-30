#include <math.h>

#include "flint.h"
#include "fmpz.h"
#include "fmpz_poly.h"
#include "fmpz_mat.h"
#include "nmod_mat.h"
#include "nmod_poly.h"

#define M_LOG2E  1.44269504088896340736  /* log2(e) */

static __inline__ long double _log2(const long double x)
{
    return log(x) * M_LOG2E;
}

slong _fmpz_mat_minpoly_small(fmpz * rop, const fmpz_mat_t op)
{
    slong len = 0;

    if (op->r == 0)
    {
        fmpz_one(rop + 0);
        len = 1;
    }
    else if (op->r == 1)
    {
        fmpz_one(rop + 1);
        fmpz_neg(rop + 0, &(op->rows[0][0]));
        len = 2;
    }

    return len;
}

slong _fmpz_mat_minpoly_modular(fmpz * rop, const fmpz_mat_t op)
{
    const slong n = op->r;
    slong len = 0;

    if (n < 2)
    {
        return _fmpz_mat_minpoly_small(rop, op);
    }
    else
    {
        /*
            If $A$ is an $n \times n$ matrix with $n \geq 4$ and 
            coefficients bounded in absolute value by $B > 1$ then 
            the coefficients of the characteristic polynomial have 
            less than $\ceil{n/2 (\log_2(n) + \log_2(B^2) + 1.6669)}$ 
            bits.
            See Lemma 4.1 in Dumas, Pernet, and Wan, "Efficient computation 
            of the characteristic polynomial", 2008.
            
            The coefficients of the minimal polynomial have at most
            n more bits. See Dumas, "Bounds on the coefficients of the
            characteristic and minimal polynomials, 2007, especially
            section 3.

            From the explicit formulae for charpolys, the absolute values
            of the coefficients in the case n = 2 and n = 3 are at most
            2B^2 and 6B^3 respectively.

            Thus we have bounds on the coefficients of the minimal
            polynomial of $2\log_2(B) + 3$ and 3\log(B) + 6$.
         */
        slong bound;

        slong pbits  = FLINT_BITS - 1;
        mp_limb_t p = (1UL << pbits);

        fmpz_t m;

        /* Determine the bound in bits */
        {
            slong i, j;
            fmpz *ptr;
            double t;

            ptr = fmpz_mat_entry(op, 0, 0);
            for (i = 0; i < n; i++)
                for (j = 0; j < n; j++)
                    if (fmpz_cmpabs(ptr, fmpz_mat_entry(op, i, j)) < 0)
                        ptr = fmpz_mat_entry(op, i, j);

            if (fmpz_bits(ptr) == 0)  /* Zero matrix */
            {
                for (i = 0; i < n; i++)
                   fmpz_zero(rop + i);
                fmpz_set_ui(rop + n, 1);
                return n + 1;
            }

            t = (fmpz_bits(ptr) <= FLINT_D_BITS) ? 
                _log2(FLINT_ABS(fmpz_get_d(ptr))) : fmpz_bits(ptr);

            if (n == 2)
               bound = ceil(2.0 * t + 3.0);
            else if (n == 3)
               bound = ceil(3.0 * t + 6.0);
            else
               bound = ceil( (n / 2.0) * (_log2(n) + 2.0 * t + 1.6669) ) + n;
        }

        fmpz_init_set_ui(m, 1);

        len = 0;

        for ( ; fmpz_bits(m) < bound; )
        {
            nmod_mat_t mat;
            nmod_poly_t poly;

            p = n_nextprime(p, 0);

            nmod_mat_init(mat, n, n, p);
            nmod_poly_init(poly, p);

            fmpz_mat_get_nmod_mat(mat, op);
            nmod_mat_minpoly(poly, mat);

            len = FLINT_MAX(len, poly->length);

            _fmpz_poly_CRT_ui(rop, rop, n + 1, m, poly->coeffs, poly->length, poly->mod.n, poly->mod.ninv, 1);

            fmpz_mul_ui(m, m, p);

            nmod_mat_clear(mat);
            nmod_poly_clear(poly);
        }

        fmpz_clear(m);
    }

    return len;
}

void fmpz_mat_minpoly_modular(fmpz_poly_t cp, const fmpz_mat_t mat)
{
    slong len;

    fmpz_poly_fit_length(cp, mat->r + 1);

    len = _fmpz_mat_minpoly_modular(cp->coeffs, mat);

    _fmpz_poly_set_length(cp, len);
}
