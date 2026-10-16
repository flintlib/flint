/*
    Copyright (C) 2026 Fredrik Johansson

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#include "acb.h"

#define SECOND_ZETA_BOUND_PREC MAG_BITS

/*
    Internal routines for the secondary zeta function

        Z(s) = sum_{n>=1} alpha_n^{-s},   rho_n = 1/2 + i alpha_n,

    following J. Arias de Reyna, "Computation of the secondary zeta
    function", arXiv:2006.04869.  With a free parameter a > 0,

        Z(s) = A(s) - P(s) + E(s) - S(s).

    Each of the four helpers below writes a rigorous enclosure of its term
    (midpoint + error radius) to res.  The error radii implement the bounds:

        A : Trudgian-majorant tail bound (companion note, Thm 5.1);
        P : tightened erfc / psi tail bound (companion note, eqs 4.1/4.2);
        S : effective Watson-remainder bound (companion note, Thm/eq 4.3);
        E : trivial tail of an absolutely convergent power series (eq 47).

    All four are unconditional.  The zeros used by A() are obtained from
    acb_dirichlet_hardy_z_zeros (rigorously isolated); for heights below the
    Platt-Trudgian verification bound 3e12 they are simple and on the line.
*/

#define SECOND_ZETA_RH_VERIFIED_HEIGHT 3e12

/* term E(s): cheap, computed first to size the working precision. */
void _acb_dirichlet_secondary_zeta_term_E(acb_t res, const acb_t s,
        const arb_t a, slong prec);

/* term S(s): singular term, asymptotic bracket of eq (53) truncated at N. */
void _acb_dirichlet_secondary_zeta_term_S(acb_t res, const acb_t s,
        const arb_t a, slong prec);

/* term P(s): prime term, eq (38) summed to the rigorous stopping index. */
void _acb_dirichlet_secondary_zeta_term_P(acb_t res, const acb_t s,
        const arb_t a, slong prec);

/* term A(s): zero term, eq (28) summed over isolated zeta zeros. */
void _acb_dirichlet_secondary_zeta_term_A(acb_t res, const acb_t s,
        const arb_t a, slong prec);

/*
    Thread-local cached accessor for the first 'num' imaginary ordinates
    gamma_1,...,gamma_num of the zeta zeros (rho = 1/2 + i gamma), each rounded
    to 'prec' bits.  Extends the cache in both length and precision on demand.
    res must have space for 'num' arb entries.
*/
void _acb_dirichlet_secondary_zeta_zeros(arb_ptr res, slong num, slong prec);

/* public entry: Z(s) = sum_n alpha_n^{-s} over zeros rho_n = 1/2 + i alpha_n */
void acb_dirichlet_secondary_zeta(acb_t res, const acb_t s, slong prec);

