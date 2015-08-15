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

    Copyright (C) 2015 Kushagra Singh

******************************************************************************/

void
mpqs_compute_A(mpqs_t mpqs_inf)
{
    mpqs_inf->current_prime = n_nextprime(mpqs_inf->A_targetprime, 0);

    fmpz_set_ui(mpqs_inf->A, mpqs_inf->A_targetprime);
    fmpz_mul(mpqs_inf->A, mpqs_inf->A, mpqs_inf->A);
}

void
mpqs_compute_B(mpqs_t mpqs_inf)
{
    mp_limb_t prime = mpqs_inf->A_targetprime;

    /* odd root of x^2 = n mod A */
    /* A is A_targetprime ^ 2 */

    mpqs_sqrtmod_psq(mpqs_inf->B, mpqs_inf->kn, prime);
}

void
mpqs_compute_C(mpqs_t mpqs_inf)
{
    /* C = (B^2 - kn)/4A */

    fmpz_set(mpqs_inf->C, mpqs_inf->B);
    fmpz_mul(mpqs_inf->C, mpqs_inf->C, mpqs_inf->C);      /* B^2 */
    fmpz_sub(mpqs_inf->C, mpqs_inf->C, mpqs_inf->kn);     /* B^2 - kn */
    fmpz_fdiv_q(mpqs_inf->C, mpqs_inf->C, mpqs_inf->A);   /* (B^2 - kn)/A */
}
