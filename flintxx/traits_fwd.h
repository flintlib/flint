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

    Copyright (C) 2013 Tom Bachmann

******************************************************************************/

#ifndef FLINTXX_TRAITS_FWD_H
#define FLINTXX_TRAITS_FWD_H

// Use this to break cyclic dependencies, by forward-declaring traits.

namespace flint {
namespace traits {
template<class T> struct is_fmpz_matxx;
template<class T> struct is_fmpq_matxx;
template<class T> struct is_fmpz_polyxx;
template<class T> struct is_nmod_matxx;
template<class T> struct is_nmod_polyxx;
} // traits
} // flint

#endif
