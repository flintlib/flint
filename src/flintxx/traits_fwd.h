/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

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
