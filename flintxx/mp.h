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

#ifndef CXX_MP_H
#define CXX_MP_H

namespace flint {
namespace mp {
/////////////////////////////////
// BASIC METAPROGRAMMING HELPERS
/////////////////////////////////
// Most of these helpers *compute* something. In this case, they have a static
// static const member "val" storing the result of the computation. The
// arguments of the computation are template parameters, and are usually *not*
// PODs, but instead themselves types with a static constant "val" member.
// See value_of, true_ and false_ for how to pass in explicit values to your
// computation.

// Wrap the boolean value "v" into the static const member "val".
template<bool v> struct value_of { static const bool val = v; };

// Boolean inputs for computations.
struct true_ : value_of<true> { };
struct false_ : value_of<false> { };

// Compute if two input *types* (not values!) are equal.
template<class T, class U>
struct equal_types : false_ { };
template<class T>
struct equal_types<T, T> : true_ { };

// Compute logical negation of the input value.
template<class T>
struct not_ : value_of<!T::val> { };

// Compute logical and of the input values.
template<class T1, class T2, class T3 = void,
    class T4 = void, class T5 = void, class T6 = void, class T7 = void, class T8 = void>
struct and_ : and_<T1, and_<T2, T3, T4, T5, T6, T7, T8> > { };

template<class T, class U>
struct and_<T, U, void, void, void, void> : value_of<T::val && U::val> { };

template<class T, bool u>
struct and_v : and_<T, value_of<u> > { };

// Compute logical or of the input values.
template<class T1, class T2, class T3 = void,
    class T4 = void, class T5 = void, class T6 = void>
struct or_ : or_<T1, or_<T2, T3, T4, T5> > { };

template<class T, class U>
struct or_<T, U, void, void, void, void> : value_of<T::val || U::val> { };

// Compute V1 or V2, depending on C
template<bool C, class V1, class V2>
struct if_v {typedef V1 type;};
template<class V1, class V2>
struct if_v<false, V1, V2> {typedef V2 type;};

template<class C, class V1, class V2>
struct if_ : if_v<C::val, V1, V2> { };

// Choose a value depending on a sequence of conditions.
// This has he same meaning as
//  int select(bool c1 = false, bool c2 = false, bool c3 = false)
//  {
//      if(c1)
//          return v1;
//      if(c2)
//          return v2;
//      if(c3)
//          return v3;
//      return d;
//  }
template<class D, class C1 = void, class V1 = void,
    class C2 = void, class V2 = void,
    class C3 = void, class V3 = void>
struct select : if_<C1, V1, typename select<D, C2, V2, C3, V3>::type> { };

template<class D>
struct select<D,  void, void,  void, void,  void, void> {typedef D type;};

// Conditional template enabling helper. (See below for explanation.)
template<bool, class U = void>
struct enable_if_v
{
    typedef U type;
    static const int val = 0;
};
template<class U>
struct enable_if_v<false, U> { };

// Conditional template enabling.
//
// These two helpers (enable_if and disable_if) can be used wherever the
// SFINAE rule applies to conditionally enable or disable template
// specialisations.
//
// If T evaluaties to true, then enable_if has a member typedef "type", which
// is the parameter U, and also a static const int member val (of zero).
// If on the other hand T evaluates to false, then enable_if is empty.
// The meaning of T is reversed for disable_if.
// See e.g. [0] or the tests for how to use this.
// 
// [0] http://www.boost.org/doc/libs/1_53_0/libs/utility/enable_if.html
//
template<class T, class U = void>
struct enable_if : public enable_if_v<T::val, U> { };
template<class T, class U = void>
struct disable_if : public enable_if<not_<T>, U> { };

struct empty { };
} // mp
} // flint

#endif
