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

#ifndef CXX_TRAITS_H
#define CXX_TRAITS_H

// only for true_/false_
#include "mp.h"

namespace flint {
namespace detail {
template<class T>
struct wrap
{
    T t;
};
} // detail

namespace traits {
///////////////////////
// BASIC TYPE TRAITS
///////////////////////
// These helpers can be used to manipulate and inquire type information.
// For example, given an arbitrary type T, one might be interested in knowing
// if it is an integral type (int, short, ulong, etc).
//
// This file contains generic traits, not specific to FLINT.

using mp::true_;
using mp::false_;

// Compute if T belongs to the signed integer types.
template<class T> struct is_signed_integer : false_ { };
template<> struct is_signed_integer<signed char> : true_ { };
template<> struct is_signed_integer<signed short> : true_ { };
template<> struct is_signed_integer<signed int> : true_ { };
template<> struct is_signed_integer<slong> : true_ { };

// Compute if T belongs to the unsigned integer types.
template<class T> struct is_unsigned_integer : false_ { };
template<> struct is_unsigned_integer<unsigned char> : true_ { };
template<> struct is_unsigned_integer<unsigned short> : true_ { };
template<> struct is_unsigned_integer<unsigned int> : true_ { };
template<> struct is_unsigned_integer<ulong> : true_ { };

// Compute if T belongs to the signed or unsigned integer types
template<class T> struct is_integer
    : mp::or_<is_unsigned_integer<T>, is_signed_integer<T> > { };

// Compute if T can always losslessly be converted into an slong
template<class T, class Enable = void> struct fits_into_slong : mp::false_ { };
template<class T>
struct fits_into_slong<T, typename mp::enable_if<traits::is_integer<T> >::type>
    : mp::or_<
          is_signed_integer<T>,
          mp::and_v<
              is_unsigned_integer<T>,
              sizeof(T) < sizeof(slong)
            >
        > { };

template<class T> struct fits_into_mp_bitcnt_t : is_unsigned_integer<T> { };

// Compute if T is like const char*
template<class T> struct is_string : mp::false_ { };
template<> struct is_string<char*> : mp::true_ { };
template<> struct is_string<const char*> : mp::true_ { };
template<int n> struct is_string<char[n]> : mp::true_ { };
template<int n> struct is_string<const char[n]> : mp::true_ { };

// Compute a type appropriate for forwarding T. This is just the appropriate
// constant reference type (but avoids things like const (int&)&, which cause
// syntax errors.
template<class T, class E = void> struct forwarding {typedef const T& type;};
template<class T> struct forwarding<T&> {typedef T& type;};
template<class T> struct forwarding<const T&> {typedef const T& type;};

template<class T>
struct forwarding<T, typename mp::enable_if<is_integer<T> >::type>
{
    typedef T type;
};
template<class T> struct forwarding<T*> {typedef T* type;};
template<> struct forwarding<bool> {typedef bool type;};

// Compute a type appropriate for referencing. Usually T&.
template<class T> struct reference {typedef T& type;};
template<class T> struct reference<T&> {typedef T& type;};
template<class T> struct reference<const T&> {typedef const T& type;};

// Add a constant qualification. In particular, turn T& into const T&.
template<class T> struct make_const {typedef const T type;};
template<class T> struct make_const<T&> {typedef const T& type;};

// Strip const and reference type annotations. This does *not* strip pointers!
template<class T> struct basetype {typedef T type;};
template<class T> struct basetype<const T> {typedef T type;};
template<class T> struct basetype<T&> {typedef T type;};
template<class T> struct basetype<const T&> {typedef T type;};

struct no { int data[2]; };
typedef int yes;
template<class T> detail::wrap<T> fakeinstance();

// use with care
template<class To, class From>
struct _is_convertible
{
private:
    static yes test(...) {return yes();}
    static no test(To) {return no();}
public:
    static const bool val = (sizeof(test(fakeinstance<From>().t)) != sizeof(yes));
};

// XXX HACK
template<class To, class A, class B>
struct _is_convertible<To, A(B)> : false_ { };
} // traits
} // flint

#endif
