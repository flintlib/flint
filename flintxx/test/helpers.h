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
#ifndef CXX_TEST_HELPERS_H
#define CXX_TEST_HELPERS_H

#include <iostream>
#include <cstdlib>

#ifndef EXIT_STATEMENT
#define EXIT_STATEMENT std::exit(1)
#endif

#define tassert(expr) do \
{ \
    if (!(expr)) \
    { \
        std::cout << "FAIL\n" __FILE__ ":" << __LINE__ << ": assertion " #expr " failed\n"; \
        EXIT_STATEMENT; \
    } \
} while (0)

// Whether or not the compiler is good enough to compile all of t-mpz.cpp
// in reasonable time and space.
#ifndef HAVE_FAST_COMPILER
#ifdef __clang__
// clang 3.4 works; let's suppose all work until someone complains
#define HAVE_FAST_COMPILER 1
#elif defined(__GNUC__)
// gcc 4.7.3 is good enough, supposedly all higher ones are too
#define HAVE_FAST_COMPILER \
    (__GNUC__ >= 4 && \
         ((__GNUC_MINOR__ == 7 && __GNUC_PATCHLEVEL__ >= 3) \
          || (__GNUC_MINOR__ > 7)))
#else
#define HAVE_FAST_COMPILER 0
#endif
#endif

// Count the number of temporaries needed to evaluate an expression of type T
template<class T>
unsigned count_temporaries(const T&)
{
    return T::ev_traits_t::rule_t::temporaries_t::len;
}

template<class T, class U>
bool typed_equals(const T&, const U&)
{
    return false;
}

template<class T>
bool typed_equals(const T& a, const T& b)
{
    return a == b;
}

#define _assert_exception(expr, type) do \
{ \
    bool exception_occurred = false; \
    try \
    { \
        expr; \
    } \
    catch(const type&) \
    { \
        exception_occurred = true; \
    } \
    if(!exception_occurred) \
    { \
        std::cout << "FAIL\n" __FILE__ ":" << __LINE__ \
            << ": expression did not cause an exception: " #expr " \n"; \
        EXIT_STATEMENT; \
    } \
} while(0)

#define assert_exception(expr) _assert_exception(expr, flint_exception)


#endif
