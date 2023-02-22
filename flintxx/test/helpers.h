/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CXX_TEST_HELPERS_H
#define CXX_TEST_HELPERS_H

#include <iostream>
#include <fstream>
#include <string>
#include <streambuf>
#include <cstdlib>

#include "flintxx/flint_classes.h"

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

// Whether or not the compiler is good enough to compile all of t-fmpzxx.cpp
// in reasonable time and space.
#ifndef HAVE_FAST_COMPILER
#ifdef __clang__
// clang 3.4 works; let's suppose all work until someone complains
#define HAVE_FAST_COMPILER 1
#elif defined(__GNUC__)
// gcc 4.7.3 is good enough, supposedly all higher ones are too
#define HAVE_FAST_COMPILER \
    (__GNUC__ > 4 || (__GNUC__ >= 4 && \
         ((__GNUC_MINOR__ == 7 && __GNUC_PATCHLEVEL__ >= 3) \
          || (__GNUC_MINOR__ > 7))))
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

template<class T>
void test_print_read(const T& t)
{
#if !defined (__WIN32) && !defined (__CYGWIN__)
    T tmp(t); tmp.set_zero();
    FILE * f = std::fopen("test_print_read", "w+");
    tassert(f);
    print(f, t+t);
    std::fflush(f);
    std::fclose(f);

    f = std::fopen("test_print_read", "r");
    typename flint::flint_classes::to_ref<T>::type tmpref(tmp);
    tassert(read(f, tmpref) > 0);
    tassert(t+t == tmp);
    tassert(std::remove("test_print_read") == 0);
#endif
}

template<class T>
void test_print_read_pretty(const T& t)
{
#if !defined (__WIN32) && !defined (__CYGWIN__)
    T tmp(t); tmp.set_zero();
    FILE * f = std::fopen("test_print_read", "w+");
    tassert(f);
    print_pretty(f, t+t, "x");
    std::fflush(f);
    std::fclose(f);

    f = std::fopen("test_print_read", "r");
    char* rvar;
    typename flint::flint_classes::to_ref<T>::type tmpref(tmp);
    tassert(read_pretty(f, tmpref, &rvar) > 0);
    tassert(t+t == tmp);
    tassert(std::remove("test_print_read") == 0);
    flint_free(rvar);
#endif
}

inline std::string readfile(const char* name)
{
    std::ifstream t(name);
    return std::string(std::istreambuf_iterator<char>(t),
                 std::istreambuf_iterator<char>());
}

template<class T>
void tassert_fprint(const T& t, const char* expected)
{
#if !defined (__WIN32) && !defined (__CYGWIN__)
    FILE * f = std::fopen("test_fprint", "w+");
    tassert(f);
    print(f, t);
    std::fflush(f);
    std::fclose(f);

    tassert(readfile("test_fprint") == expected);
    tassert(std::remove("test_fprint") == 0);
#endif
}

template<class T, class U>
void tassert_fprint_pretty(const T& t, const U& extra, const char* expected)
{
#if !defined (__WIN32) && !defined (__CYGWIN__)
    FILE * f = std::fopen("test_fprint_pretty", "w+");
    tassert(f);
    print_pretty(f, t, extra);
    std::fflush(f);
    std::fclose(f);

    tassert(readfile("test_fprint_pretty") == expected);
    tassert(std::remove("test_fprint_pretty") == 0);
#endif
}

template<class T>
void tassert_fprint_pretty(const T& t, const char* expected)
{
#if !defined (__WIN32) && !defined (__CYGWIN__)
    FILE * f = std::fopen("test_fprint_pretty", "w+");
    tassert(f);
    print_pretty(f, t);
    std::fflush(f);
    std::fclose(f);

    tassert(readfile("test_fprint_pretty") == expected);
    tassert(std::remove("test_fprint_pretty") == 0);
#endif
}


#endif
