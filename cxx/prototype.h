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

// TODO
// * formatting
// * static asserts
// * evaluation

#ifndef CXX_PROTOTYPE_H
#define CXX_PROTOTYPE_H

#include <iosfwd>
#include <cstdlib>

namespace flint {

namespace operations {
struct immediate {};
struct plus {};
} // operations

namespace traits {
template<class T>
struct no_op
{
    template<class U>
    static void doit(const U&) {}
};

struct UNIMPLEMENTED
{
    static const bool unimplemented_marker = true;
};

template<class T, class From>
struct initialization
    : UNIMPLEMENTED
{ };

template<class T>
struct print
    : UNIMPLEMENTED
{ };

template<class T>
struct destruction
    : public no_op<T>
{};
} // traits

namespace mp {
struct no { int data[2]; };
typedef int yes;
template<class T> T fakeinstance();

// use with care
template<class To, class From>
struct _is_convertible
{
private:
    static yes test(...) {return yes();}
    static no test(To) {return no();}
public:
    static const bool val = (sizeof(test(fakeinstance<From>())) != sizeof(yes));
};

template<class T>
struct not_
{
    static const bool val = !T::val;
};

template<class T>
struct is_implemented
    : public not_<_is_convertible<traits::UNIMPLEMENTED, T> >
{ };

template<bool, class U = void>
struct enable_if_v
{
    typedef U type;
    static const int val = 0;
};
template<class U>
struct enable_if_v<false, U>
{ };

template<class T, class U = void>
struct enable_if
    : public enable_if_v<T::val, U>
{ };
}

namespace detail {
struct EXPRESSION
{ };
}

template<class Derived, class Operation, class Data>
class expression
    : public detail::EXPRESSION
{
protected:
    Data data;
    typedef typename Derived::template type<Operation, Data>::result derived_t;
    typedef derived_t evaluated_t;

private:
    derived_t& downcast() {return *static_cast<derived_t*>(this);}
    const derived_t& downcast() const {return *static_cast<const derived_t*>(this);}

public:
    explicit expression(const Data & d) : data(d) {}

    // TODO strip qualifiers?
    template<class T>
    explicit expression(const T& t, int enable = mp::enable_if<mp::is_implemented<traits::initialization<derived_t, T> > >::val)
    {
        traits::initialization<derived_t, T>::doit(downcast(), t);
    }

    ~expression() {traits::destruction<derived_t>::doit(downcast());}

    Data& _data() {return data;}
    const Data& _data() const {return data;}

    // XXX move this stuff below the traits?
    void print(std::ostream& o) const
    {
        traits::print<evaluated_t>::doit(evaluate(), o);
    }

    const evaluated_t& evaluate() const {return downcast(); /* TODO */ }
};

namespace mp {
template<class T>
struct is_expression
    : public _is_convertible<detail::EXPRESSION, T>
{ };
}

template<template<class O, class D> class Derived>
struct derived_wrapper
{
    template<class Operation, class Data>
    struct type
    {
        typedef Derived<Operation, Data> result;
    };
};


// operators

template<class Expr>
typename mp::enable_if<mp::is_expression<Expr>, std::ostream&>::type
operator<<(std::ostream& o, const Expr& e)
{
    e.print(o);
    return o;
}


// "concrete" expression classes

template<class Operation, class Data>
class mpz_expression
    : public expression<derived_wrapper<mpz_expression>, Operation, Data>
{
public:
    template<class T>
    explicit mpz_expression(const T& t) : mpz_expression::expression(t) {}
};

} // flint

#include "fmpz.h"

namespace flint {
typedef mpz_expression<operations::immediate, fmpz_t> mpz;

namespace traits {
template<>
struct initialization<mpz, mpz>
{
    static void doit(mpz& target, const mpz& source)
    {
        fmpz_init_set(target._data(), source._data());
    }
};

template<>
struct initialization<mpz, unsigned long>
{
    static void doit(mpz& target, unsigned long source)
    {
        fmpz_init_set_ui(target._data(), source);
    }
};

template<>
struct destruction<mpz>
{
    static void doit(mpz& v)
    {
        fmpz_clear(v._data());
    }
};

template<>
struct print<mpz>
{
    static void doit(const mpz& v, std::ostream& o)
    {
        char* str = fmpz_get_str(0, 10, v._data());
        o << str;
        std::free(str);
    }
};
} // traits
} // flint

#endif
