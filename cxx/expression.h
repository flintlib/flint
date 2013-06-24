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

#ifndef CXX_EXPRESSION_H
#define CXX_EXPRESSION_H

// TODO
// * static asserts

#include <iosfwd>

#include "cxx/mp.h"
#include "cxx/traits.h"
#include "cxx/tuple.h"

namespace flint {

namespace operations {
// TODO document these
struct immediate { };
struct plus { };
} // operations

namespace rules {
// TODO document these
struct no_op
{
    template<class U>
    static void doit(const U&) {}
};

struct UNIMPLEMENTED
{
    static const bool unimplemented_marker = true;
};

// NB: this is mostly internal
template<class T, class Enable = void>
struct copy_initialization : UNIMPLEMENTED { };

template<class T, class From, class Enable = void>
struct initialization : UNIMPLEMENTED { };

template<class T, class Enable = void>
struct print : UNIMPLEMENTED { };

template<class T, class Enable = void>
struct destruction : no_op { };

template<class T, class Enable = void>
struct empty_initialization : UNIMPLEMENTED { };

template<class T, class U, class Enable = void>
struct assignment : UNIMPLEMENTED { };

template<class T, class U, class Enable = void>
struct equals : UNIMPLEMENTED { };

// If result_is_temporary is true, then the result coincides with the
// first temporary (provided these have the same type)
// Priorities 2, 1, 0 can be used to resolve conflicts.
template<
    class T,
    bool result_is_temporary,
    unsigned priority,
    class Enable = void>
struct evaluation : UNIMPLEMENTED { };
//{
//    typedef X return_t;
//    typedef Y temporaries_t; // a tuple of *pointers*
//    static void doit(const T& input, temporaries_t temps, return_t* output);
//};
} // rules

namespace traits {
template<class T>
struct is_implemented : mp::not_<_is_convertible<rules::UNIMPLEMENTED, T> > { };
} // traits

namespace mp {
template<class T, bool result_is_temporary, unsigned min_prio = 0>
struct find_evaluation
{
private:
    typedef rules::evaluation<T, result_is_temporary, 2> r2;
    typedef rules::evaluation<T, result_is_temporary, 1> r1;
    typedef rules::evaluation<T, result_is_temporary, 0> r0;

    typedef traits::is_implemented<r2> i2;
    typedef traits::is_implemented<r1> i1;
    typedef traits::is_implemented<r0> i0;

public:
    typedef typename mp::select<void, // TODO
        mp::and_c<i2, min_prio <= 2>, r2,
        mp::and_c<i1, min_prio <= 1>, r1,
        mp::and_c<i0, min_prio <= 0>, r0
      >::type type;
};
} // mp

namespace detail {
struct EXPRESSION { };
}

namespace traits {
template<class T>
struct is_expression : _is_convertible<detail::EXPRESSION, T> { };

template<class T>
struct _is_immediate_expr
    : _is_convertible<
          typename basetype<T>::type::operation_t,
          operations::immediate
        >
{ };

// Compute if T is an expression, with operation "immediate"
template<class T, class Enable = void>
struct is_immediate_expr : _is_immediate_expr<T> { };
template<class T>
struct is_immediate_expr<T,
  typename mp::enable_if<mp::not_<is_expression<T> > >::type>
    : mp::false_ { };

// Compute if T is an immediate expression, *or not an expression at all*
template<class T>
struct is_immediate
    : mp::or_<mp::not_<is_expression<T> >, is_immediate_expr<T> > { };

// Compute if T is a non-immediate expression
template<class T>
struct is_lazy_expr
    : mp::and_<is_expression<T>, mp::not_<is_immediate_expr<T> > > { };
} // traits

namespace detail {
template<class Operation, class Expr>
struct evaluation_traits
{
    typedef typename Expr::derived_t derived_t;
    typedef typename mp::find_evaluation<derived_t, false>::type rule_t;
    typedef typename mp::find_evaluation<derived_t, true>::type temp_rule_t;
    typedef typename rule_t::return_t evaluation_return_t;
    typedef evaluation_return_t evaluated_t;

    static evaluation_return_t evaluate(const derived_t& from)
    {
        evaluation_return_t res;
        evaluate_into_fresh(res, from);
        return res;
    }

    static void evaluate_into(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<typename rule_t::temporaries_t> back_t;
        typename back_t::type temps_backing;
        typename rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, 0);
        rule_t::doit(from, temps, &to);
    }

    static void evaluate_into_fresh(evaluation_return_t& to, const derived_t& from)
    {
        typedef mp::back_tuple<
            typename temp_rule_t::temporaries_t,
            evaluation_return_t
          > back_t;
        typename back_t::type temps_backing;
        typename temp_rule_t::temporaries_t temps;
        back_t::init(temps, temps_backing, &to);
        temp_rule_t::doit(from, temps, &to);
    }
};

template<class Expr>
struct evaluation_traits<operations::immediate, Expr>
{
    typedef typename Expr::derived_t derived_t;
    typedef typename Expr::derived_t evaluated_t;
    typedef evaluated_t& evaluation_return_t;

    static evaluated_t& evaluate(derived_t& d) {return d;}
    static const evaluated_t& evaluate(const derived_t& d) {return d;}

    static void evaluate_into(derived_t& to, const derived_t& from)
    {
        rules::assignment<derived_t, derived_t>::doit(to, from);
    }

    static void evaluate_into_fresh(derived_t& to, const derived_t& from)
    {
        evaluate_into(to, from);
    }
};

template<class Expr>
struct storage_traits
    : mp::if_<
          traits::is_immediate<Expr>,
          typename traits::forwarding<Expr>::type,
          Expr
        > { };
} // detail

template<class Derived, class Operation, class Data>
class expression
    : public detail::EXPRESSION
{
private:
    Data data;

protected:
    explicit expression(const Data & d) : data(d) {}

public:
    typedef detail::evaluation_traits<Operation, expression> ev_traits_t;
    typedef typename Derived::template type<Operation, Data>::result derived_t;
    typedef typename ev_traits_t::evaluated_t evaluated_t;
    typedef typename ev_traits_t::evaluation_return_t evaluation_return_t;
    typedef Data data_t;
    typedef Operation operation_t;

private:
    derived_t& downcast() {return *static_cast<derived_t*>(this);}
    const derived_t& downcast() const
    {
        return *static_cast<const derived_t*>(this);
    }

public:
    // TODO strip qualifiers?
    template<class T>
    explicit expression(const T& t,
            typename mp::enable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0)
    {
        rules::initialization<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    explicit expression(const T& t,
            typename mp::disable_if<traits::is_implemented<
                rules::initialization<derived_t, T> > >::type* = 0,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into_fresh(downcast(), t);
    }

    expression(const expression& e)
        : data(rules::copy_initialization<derived_t>::get(e.downcast()))
    {
    }

    // NB: the compiler is not allowed to eagerly instantiate!
    expression() {rules::empty_initialization<derived_t>::doit(downcast());}

    ~expression() {rules::destruction<derived_t>::doit(downcast());}

    Data& _data() {return data;}
    const Data& _data() const {return data;}

    void print(std::ostream& o) const
    {
        rules::print<evaluated_t>::doit(evaluate(), o);
    }

    typename traits::make_const<evaluation_return_t>::type evaluate() const
    {
        return ev_traits_t::evaluate(downcast());
    }
    evaluation_return_t evaluate() {return ev_traits_t::evaluate(downcast());}

    template<class T>
    void set(const T& t,
            typename mp::enable_if<traits::is_expression<T> >::type* = 0)
    {
        T::ev_traits_t::evaluate_into(downcast(), t);
    }
    template<class T>
    void set(const T& t,
            typename mp::disable_if<traits::is_expression<T> >::type* = 0)
    {
        rules::assignment<derived_t, T>::doit(downcast(), t);
    }

    template<class T>
    bool equals(const T& t,
            typename mp::enable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return equals(t.evaluate());
    }
    template<class T>
    bool equals(const T& t,
            typename mp::disable_if<traits::is_lazy_expr<T> >::type* = 0) const
    {
        return rules::equals<evaluated_t, T>::get(evaluate(), t);
    }

protected:
    template<class T, class Op>
    struct helper
    {
        typedef mp::make_tuple<
            typename detail::storage_traits<derived_t>::type,
            typename detail::storage_traits<T>::type
          > maker;
        typedef typename maker::type type;
        typedef typename Derived::template type<
            Op, typename maker::type>::result new_derived_t;

        static new_derived_t
        make(const expression& self, const T& other)
        {
            return new_derived_t(maker::make(self.downcast(), other));
        }
    };

public:
    template<class T>
    typename helper<T, operations::plus>::new_derived_t plus(const T& t) const
    {
        return helper<T, operations::plus>::make(*this, t);
    }

    template<class T>
    struct return_types
    {
        typedef typename helper<T, operations::plus>::new_derived_t plus_return_t;
    };
};

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
inline typename mp::enable_if<traits::is_expression<Expr>, std::ostream&>::type
operator<<(std::ostream& o, const Expr& e)
{
    e.print(o);
    return o;
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<traits::is_expression<Expr1>, bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e1.equals(e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::and_<
        mp::not_<traits::is_expression<Expr1> >,
        traits::is_expression<Expr2> >,
    bool>::type
operator==(const Expr1& e1, const Expr2& e2)
{
  return e2.equals(e1);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::or_<
        traits::is_expression<Expr1>,
        traits::is_expression<Expr2> >,
    bool>::type
operator!=(const Expr1& e1, const Expr2& e2)
{
  return !(e1 == e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<
    traits::is_expression<Expr1>,
    typename Expr1::template return_types<Expr2>::plus_return_t>::type
operator+(const Expr1& e1, const Expr2& e2)
{
    return e1.plus(e2);
}

template<class Expr1, class Expr2>
inline typename mp::enable_if<mp::and_<
        mp::not_<traits::is_expression<Expr1> >,
        traits::is_expression<Expr2> >,
    typename Expr2::template return_types<Expr1>::plus_return_t>::type
operator+(const Expr1& e1, const Expr2& e2)
{
  return e2.plus(e1);
}

// default rules

namespace rules {
// Copy initialization is special in that data does not necessarily have
// a default constructor.
// TODO could store pointers instead of references.
template<class T>
struct copy_initialization<T,
    typename mp::disable_if<traits::is_immediate<T> >::type>
{
    static typename T::data_t get(const T& from)
    {
        return typename T::data_t(from._data());
    }
};

template<class T>
struct copy_initialization<T,
    typename mp::enable_if<traits::is_immediate<T> >::type>
{
    static typename T::data_t get(const T& from)
    {
        T tmp;
        initialization<T, T>::doit(tmp, from);
        return tmp._data();
    }
};

namespace rdetail {
template<class T>
struct is_evaluation_2_2 : mp::false_ { };
template<class Data1, class Data2>
struct is_evaluation_2_2<tuple<Data1, tuple<Data2, empty_tuple> > >
    : mp::and_<
          traits::is_lazy_expr<Data1>,
          traits::is_lazy_expr<Data2>
        > { };
} // rdetail

template<class T, bool result_is_temporary>
struct evaluation<
    T, result_is_temporary, 0,
    typename mp::enable_if<mp::and_<
        traits::is_expression<T>,
        rdetail::is_evaluation_2_2<typename T::data_t>
      > >::type
  >
{
};
} // rules
} // flint

#endif
