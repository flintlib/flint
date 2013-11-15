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

// Lazy tuple class (for use in expression templates)
// Note that assignment and comparison are performed elementwise, and types
// need not match (even for equality) as long as the operations can be performed
// on underlying types.

#ifndef FLINTXX_LTUPLE_H
#define FLINTXX_LTUPLE_H

#ifndef FLINT_LTUPLE_PLACEHOLDER_NAME
#define FLINT_LTUPLE_PLACEHOLDER_NAME _
#endif

#include "expression.h"
#include "tuple.h"

namespace flint {
// For lazy get<n>, this operation type is created.
namespace operations {
template<unsigned n> struct ltuple_get_op { };
} // operations

namespace detail {
// Empty marker type
struct INSTANTIATE_FROM_TUPLE { };

// Traits for the ltuple expression template get<> operation. If the ltuple is
// an immediate, return references. Else if the return type is an expression
// template, return an expression template. Otherwise, evaluate the ltuple and
// return a copy of the entry.
template<unsigned n, class Underlying, class Operation, class Data, class Expr,
    class Enable = void>
struct ltuple_get_traits
{
    typedef unary_op_helper<operations::ltuple_get_op<n>, Expr> uoh;
    typedef typename uoh::return_t type;
    typedef type ctype;
    static type get(const Expr& e)
    {
        return uoh::make(e);
    }
};

template<unsigned n, class Underlying, class Data, class Expr>
struct ltuple_get_traits<n, Underlying, operations::immediate, Data, Expr>
{
    typedef mp::tuple_get<Underlying, n> getter;
    typedef typename getter::type btype;
    typedef typename traits::forwarding<btype>::type ctype;
    typedef typename traits::reference<btype>::type type;

    static type get(Expr& t)
    {
        return getter::get(t._data().inner);
    }
    static ctype get(const Expr& t)
    {
        return getter::get(t._data().inner);
    }
};

template<unsigned n, class Underlying, class Operation, class Data, class Expr>
struct ltuple_get_traits<n, Underlying, Operation, Data, Expr,
    typename mp::enable_if<mp::and_<
        mp::not_<mp::equal_types<Operation, operations::immediate> >,
        mp::not_<traits::is_expression<typename mp::tuple_get<Underlying, n>::type> >
      > >::type>
{
    typedef mp::tuple_get<Underlying, n> getter;
    typedef typename getter::type type;
    typedef type ctype;

    static type get(const Expr& t)
    {
        return getter::get(t.evaluate()._data().inner);
    }
};

// Instances of this can be passed to ltuple[ref]() and will be replaced by
// temporaries of the right type before assigment.
struct IGNORED_TYPE { };
template<class Ltuple, class To, class Enable = void>
struct ltuple_instantiate_ignored_types;
} // detail

// The ltuple expression template class. Underlying is a tuple type.
template<class Underlying, class Operation, class Data>
class ltuple_expression
    : public expression<derived_wrapper2<ltuple_expression, Underlying>,
                        Operation, Data>
{
public:
    typedef expression<derived_wrapper2< ::flint::ltuple_expression, Underlying>,
                Operation, Data> base_t;
    // internal
    typedef void IS_LTUPLE_EXPRESSION;

    typedef Underlying underlying_t;

    ltuple_expression() {}

    template<class T>
    ltuple_expression(detail::INSTANTIATE_FROM_TUPLE i, const T& t)
        : base_t(i, t) {}

    template<class T>
    ltuple_expression& operator=(const T& t)
    {
        detail::ltuple_instantiate_ignored_types<
            ltuple_expression, T> inst(*this, t);
        inst.set(t);
        return *this;
    }

    template<unsigned n>
    typename detail::ltuple_get_traits<n, Underlying,
             Operation, Data, ltuple_expression>::type get()
    {
        return detail::ltuple_get_traits<
            n, Underlying, Operation, Data, ltuple_expression>::get(*this);
    }
    template<unsigned n>
    typename detail::ltuple_get_traits<n, Underlying,
             Operation, Data, ltuple_expression>::ctype
    get() const
    {
        return detail::ltuple_get_traits<
            n, Underlying, Operation, Data, ltuple_expression>::get(*this);
    }

    typename base_t::evaluated_t create_temporary() const
    {
        return typename base_t::evaluated_t(detail::INSTANTIATE_FROM_TUPLE(),
                mp::htuples::fill<underlying_t>(tools::temporaries_filler(*this)));
    }

protected:
    explicit ltuple_expression(const Data& d) : base_t(d) {}

    template<class D, class O, class Da>
    friend class expression;
};

namespace detail {
template<class Underlying>
struct ltuple_data
{
    Underlying inner;

    ltuple_data() {}

    template<class T>
    ltuple_data(INSTANTIATE_FROM_TUPLE, const T& t) : inner(t) {}
};

template<class T> struct to_ref : traits::reference<T> { };
template<class T> struct to_srcref
    : traits::reference<typename traits::make_const<T>::type> { };

template<template<class>class Transform, class Tuple>
struct transform_tuple
{
    typedef tuple<typename Transform<typename Tuple::head_t>::type,
              typename transform_tuple<Transform, typename Tuple::tail_t>::type>
                  type;
};
template<template<class>class Transform>
struct transform_tuple<Transform, empty_tuple>
{
    typedef empty_tuple type;
};
} // detail

// Helper for building ltuple types.
template<class Underlying>
struct make_ltuple
{
    typedef ltuple_expression<Underlying, operations::immediate,
              detail::ltuple_data<Underlying> > type;

    typedef typename detail::transform_tuple<detail::to_ref, Underlying>::type
        Underlying_ref;
    typedef typename detail::transform_tuple<detail::to_srcref, Underlying>::type
        Underlying_srcref;

    typedef ltuple_expression<Underlying_ref, operations::immediate,
              detail::ltuple_data<Underlying_ref> > ref_type;
    typedef ltuple_expression<Underlying_srcref, operations::immediate,
              detail::ltuple_data<Underlying_srcref> > srcref_type;
};

namespace traits {
template<class Tuple, class Enable = void>
struct is_ltuple_expr : mp::false_ { };
template<class Tuple>
struct is_ltuple_expr<Tuple, typename Tuple::IS_LTUPLE_EXPRESSION>
    : mp::true_ { };

// enable evaluation directly into tuple
template<class To, class From>
struct can_evaluate_into<To, From,
    typename mp::enable_if<mp::and_<is_ltuple_expr<From>,
        is_ltuple_expr<To>, mp::not_<mp::equal_types<To, From> >
      > >::type> : mp::true_ { };
} // traits

namespace detail {
template<class Ltuple, class To, class Enable>
struct ltuple_instantiate_ignored_types
{
    // degenerate case: To is not a tuple
    Ltuple& saved;
    ltuple_instantiate_ignored_types(Ltuple& s, const To&) : saved(s) {}
    void set(const To& to) {saved.set(to);}
};
template<class ToTuple, class FromTuple>
struct tuple_instantiate_ignored
{
    typedef tuple_instantiate_ignored<typename ToTuple::tail_t,
                typename FromTuple::tail_t> next_t;
    typedef typename traits::reference<typename ToTuple::head_t>::type ref_t;
    typedef tuple<ref_t, typename next_t::type> type;
    next_t next;
    ref_t ref;

    template<class From>
    tuple_instantiate_ignored(ToTuple& to, const From& from)
        : next(to.tail, from), ref(to.head) {}

    type get()
    {
        return type(ref, next.get());
    }
};
template<>
struct tuple_instantiate_ignored<empty_tuple, empty_tuple>
{
    typedef empty_tuple type;
    template<class F>
    tuple_instantiate_ignored(empty_tuple, const F&) {}
    empty_tuple get() {return empty_tuple();}
};
template<class ToTail, class From, class FromTail>
struct tuple_instantiate_ignored<
    tuple<IGNORED_TYPE&, ToTail>, tuple<From, FromTail> >
{
    typedef tuple_instantiate_ignored<ToTail, FromTail> next_t;
    typedef typename traits::reference<From>::type ref_t;
    typedef tuple<ref_t, typename next_t::type> type;
    next_t next;
    From tmp;

    template<class ToTuple, class FromExpr>
    tuple_instantiate_ignored(ToTuple& to, const FromExpr& from)
        : next(to.tail, from),
          tmp(rules::instantiate_temporaries<FromExpr, From>::get(from)) {}

    type get()
    {
        return type(tmp, next.get());
    }
};
template<class Ltuple, class T>
struct ltuple_instantiate_ignored_types<Ltuple, T,
    typename mp::enable_if<traits::is_ltuple_expr<T> >::type>
{
    typedef tuple_instantiate_ignored<
        typename Ltuple::underlying_t, typename T::underlying_t> tii_t;
    tii_t tii;
    ltuple_instantiate_ignored_types(Ltuple& l, const T& t)
        : tii(l._data().inner, t) {}
    void set(const T& t)
    {
        typename make_ltuple<typename tii_t::type>::type(INSTANTIATE_FROM_TUPLE(),
            tii.get()).set(t);
    }
};
}

namespace rules {
template<class Tuple1, class Tuple2>
struct assignment<Tuple1, Tuple2,
    typename mp::enable_if<mp::and_<
        traits::is_ltuple_expr<Tuple2>,
        traits::is_ltuple_expr<Tuple1> > >::type>
{
    static void doit(Tuple1& to, const Tuple2& from)
    {
        to._data().inner.set(from._data().inner);
    }
};

template<class Tuple1, class Tuple2>
struct equals<Tuple1, Tuple2,
    typename mp::enable_if<mp::and_<
        traits::is_ltuple_expr<Tuple2>,
        traits::is_ltuple_expr<Tuple1> > >::type>
{
    static bool get(const Tuple1& to, const Tuple2& from)
    {
        return to._data().inner.equals_elementwise(from._data().inner);
    }
};

template<class Tuple, unsigned n>
struct unary_expression<operations::ltuple_get_op<n>, Tuple,
    typename mp::enable_if<mp::and_<
        traits::is_ltuple_expr<Tuple>,
        traits::is_immediate<Tuple> > >::type>
{
    typedef typename mp::tuple_get<typename Tuple::underlying_t, n>::type
        return_t;
    template<class R>
    static void doit(R& to, const Tuple& from)
    {
        to = from.template get<n>();
    }
};
} // rules

// TODO we would really like variadic templates / lvalue references here

// Helpers to build ltuples from (references to) arguments.

template<class T>
inline typename make_ltuple<typename mp::make_tuple<T>::type>::type
ltuple(const T& t)
{
    return typename make_ltuple<typename mp::make_tuple<T>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T>::make(t));
}
template<class T, class U>
inline typename make_ltuple<typename mp::make_tuple<T, U>::type>::type
ltuple(const T& t, const U& u)
{
    return typename make_ltuple<typename mp::make_tuple<T, U>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T, U>::make(t, u));
}

template<class T, class U, class V>
inline typename make_ltuple<typename mp::make_tuple<T, U, V>::type>::type
ltuple(const T& t, const U& u, const V& v)
{
    return typename make_ltuple<typename mp::make_tuple<T, U, V>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T, U, V>::make(t, u, v));
}
template<class T, class U, class V, class W>
inline typename make_ltuple<typename mp::make_tuple<T, U, V, W>::type>::type
ltuple(const T& t, const U& u, const V& v, const W& w)
{
    return typename make_ltuple<typename mp::make_tuple<T, U, V, W>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T, U, V, W>::make(t, u, v, w));
}
template<class T>
inline typename make_ltuple<typename mp::make_tuple<T&>::type>::type
ltupleref(T& t)
{
    return typename make_ltuple<typename mp::make_tuple<T&>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T&>::make(t));
}
template<class T, class U>
inline typename make_ltuple<typename mp::make_tuple<T&, U&>::type>::type
ltupleref(T& t, U& u)
{
    return typename make_ltuple<typename mp::make_tuple<T&, U&>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T&, U&>::make(t, u));
}
template<class T, class U, class V>
inline typename make_ltuple<typename mp::make_tuple<T&, U&, V&>::type>::type
ltupleref(T& t, U& u, V& v)
{
    return typename make_ltuple<typename mp::make_tuple<T&, U&, V&>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T&, U&, V&>::make(t, u, v));
}
template<class T, class U, class V, class W>
inline typename make_ltuple<typename mp::make_tuple<T&, U&, V&, W&>::type>::type
ltupleref(T& t, U& u, V& v, W& w)
{
    return typename make_ltuple<typename mp::make_tuple<T&, U&, V&, W&>::type>::type(
            detail::INSTANTIATE_FROM_TUPLE(),
            mp::make_tuple<T&, U&, V&, W&>::make(t, u, v, w));
}

// static placeholder
static detail::IGNORED_TYPE FLINT_LTUPLE_PLACEHOLDER_NAME;

namespace detail {
void remove_compiler_warning(
        detail::IGNORED_TYPE* = &FLINT_LTUPLE_PLACEHOLDER_NAME);
} // detail
} // flint

#endif
