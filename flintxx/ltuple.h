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

#include "expression.h"
#include "flint_classes.h"
#include "tuple.h"

namespace flint {
namespace detail {
struct ltuple_c_base_t { };
struct INSTANTIATE_FROM_TUPLE { };
} // detail

namespace ltdetail {
template<class T, class Enable = void>
struct to_ref : flint_classes::to_ref<T> { };
template<class T>
struct to_ref<T, typename mp::disable_if<flint_classes::is_flint_class<T> >::type>
{
    typedef T& type;
};
template<> struct to_ref<void> {typedef void type;};
template<class T, class Enable = void>
struct to_srcref : flint_classes::to_srcref<T> { };
template<class T>
struct to_srcref<T, typename mp::disable_if<flint_classes::is_flint_class<T> >::type>
{
    typedef const T& type;
};

template<class T, class U = void, class V = void>
struct make_reftuple : mp::make_tuple<typename to_ref<T>::type,
    typename to_ref<U>::type, typename to_ref<V>::type> { };

// TODO make more generic and move to tuple.h?
template<template<class, class> class Transform, class Tuple>
struct transform_tuple
{
    typedef tuple<typename Transform<typename Tuple::head_t, void>::type,
              typename transform_tuple<
                  Transform, typename Tuple::tail_t>::type> type;
};
template<template<class, class> class Transform>
struct transform_tuple<Transform, empty_tuple>
{
    typedef empty_tuple type;
};

template<class Tuple> struct ref_tuple : transform_tuple<to_ref, Tuple> { };
template<class Tuple> struct srcref_tuple : transform_tuple<to_srcref, Tuple> { };
} // ltdetail

template<class Underlying, class Operation, class Data>
class ltuple_expression
    : public expression<derived_wrapper2<ltuple_expression, Underlying>,
                        Operation, Data>
{
public:
    typedef expression<derived_wrapper2< ::flint::ltuple_expression, Underlying>,
                Operation, Data> base_t;
    typedef detail::ltuple_c_base_t c_base_t;
    // internal
    typedef void IS_FLINT_CLASS;
    typedef void IS_LTUPLE_EXPRESSION;

    typedef Underlying underlying_t;

    ltuple_expression() {}

    template<class T>
    ltuple_expression(detail::INSTANTIATE_FROM_TUPLE i, const T& t)
        : base_t(i, t) {}

    template<class T>
    ltuple_expression& operator=(const T& t)
    {
        this->set(t);
        return *this;
    }

    // Only makes sense for reference types
    template<class T>
    static ltuple_expression _make(const T& t)
    {
        return ltuple_expression(Data::make(t));
    }

protected:
    explicit ltuple_expression(const Data& d) : base_t(d) {}

    template<class D, class O, class Da>
    friend class expression;
};

namespace detail {
template<class Underlying>
struct ltuple_data;
}

template<class Underlying>
struct make_ltuple
{
    typedef ltuple_expression<Underlying, operations::immediate,
              detail::ltuple_data<Underlying> > type;
    typedef ltuple_expression<Underlying, operations::immediate,
              flint_classes::ref_data<type, detail::ltuple_c_base_t> >
                  ref_type;
    typedef ltuple_expression<Underlying, operations::immediate,
              flint_classes::srcref_data<type, ref_type, detail::ltuple_c_base_t> >
                  srcref_type;

    template<class T> struct is_source : mp::or_<
        mp::equal_types<T, type>,
        mp::equal_types<T, ref_type>,
        mp::equal_types<T, srcref_type> > { };
    template<class T> struct is_target : mp::or_<
        mp::equal_types<T, type>,
        mp::equal_types<T, ref_type> > { };
};

namespace flint_classes {
template<class Ltuple>
struct ref_data<Ltuple, detail::ltuple_c_base_t>
{
    typedef void IS_REF_OR_CREF;
    typedef Ltuple wrapped_t;

    typedef typename ltdetail::ref_tuple<
        typename wrapped_t::underlying_t>::type underlying_t;
    underlying_t inner;

    static ref_data make(const underlying_t& i) {return ref_data(i);}

private:
    ref_data(const underlying_t& i) : inner(i) {}
};
template<class Ltuple, class Ref>
struct srcref_data<Ltuple, Ref, detail::ltuple_c_base_t>
{
    typedef void IS_REF_OR_CREF;
    typedef Ltuple wrapped_t;

    typedef typename ltdetail::srcref_tuple<
        typename wrapped_t::underlying_t>::type underlying_t;
    underlying_t inner;

    static srcref_data make(const underlying_t& i) {return srcref_data(i);}

private:
    srcref_data(const underlying_t& i) : inner(i) {}
};
} // flint_classes

namespace detail {
template<class Underlying>
struct ltuple_data
{
    Underlying inner;

    ltuple_data() {}

    template<class T>
    ltuple_data(INSTANTIATE_FROM_TUPLE, const T& t) : inner(t) {}
};
} // detail

namespace traits {
template<class Tuple, class Enable = void>
struct is_ltuple_expr : mp::false_ { };
template<class Tuple>
struct is_ltuple_expr<Tuple, typename Tuple::IS_LTUPLE_EXPRESSION>
    : mp::true_ { };

template<class Tuple, class Enable = void>
struct is_ltuple_source : mp::false_ { };
template<class Tuple>
struct is_ltuple_source<Tuple,
    typename mp::enable_if<is_ltuple_expr<Tuple> >::type>
    : make_ltuple<typename Tuple::underlying_t>::template is_source<Tuple> { };

template<class Tuple, class Enable = void>
struct is_ltuple_target : mp::false_ { };
template<class Tuple>
struct is_ltuple_target<Tuple,
    typename mp::enable_if<is_ltuple_expr<Tuple> >::type>
    : make_ltuple<typename Tuple::underlying_t>::template is_target<Tuple> { };
} // traits

namespace rules {
template<class Tuple1, class Tuple2>
struct assignment<Tuple1, Tuple2,
    typename mp::enable_if<mp::and_<
        traits::is_ltuple_source<Tuple2>,
        traits::is_ltuple_target<Tuple1> > >::type>
{
    static void doit(Tuple1& to, const Tuple2& from)
    {
        to._data().inner.set(from._data().inner);
    }
};

template<class Tuple1, class Tuple2>
struct equals<Tuple1, Tuple2,
    typename mp::enable_if<mp::and_<
        traits::is_ltuple_source<Tuple2>,
        traits::is_ltuple_source<Tuple1> > >::type>
{
    static bool get(const Tuple1& to, const Tuple2& from)
    {
        return to._data().inner.equals_elementwise(from._data().inner);
    }
};
} // rules

// TODO we would really like variadic templates / lvalue references here

template<class T>
inline typename make_ltuple<typename mp::make_tuple<T>::type>::ref_type
ltupleref(T& t)
{
    return make_ltuple<typename mp::make_tuple<T>::type>::ref_type::_make(
            ltdetail::make_reftuple<T>::make(t));
}
template<class T, class U>
inline typename make_ltuple<typename mp::make_tuple<T, U>::type>::ref_type
ltupleref(T& t, U& u)
{
    return make_ltuple<typename mp::make_tuple<T, U>::type>::ref_type::_make(
            ltdetail::make_reftuple<T, U>::make(t, u));
}
template<class T, class U, class V>
inline typename make_ltuple<typename mp::make_tuple<T, U, V>::type>::ref_type
ltupleref(T& t, U& u, V& v)
{
    return make_ltuple<typename mp::make_tuple<T, U, V>::type>::ref_type::_make(
            ltdetail::make_reftuple<T, U, V>::make(t, u, v));
}

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
            mp::make_tuple<T&, U&, V&>::make(t, u, v));
}

} // flint
