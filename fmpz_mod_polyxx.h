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

#ifndef FMPZ_MOD_POLYXX_H
#define FMPZ_MOD_POLYXX_H

#include "fmpz_mod_poly.h"

#include "fmpzxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"
#include "flintxx/frandxx.h"
#include "flintxx/stdmath.h"

// TODO create/use fmpz_modxx class?

namespace flint {
FLINT_DEFINE_BINOP(fmpz_mod_polyxx_get_coeff) // TODO standardise?

template<class Operation, class Data>
class fmpz_mod_polyxx_expression
    : public expression<derived_wrapper<fmpz_mod_polyxx_expression>,
                            Operation, Data>
{
    typedef expression<derived_wrapper< ::flint::fmpz_mod_polyxx_expression>,
              Operation, Data> base_t;

    FLINTXX_DEFINE_BASICS(fmpz_mod_polyxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_mod_polyxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_mod_polyxx_expression, fmpz_mod_poly_struct, _poly)

    // these only make sense with immediates
    fmpzxx_srcref _mod() const
        {return fmpzxx_srcref::make(fmpz_mod_poly_modulus(_poly()));}

    // These only make sense with target immediates
    void realloc(slong alloc) {fmpz_mod_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {fmpz_mod_poly_fit_length(_poly(), len);}
    void _normalise() {_fmpz_mod_poly_normalise(_poly());}
    void set_coeff(slong n, ulong c) {fmpz_mod_poly_set_coeff_ui(_poly(), n, c);}
    template<class Fmpz>
    typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type
    set_coeff(slong j, const Fmpz& c)
    {
        fmpz_mod_poly_set_coeff_fmpz(_poly(), j, c.evaluate()._fmpz());
    }
    void truncate(slong n) {fmpz_mod_poly_truncate(_poly(), n);}

    // this works without evaluation
    fmpzxx_srcref modulus() const;

    evaluated_t create_temporary() const
    {
        return evaluated_t(modulus());
    }

    // These cause evaluation
    slong length() const {return fmpz_mod_poly_length(this->evaluate()._poly());}
    slong degree() const {return fmpz_mod_poly_degree(this->evaluate()._poly());}
    bool is_zero() const {return fmpz_mod_poly_is_zero(this->evaluate()._poly());}

    // Lazy members
    FLINTXX_DEFINE_MEMBER_BINOP_(get_coeff, fmpz_mod_polyxx_get_coeff)
    FLINTXX_DEFINE_MEMBER_BINOP_(operator(), compeval)
};

namespace detail {
struct fmpz_mod_poly_data;
}

typedef fmpz_mod_polyxx_expression<operations::immediate, detail::fmpz_mod_poly_data>
           fmpz_mod_polyxx;
typedef fmpz_mod_polyxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_mod_polyxx, fmpz_mod_poly_struct> >
           fmpz_mod_polyxx_ref;
typedef fmpz_mod_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_mod_polyxx, fmpz_mod_polyxx_ref, fmpz_mod_poly_struct> >
           fmpz_mod_polyxx_srcref;

#define FMPZ_MOD_POLYXX_COND_S FLINTXX_COND_S(fmpz_mod_polyxx)
#define FMPZ_MOD_POLYXX_COND_T FLINTXX_COND_T(fmpz_mod_polyxx)

namespace detail {
struct fmpz_mod_poly_data
{
    fmpz_mod_poly_t inner;
    typedef fmpz_mod_poly_t& data_ref_t;
    typedef const fmpz_mod_poly_t& data_srcref_t;

    template<class Fmpz>
    fmpz_mod_poly_data(const Fmpz& n,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpz_mod_poly_init(inner, n.evaluate()._fmpz());
    }
    template<class Fmpz>
    fmpz_mod_poly_data(const Fmpz& n, slong alloc,
            typename mp::enable_if<traits::is_fmpzxx<Fmpz> >::type* = 0)
    {
        fmpz_mod_poly_init2(inner, n.evaluate()._fmpz(), alloc);
    }
    ~fmpz_mod_poly_data() {fmpz_mod_poly_clear(inner);}

    fmpz_mod_poly_data(const fmpz_mod_poly_data& o)
    {
        fmpz_mod_poly_init(inner, fmpz_mod_poly_modulus(o.inner));
        fmpz_mod_poly_set(inner, o.inner);
    }

    fmpz_mod_poly_data(fmpz_mod_polyxx_srcref r)
    {
        fmpz_mod_poly_init(inner, r.modulus()._fmpz());
        fmpz_mod_poly_set(inner, r._poly());
    }
};

struct is_fmpz_mod_polyxx_predicate
{
    template<class T> struct type : FMPZ_MOD_POLYXX_COND_S<T> { };
};
}
template<class Operation, class Data>
inline fmpzxx_srcref
fmpz_mod_polyxx_expression<Operation, Data>::modulus() const
{
    return tools::find_subexpr<detail::is_fmpz_mod_polyxx_predicate>(
            *this)._mod();
}

namespace rules {
FLINT_DEFINE_DOIT_COND2(assignment, FMPZ_MOD_POLYXX_COND_T,
        FMPZ_MOD_POLYXX_COND_S, fmpz_mod_poly_set(to._poly(), from._poly()))

FLINTXX_DEFINE_SWAP(fmpz_mod_polyxx, fmpz_mod_poly_swap(e1._poly(), e2._poly()))

FLINTXX_DEFINE_EQUALS(fmpz_mod_polyxx, fmpz_mod_poly_equal(e1._poly(), e2._poly()))

FLINT_DEFINE_BINARY_EXPR_COND2(fmpz_mod_polyxx_get_coeff_op, fmpzxx,
        FMPZ_MOD_POLYXX_COND_S, traits::fits_into_slong,
        fmpz_mod_poly_get_coeff_fmpz(to._fmpz(), e1._poly(), e2))

FLINT_DEFINE_BINARY_EXPR_COND2(plus, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_add(to._poly(), e1._poly(), e2._poly()))
FLINT_DEFINE_BINARY_EXPR_COND2(minus, fmpz_mod_polyxx,
        FMPZ_MOD_POLYXX_COND_S, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_sub(to._poly(), e1._poly(), e2._poly()))

FLINT_DEFINE_UNARY_EXPR_COND(negate, fmpz_mod_polyxx, FMPZ_MOD_POLYXX_COND_S,
        fmpz_mod_poly_neg(to._poly(), from._poly()))
} // rules
} // flint

#endif
