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

#ifndef FMPZ_POLYXX_H
#define FMPZ_POLYXX_H

#include "fmpz_poly.h"

#include "fmpzxx.h"

#include "flintxx/expression.h"
#include "flintxx/flint_classes.h"

// TODO exhibit this as a specialisation of a generic poly<fmpzxx>

namespace flint {
namespace detail {
template<class Poly>
struct fmpz_poly_traits
{
    struct coeff_ref_t { };
    struct coeff_srcref_t { };
};
}

template<class Operation, class Data>
class fmpz_polyxx_expression
    : public expression<derived_wrapper<fmpz_polyxx_expression>,
                            Operation, Data>
{
public:
    typedef expression<derived_wrapper< ::flint::fmpz_polyxx_expression>,
              Operation, Data> base_t;
    typedef detail::fmpz_poly_traits<fmpz_polyxx_expression> poly_traits_t;
    typedef typename poly_traits_t::coeff_ref_t coeff_ref_t;
    typedef typename poly_traits_t::coeff_srcref_t coeff_srcref_t;

    FLINTXX_DEFINE_BASICS(fmpz_polyxx_expression)
    FLINTXX_DEFINE_CTORS(fmpz_polyxx_expression)
    FLINTXX_DEFINE_C_REF(fmpz_polyxx_expression, fmpz_poly_struct, _poly)

    // These only make sense with immediates
    void realloc(slong alloc) {fmpz_poly_realloc(_poly(), alloc);}
    void fit_length(slong len) {fmpz_poly_fit_length(_poly(), len);}
    void _normalise() {_fmpz_poly_normalise(_poly());}
    void _set_length(slong len) {_fmpz_poly_set_length(_poly(), len);}
    void zero_coeffs(slong i, slong j) {fmpz_poly_zero_coeffs(_poly(), i, j);}

    // The result of these are undefined if n is >= length
    // You also may have to call _normalise().
    coeff_ref_t get_coeff(slong n)
    {
        return coeff_ref_t::make(fmpz_poly_get_coeff_ptr(_poly(), n));
    }
    coeff_srcref_t get_coeff(slong n) const
    {
        return coeff_srcref_t::make(fmpz_poly_get_coeff_ptr(_poly(), n));
    }
    coeff_ref_t lead()
    {
        return coeff_ref_t::make(fmpz_poly_lead(_poly()));
    }
    coeff_srcref_t lead() const
    {
        return coeff_srcref_t::make(fmpz_poly_lead(_poly()));
    }

    // These only make sense with target immediates
    template<class Fmpz>
    typename mp::enable_if<FMPZXX_COND_S<Fmpz> >::type
    set_coeff(slong n, const Fmpz& x)
    {
        fmpz_poly_set_coeff_fmpz(_poly(), n, x._fmpz());
    }
    template<class T>
    typename mp::enable_if<traits::is_signed_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpz_poly_set_coeff_si(_poly(), n, x);
    }
    template<class T>
    typename mp::enable_if<traits::is_unsigned_integer<T> >::type
    set_coeff(slong n, T x)
    {
        fmpz_poly_set_coeff_ui(_poly(), n, x);
    }

    // These cause evaluation
    slong length() const {return fmpz_poly_length(this->evaluate()._poly());}
    slong degree() const {return fmpz_poly_degree(this->evaluate()._poly());}
    bool is_one() const {return fmpz_poly_is_one(this->evaluate()._poly());}
    bool is_zero() const {return fmpz_poly_is_zero(this->evaluate()._poly());}
    bool is_unit() const {return fmpz_poly_is_unit(this->evaluate()._poly());}
    // TODO get_coeff with copies?
};

namespace detail {
struct fmpz_poly_data;
}

typedef fmpz_polyxx_expression<operations::immediate, detail::fmpz_poly_data>
           fmpz_polyxx;
typedef fmpz_polyxx_expression<operations::immediate,
            flint_classes::ref_data<fmpz_polyxx, fmpz_poly_struct> >
           fmpz_polyxx_ref;
typedef fmpz_polyxx_expression<operations::immediate,
            flint_classes::srcref_data<
                fmpz_polyxx, fmpz_polyxx_ref, fmpz_poly_struct> >
           fmpz_polyxx_srcref;

namespace detail {
template<>
struct fmpz_poly_traits<fmpz_polyxx>
{
    typedef fmpzxx_ref coeff_ref_t;
    typedef fmpzxx_srcref coeff_srcref_t;
};
template<>
struct fmpz_poly_traits<fmpz_polyxx_ref>
{
    typedef fmpzxx_ref coeff_ref_t;
    typedef fmpzxx_srcref coeff_srcref_t;
};
template<>
struct fmpz_poly_traits<fmpz_polyxx_srcref>
{
    typedef fmpzxx_srcref coeff_ref_t;
    typedef fmpzxx_srcref coeff_srcref_t;
};

struct fmpz_poly_data
{
    fmpz_poly_t inner;
    typedef fmpz_poly_t& data_ref_t;
    typedef const fmpz_poly_t& data_srcref_t;

    fmpz_poly_data() {fmpz_poly_init(inner);}
    ~fmpz_poly_data() {fmpz_poly_clear(inner);}

    fmpz_poly_data(const fmpz_poly_data& o)
    {
        fmpz_poly_init(inner);
        fmpz_poly_set(inner, o.inner);
    }

    fmpz_poly_data(fmpz_polyxx_srcref r)
    {
        fmpz_poly_init(inner);
        fmpz_poly_set(inner, r._poly());
    }

    fmpz_poly_data(slong alloc)
    {
        fmpz_poly_init2(inner, alloc);
    }
};
} // detail

namespace rules {
#define FMPZ_POLYXX_COND_S FLINTXX_COND_S(fmpz_polyxx)
#define FMPZ_POLYXX_COND_T FLINTXX_COND_T(fmpz_polyxx)
FLINTXX_DEFINE_EQUALS(fmpz_polyxx, fmpz_poly_equal(e1._poly(), e2._poly()))
} // rules
} // flint

#endif
