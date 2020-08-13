/*
    Copyright (C) 2013 Tom Bachmann

    This file is part of FLINT.

    FLINT is free software: you can redistribute it and/or modify it under
    the terms of the GNU Lesser General Public License (LGPL) as published
    by the Free Software Foundation; either version 2.1 of the License, or
    (at your option) any later version.  See <https://www.gnu.org/licenses/>.
*/

#ifndef CXX_FRANDXX_H
#define CXX_FRANDXX_H

#include "../flint.h"

// This class contains a first-class wrapper of flint_rand_t.
// Note that frandxx is not copyable.

namespace flint {
class frandxx
{
private:
    flint_rand_t inner;

    // not copyable
    frandxx(const frandxx&);

public:
    frandxx() {flint_randinit(inner);}
    ~frandxx() {flint_randclear(inner);}

    flint_rand_t& _data() {return inner;}
    const flint_rand_t& _data() const {return inner;}
};
} // flint

#endif
