.. _nf:

**nf.h** -- number fields
========================================================================================

.. type:: nf_struct

.. type:: nf_t

    Represents a number field.

.. function:: void nf_init(nf_t nf, const fmpq_poly_t pol)

    Perform basic initialisation of a number field (for element arithmetic)
    given a defining polynomial over `\mathbb{Q}`. 

.. function:: void nf_clear(nf_t nf)

    Release resources used by a number field object. The object will need
    initialisation again before it can be used.

