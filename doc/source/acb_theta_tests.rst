.. _acb-theta-test:

**Tests for acb_theta.h**
===============================================================================

We document the tests performed in that module. This file will later be removed
from the Arb source tree. All the tests have been run under Valgrind.

.. function:: t-eld_interval.c

Generate random *a*, *ctr* and *rad*. Work at default precision. 50% of the
time, *ctr* is set to *a* and *rad* is set to a positive number so that the
presence of at least one point in the interval is guaranteed. Compute *min*,
*mid*, *max and check that
* If the existence of a point is guaranteed, then *min* is at most *max*.
* If *min* is at most *max*, then both are congruent to *a* mod 2, and *mid* as
  well; moreover *mid* is between *min* and *max*.
* If *min* is at most *max*, then max+3-rad is greater than *ctr*, and
  min-3+rad is smaller than *ctr*.

.. function:: t-eld_points.c

Generate random integers *d, g*, Cholesky matrix *Y*, random (positive) radius
*R2*, *a* last coordinates (congruent to the bits of *a* mod 2), and
offset. Fill the associated ellipsoid *E*, and list its points. Check that:
* All ellipsoid points are within the ellipsoid box.
* All ellipsoid points have the correct last coordinates.
Then, generate random points in box; check that
* Points inside the ellipsoid must appear in the list of all points.
* Points outside the ellipsoid cannot have norm smaller than *R2*.

This indirectly tests the following functions: :func:`eld_init, eld_clear,
eld_fill, eld_points, eld_contains`.

.. function:: t-naive_radius.c

Generate random choices of precision and Cholesky matrix. Run
`acb_theta_naive_radius`, then `acb_theta_naive_tail`, and check that this is
indeed not greater than the desired error.

.. function:: t-naive_ind_const.c

Generate random *tau* in the upper half plane with `g=1`. Check that the
results, at various precisions, overlap with the results of
:func:`acb_modular_theta`.

This indirectly tests the generic theta machinery for `g=1, z=0`.

.. function:: t-naive_const.c

Generate random *tau* in the Siegel half space. Check that the results agree
with :func:`naive_ind_const` and the duplication formula.

This indirectly tests the generic theta machinery for any *g* and `z=0`, as
well as :func:`acb_theta_duplication`.

.. function:: t-naive_all_const.c

Generate random *tau* in the Siegel half space. Check that the results agree
with the duplication formula. This indirectly tests
`acb_theta_duplication_all`.

.. function:: t-naive.c


