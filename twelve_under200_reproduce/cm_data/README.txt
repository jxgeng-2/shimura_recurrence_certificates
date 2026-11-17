CM input and optional Magma regeneration
=======================================

For each N in 14,15,21,26,33,34,38,46,58:

  xN_cm_data.sage   exact fixed CM representatives for the supplied table
  xN_cm_under200.m  standalone print-only Magma exporter, full 3 <= |D| <= 199

The scope is negative fundamental D, with no class-number cutoff, satisfying
Kronecker(D,p) != 1 for every prime p dividing N. All h(D) classes and all
local orientations are included. The three nonidentity Atkin-Lehner lifts
are supplied in CM_AL_COORDS. Coordinates are in the basis (1,i,j,i*j),
with i^2=CM_ALGEBRA_I2 and j^2=CM_ALGEBRA_J2.

Load the supplied Sage file to reproduce the actual table representatives.
The flattened files retain the validated order, base point, AL lifts, record
ordering and indices exactly. No previous output, cache or historical file
is required. They contain no incremental or historical wrapper.

To regenerate coordinates independently, paste the corresponding .m file
into the Magma online calculator (https://magma.maths.usyd.edu.au/calc/).
It writes nothing: copy its printed output into a replacement .sage file.
Its algebra, maximal order, base point and AL lifts are fixed. However,
Magma's embedding/principal-ideal routines may select equivalent different
CM representatives or a different class enumeration. Such regenerated
output is mathematically suitable but need not reproduce the supplied
record-by-record numerical table without rematching its classes. Never mix
new representative indices with old computed image caches.

The .m files cover the complete range, not only an extension. They were
adapted from the tested exporter; this lean-package preparation checked
their range and print-only structure but did not run a fresh Magma batch.

Input contract:
  CM_DB, CM_RAMIFIED_PRIMES, CM_ABS_D_MIN=3, CM_ABS_D_MAX=199;
  CM_ALGEBRA_I2, CM_ALGEBRA_J2, CM_ORDER_BASIS, CM_BASE_COORDS,
  CM_ORIENTATION_UNIT, CM_AL_COORDS;
  CM_SELECTED_DS, CM_EXPECTED_CLASS_REPS, CM_EXPECTED_ALL_ORIENTATIONS;
  CM_CLASS_FIELDS/CM_CLASS_RECORDS and cm_points;
  CM_ALL_FIELDS/CM_ALL_RECORDS and cm_points_all_orientations.

Each all-orientation record has D,h,orbit_index,chi,mu_coords,
orientation_translate,orientation_index,class_index. The orbit_index is
(orientation_index-1)*h+class_index within each discriminant. The list is
ordered as exported (not necessarily sorted by orbit_index).

Optional light sanity check, from this directory in Sage:
  load('check_cm_inputs.sage')
  check_cm_inputs()
This checks exact coordinates, norms, local coverage and the source bounds.
It does not run Magma, construct periods, integrate CM divisors or run algdep.
