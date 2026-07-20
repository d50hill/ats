/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Devon Hill
*/

/*!

Computes the depth to the supra-permafrost water table: the top of the first
zone that is both unfrozen and saturated, scanning from the surface downward.

This is definition D4 in the water-table-definition family, and the target of
the experiment -- the flow-relevant, active-layer / talik water table that gates
lateral flow and carbon export.  Frozen cells are masked out (temperature at or
below the freezing transition), so both the ice-saturated zone and the
cryosuction gas whisper that confound the gas-onset `"water table depth`" are
excluded by construction.  When the active layer is fully frozen (no
unfrozen-saturated cell), the result is NaN -- "no liquid water table".

`"evaluator type`" = `"supra-permafrost water table depth`"

.. _supra-permafrost-water-table-depth-spec:
.. admonition:: supra-permafrost-water-table-depth-spec

    * `"transition width [K]`" ``[double]`` **0.2** Width of the freezing
      transition; a cell is unfrozen when T > 273.15 + 0.5 * width (anchored to
      `"thaw depth`" for consistency).

    KEYS:

    - `"temperature`" **SUBSURFACE_DOMAIN-temperature**
    - `"saturation of gas`" **SUBSURFACE_DOMAIN-saturation_gas**
    - `"subsurface cell volume`" **SUBSURFACE_DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN-cell_volume**


*/

#pragma once

#include "Factory.hh"
#include "SupraPermafrostWaterTableColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserSupraPermafrostWaterTableDepth {
  ParserSupraPermafrostWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorSupraPermafrostWaterTableDepth {
 public:
  IntegratorSupraPermafrostWaterTableDepth(Teuchos::ParameterList& plist,
                                           std::vector<const Epetra_MultiVector*>& deps,
                                           const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col);

 private:
  const Epetra_MultiVector* temp_;
  const Epetra_MultiVector* sat_;
  const Epetra_MultiVector* cv_;
  const Epetra_MultiVector* surf_cv_;
  const AmanziMesh::Mesh* mesh_;
  double trans_temp_;
};

using SupraPermafrostWaterTableDepthEvaluator =
  SupraPermafrostWaterTableColumnIntegrator<ParserSupraPermafrostWaterTableDepth,
                                            IntegratorSupraPermafrostWaterTableDepth>;

} //namespace Relations
} //namespace Amanzi
