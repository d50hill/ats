/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Devon Hill
*/

/*!

Computes the depth to the water table using a liquid-saturation threshold
instead of the gas-onset test used by `"water table depth`".

This is definition D3 in the water-table-definition family: scan bottom-up for
the top of the continuously saturated zone, where "saturated" means the liquid
saturation is at or above a cutoff.  Because an ice-filled pore has low liquid
saturation, the frozen zone is excluded by construction -- unlike the gas-onset
`"water table depth`", which counts an ice-saturated pore as saturated.

`"evaluator type`" = `"liquid saturation water table depth`"

.. _liquid-saturation-water-table-depth-spec:
.. admonition:: liquid-saturation-water-table-depth-spec

    * `"saturation cutoff [-]`" ``[double]`` **0.9** Liquid saturation at or
      above which a cell is considered saturated.

    * `"interpolate depth from pressure`" ``[bool]`` **false** Default to locate
      the top face of the last continuously saturated cell from bottom upward.
      If true, interpolate the p = p_atm crossing between centroids (note: this
      interpolates the pressure datum, not the saturation crossing).

    KEYS:

    - `"saturation of liquid`" **SUBSURFACE_DOMAIN-saturation_liquid**
    - `"pressure`" **SUBSURFACE_DOMAIN-pressure**
    - `"subsurface cell volume`" **SUBSURFACE_DOMAIN-cell_volume**
    - `"surface cell volume`" **DOMAIN-cell_volume**


*/

#pragma once

#include "Factory.hh"
#include "WaterTableColumnIntegrator.hh"

namespace Amanzi {
namespace Relations {

struct ParserLiquidSaturationWaterTableDepth {
  ParserLiquidSaturationWaterTableDepth(Teuchos::ParameterList& plist, const KeyTag& key_tag);
  KeyTagSet dependencies;
};


class IntegratorLiquidSaturationWaterTableDepth {
 public:
  IntegratorLiquidSaturationWaterTableDepth(Teuchos::ParameterList& plist,
                                            std::vector<const Epetra_MultiVector*>& deps,
                                            const AmanziMesh::Mesh* mesh);
  int scan(AmanziMesh::Entity_ID col, AmanziMesh::Entity_ID c, AmanziGeometry::Point& p);
  double coefficient(AmanziMesh::Entity_ID col);

 private:
  const Epetra_MultiVector* sat_;
  const Epetra_MultiVector* pres_;
  const Epetra_MultiVector* cv_;
  const Epetra_MultiVector* surf_cv_;
  const AmanziMesh::Mesh* mesh_;
  bool is_interp_;
  double s_cut_;
};

using LiquidSaturationWaterTableDepthEvaluator =
  WaterTableColumnIntegrator<ParserLiquidSaturationWaterTableDepth,
                             IntegratorLiquidSaturationWaterTableDepth>;

} //namespace Relations
} //namespace Amanzi
