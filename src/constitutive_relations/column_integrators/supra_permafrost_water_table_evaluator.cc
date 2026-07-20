/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Devon Hill
*/

//! Calculates the supra-permafrost (unfrozen and saturated) water table depth (D4).
#include "supra_permafrost_water_table_evaluator.hh"

namespace Amanzi {
namespace Relations {

ParserSupraPermafrostWaterTableDepth::ParserSupraPermafrostWaterTableDepth(
  Teuchos::ParameterList& plist,
  const KeyTag& key_tag)
{
  Key domain = Keys::getDomain(key_tag.first);
  Tag tag = key_tag.second;

  Key domain_ss = Keys::readDomainHint(plist, domain, "surface", "subsurface");
  Key temp_key = Keys::readKey(plist, domain_ss, "temperature", "temperature");
  dependencies.insert(KeyTag{ temp_key, tag });

  Key sat_key = Keys::readKey(plist, domain_ss, "saturation of gas", "saturation_gas");
  dependencies.insert(KeyTag{ sat_key, tag });

  Key cv_key = Keys::readKey(plist, domain_ss, "subsurface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_key, tag });

  Key cv_surf_key = Keys::readKey(plist, domain, "surface cell volume", "cell_volume");
  dependencies.insert(KeyTag{ cv_surf_key, tag });
}


IntegratorSupraPermafrostWaterTableDepth::IntegratorSupraPermafrostWaterTableDepth(
  Teuchos::ParameterList& plist,
  std::vector<const Epetra_MultiVector*>& deps,
  const AmanziMesh::Mesh* mesh)
  : mesh_(mesh)
{
  AMANZI_ASSERT(deps.size() == 4);
  temp_ = deps[0];
  sat_ = deps[1];
  cv_ = deps[2];
  surf_cv_ = deps[3];

  double trans_width = plist.get<double>("transition width [K]", 0.2);
  trans_temp_ = 273.15 + 0.5 * trans_width;
}

int
IntegratorSupraPermafrostWaterTableDepth::scan(AmanziMesh::Entity_ID col,
                                               AmanziMesh::Entity_ID c,
                                               AmanziGeometry::Point& p)
{
  bool unfrozen = (*temp_)[0][c] > trans_temp_;
  bool saturated = (*sat_)[0][c] == 0.0;
  if (!(unfrozen && saturated)) {
    // frozen, or unfrozen-but-unsaturated: not the supra-permafrost water table.
    // Keep scanning down and accumulate the thickness above the zone.
    p[2] = mesh_->getCellCentroid(c)[2]; // last non-(unfrozen+saturated) cell centroid
    p[0] += (*cv_)[0][c];                // cumulative volume above the zone
    return false;
  }
  // first unfrozen + saturated cell from the top: top of the supra-permafrost zone
  return true;
}

double
IntegratorSupraPermafrostWaterTableDepth::coefficient(AmanziMesh::Entity_ID col)
{
  return 1. / (*surf_cv_)[0][col];
}

} // namespace Relations
} // namespace Amanzi
