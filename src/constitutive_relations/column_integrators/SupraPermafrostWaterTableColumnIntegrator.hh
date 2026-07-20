/*
  Copyright 2010-202x held jointly by participating institutions.
  ATS is released under the three-clause BSD License.
  The terms of use and "as is" disclaimer for this license are
  provided in the top-level COPYRIGHT file.

  Authors: Devon Hill
*/

/*

Column integrator for the supra-permafrost (perched, flow-relevant) water
table -- the top of the first unfrozen AND saturated zone found scanning from
the surface downward.

This mirrors PerchedWaterTableColumnIntegrator (top-down) with two deliberate
differences:

  1. When no unfrozen-saturated cell is found in the column (e.g. the active
     layer is fully frozen, or fully dry), the result is NaN -- a "no liquid
     water table" sentinel -- rather than the full column depth.  This lets a
     winter/frozen column read as an *absence* of a water table instead of a
     spurious metres-deep value.

  2. There is no pressure-interpolation branch; the depth is always the top
     face of the accumulated (masked) cells above the zone.

The masking (unfrozen vs frozen) lives entirely in the integrator's scan()
predicate; this driver only owns the column geometry and the not-found
sentinel.

*/

#pragma once

#include "Factory.hh"
#include "EvaluatorColumnIntegrator.hh"
#include <math.h>

namespace Amanzi {
namespace Relations {

template <class Parser, class Integrator>
class SupraPermafrostWaterTableColumnIntegrator
  : public EvaluatorColumnIntegrator<Parser, Integrator> {
 public:
  using EvaluatorColumnIntegrator<Parser, Integrator>::EvaluatorColumnIntegrator;
  SupraPermafrostWaterTableColumnIntegrator(const SupraPermafrostWaterTableColumnIntegrator& other) =
    default;
  Teuchos::RCP<Evaluator> Clone() const override;

 protected:
  // Required methods from EvaluatorColumnIntegrator to overide
  virtual void Evaluate_(const State& S, const std::vector<CompositeVector*>& result) override;

  using EvaluatorColumnIntegrator<Parser, Integrator>::plist_;
  using EvaluatorColumnIntegrator<Parser, Integrator>::dependencies_;

 private:
  static Utils::RegisteredFactory<Evaluator, SupraPermafrostWaterTableColumnIntegrator<Parser, Integrator>>
    reg_;
};


template <class Parser, class Integrator>
Teuchos::RCP<Evaluator>
SupraPermafrostWaterTableColumnIntegrator<Parser, Integrator>::Clone() const
{
  return Teuchos::rcp(new SupraPermafrostWaterTableColumnIntegrator<Parser, Integrator>(*this));
}


// Required methods from EvaluatorColumnIntegrator
template <class Parser, class Integrator>
void
SupraPermafrostWaterTableColumnIntegrator<Parser, Integrator>::Evaluate_(
  const State& S,
  const std::vector<CompositeVector*>& result)
{
  // collect the dependencies and mesh, and instantiate the integrator functor
  std::vector<const Epetra_MultiVector*> deps;
  for (const auto& dep : dependencies_) {
    deps.emplace_back(
      S.Get<CompositeVector>(dep.first, dep.second).ViewComponent("cell", false).get());
  }
  auto mesh = result[0]->Mesh()->getParentMesh();
  Integrator integrator(plist_, deps, &*mesh);

  Epetra_MultiVector& res = *result[0]->ViewComponent("cell", false);

  for (int col = 0; col != res.MyLength(); ++col) {
    // for each column, loop over cells calling the integrator until stop is
    // requested or the column is complete
    AmanziGeometry::Point val(0., 0., NAN);
    auto col_cell = mesh->columns.getCells(col);
    double h_top = mesh->getCellCentroid(col_cell[0])[2]
        + mesh->getCellVolume(col_cell[0]) * integrator.coefficient(col) / 2;
    double h_end, h_half0;
    h_end = mesh->getCellCentroid(col_cell[col_cell.size() - 1])[2]; // default at bottom centroid
    h_half0 = mesh->getCellVolume(col_cell[0]) * integrator.coefficient(col) / 2 * (-1);

    for (int i = 0; i != col_cell.size(); ++i) { // loop from top down looking for the first unfrozen+saturated cell
      bool completed = integrator.scan(col, col_cell[i], val);
      if (completed) {
        h_end = mesh->getCellCentroid(col_cell[i])[2]; // first unfrozen+saturated cell centroid from top down
        break;
      }
    }

    if (std::isnan(val[2])) {
      // completed at the very first (top) cell: supra-permafrost table is at the surface
      res[0][col] = h_top - h_end + h_half0;
    } else if (val[2] == h_end) {
      // no unfrozen+saturated cell anywhere: no liquid water table (frozen or dry)
      res[0][col] = NAN;
    } else {
      // depth to the top of the unfrozen+saturated zone
      res[0][col] = integrator.coefficient(col) * val[0];
    }
  }
}

} // namespace Relations
} // namespace Amanzi
