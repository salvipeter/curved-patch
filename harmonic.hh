#pragma once

#include "parameterization.hh"

using namespace Geometry;
using Transfinite::Parameterization;

class ConstrainedHarmonic : public Parameterization {
public:
  ConstrainedHarmonic();
  virtual ~ConstrainedHarmonic();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
  virtual void update() override;
};
