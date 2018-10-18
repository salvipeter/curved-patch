#pragma once

#include <parameterization.hh>

using namespace Geometry;
using Transfinite::Parameterization;

class CurvedMean : public Parameterization {
public:
  virtual ~CurvedMean();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};
