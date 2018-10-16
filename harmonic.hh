#pragma once

#include <parameterization.hh>

using namespace Geometry;
using Transfinite::Parameterization;

struct GridValue {
  bool boundary;
  double value;
};
using HarmonicMap = std::vector<GridValue>;

class ConstrainedHarmonic : public Parameterization {
public:
  ConstrainedHarmonic(size_t levels);
  virtual ~ConstrainedHarmonic();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
  virtual void update() override;
private:
  Point2D harmonicMap(size_t i, const Point2D &uv) const;
  size_t levels_, size_;
  std::vector<HarmonicMap> maps_;
};
