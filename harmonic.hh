#pragma once

#include <parameterization.hh>

using namespace Geometry;
using Transfinite::Parameterization;

struct GridValue {
  bool boundary;
  double value;
};
using HarmonicMap = std::vector<GridValue>;

class Harmonic : public Parameterization {
public:
  Harmonic(size_t levels);
  virtual ~Harmonic();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
  virtual void update() override;
private:
  size_t levels_, size_;
  std::vector<HarmonicMap> maps_;
};
