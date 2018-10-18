#pragma once

#include "harmonic.hh"

class ConstrainedHarmonic : public Harmonic {
public:
  ConstrainedHarmonic(size_t levels);
  virtual ~ConstrainedHarmonic();
  virtual Point2D mapToRibbon(size_t i, const Point2D &uv) const override;
};
