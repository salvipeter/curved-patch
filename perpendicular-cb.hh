#pragma once

#include "surface.hh"

using namespace Geometry;
using Transfinite::Ribbon;
using Transfinite::Surface;

class PerpCB : public Surface {
public:
  PerpCB();
  PerpCB(const PerpCB &) = default;
  virtual ~PerpCB();
  PerpCB &operator=(const PerpCB &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};
