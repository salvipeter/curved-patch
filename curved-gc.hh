#pragma once

#include <surface.hh>

using namespace Geometry;
using Transfinite::Ribbon;
using Transfinite::Surface;

class CurvedGC : public Surface {
public:
  CurvedGC();
  CurvedGC(const CurvedGC &) = default;
  virtual ~CurvedGC();
  CurvedGC &operator=(const CurvedGC &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};
