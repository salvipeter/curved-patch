#pragma once

#include <surface.hh>

using namespace Geometry;
using Transfinite::Ribbon;
using Transfinite::Surface;

class CurvedCB : public Surface {
public:
  CurvedCB();
  CurvedCB(const CurvedCB &) = default;
  virtual ~CurvedCB();
  CurvedCB &operator=(const CurvedCB &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
};
