#pragma once

#include <surface.hh>

using namespace Geometry;
using Transfinite::Ribbon;
using Transfinite::Surface;

class CurvedCR : public Surface {
public:
  CurvedCR();
  CurvedCR(const CurvedCR &) = default;
  virtual ~CurvedCR();
  CurvedCR &operator=(const CurvedCR &) = default;
  virtual Point3D eval(const Point2D &uv) const override;
  using Surface::eval;

protected:
  virtual std::shared_ptr<Ribbon> newRibbon() const override;
  Point3D compositeRibbon(size_t i, const Point2D &sd) const;
};
