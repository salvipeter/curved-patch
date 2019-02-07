#include "curved-cr.hh"

#include <ribbon-perpendicular.hh>
#include <utilities.hh>

#include "curved-domain.hh"
#include "harmonic.hh"

using DomainType = CurvedDomain;
using ParamType = Harmonic;
using RibbonType = Transfinite::RibbonPerpendicular;

CurvedCR::CurvedCR() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>(10); // 2^k x 2^k grid
  param_->setDomain(domain_);
}

CurvedCR::~CurvedCR() {
}

Point3D
CurvedCR::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += compositeRibbon(i, sds[i]) * (blends[i] + blends[prev(i)]);
  return p * 0.5;
}

std::shared_ptr<Ribbon>
CurvedCR::newRibbon() const {
  return std::make_shared<RibbonType>();
}

Point3D
CurvedCR::compositeRibbon(size_t i, const Point2D &sd) const {
  double s = sd[0], d = sd[1], s1 = 1.0 - s, d1 = 1.0 - d;
  double Hs = Transfinite::hermite(0, s), Hd = Transfinite::hermite(0, d), Hs1 = 1.0 - Hs;
  return sideInterpolant(i, s, d) * Hd
    + sideInterpolant(prev(i), d1, s) * Hs
    + sideInterpolant(next(i), d, s1) * Hs1
    - cornerCorrection(prev(i), d, s) * Hs * Hd
    - cornerCorrection(i, s1, d) * Hs1 * Hd;
}
