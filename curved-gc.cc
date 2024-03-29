#include "curved-gc.hh"

#include <ribbon-perpendicular.hh>

#include "curved-domain.hh"
#include "constrained-harmonic.hh"

using DomainType = CurvedDomain;
using ParamType = ConstrainedHarmonic;
using RibbonType = Transfinite::RibbonPerpendicular;

CurvedGC::CurvedGC() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>(10); // 2^k x 2^k grid
  param_->setDomain(domain_);
}

CurvedGC::~CurvedGC() {
}

Point3D
CurvedGC::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i) {
    double s = sds[i][0], d = sds[i][1], s1 = sds[next(i)][0];
    p += sideInterpolant(i, s, d) * (blends[i] + blends[prev(i)]);
    p -= cornerCorrection(i, 1.0 - s, s1) * blends[i];
  }
  return p;
}

std::shared_ptr<Ribbon>
CurvedGC::newRibbon() const {
  return std::make_shared<RibbonType>();
}
