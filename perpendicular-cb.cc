#include "domain-regular.hh"
#include "parameterization-bilinear.hh"
#include "ribbon-perpendicular.hh"
#include "perpendicular-cb.hh"

using DomainType = Transfinite::DomainRegular;
using ParamType = Transfinite::ParameterizationBilinear;
using RibbonType = Transfinite::RibbonPerpendicular;

PerpCB::PerpCB() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

PerpCB::~PerpCB() {
}

Point3D
PerpCB::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += cornerInterpolant(i, sds) * blends[i];
  return p;
}

std::shared_ptr<Ribbon>
PerpCB::newRibbon() const {
  return std::make_shared<RibbonType>();
}
