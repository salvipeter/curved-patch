#include "curved-cb.hh"

#include <ribbon-compatible-with-handler.hh>

#include "curved-domain.hh"
#include "curved-mean.hh"

using DomainType = CurvedDomain;
using ParamType = CurvedMean;
using RibbonType = Transfinite::RibbonCompatibleWithHandler;

CurvedCB::CurvedCB() {
  domain_ = std::make_shared<DomainType>();
  param_ = std::make_shared<ParamType>();
  param_->setDomain(domain_);
}

CurvedCB::~CurvedCB() {
}

Point3D
CurvedCB::eval(const Point2D &uv) const {
  Point2DVector sds = param_->mapToRibbons(uv);
  // return { uv[0], uv[1], sds[0][1] };
  DoubleVector blends = blendCorner(sds);
  Point3D p(0,0,0);
  for (size_t i = 0; i < n_; ++i)
    p += cornerInterpolantD(i, sds) * blends[i];
  return p;
}

std::shared_ptr<Ribbon>
CurvedCB::newRibbon() const {
  return std::make_shared<RibbonType>();
}
