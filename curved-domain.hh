#pragma once

#include "domain.hh"

using namespace Geometry;
using Transfinite::Domain;

class CurvedDomain : public Domain {
public:
  CurvedDomain();
  virtual ~CurvedDomain();
  virtual bool update() override;
  virtual const Point2DVector &parameters(size_t resolution) const override;
  virtual TriMesh meshTopology(size_t resolution) const override;
private:
  Point2DVector parameters_;
};
