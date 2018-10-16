#pragma once

#include <domain.hh>

using namespace Geometry;
using Transfinite::Domain;

class CurvedDomain : public Domain {
public:
  CurvedDomain();
  virtual ~CurvedDomain();
  virtual bool update() override;
  virtual const Point2DVector &parameters(size_t resolution) const override;
  virtual TriMesh meshTopology(size_t resolution) const override;
  const std::vector<BSCurve> &boundaries() const;
private:
  void updateMesh(size_t resolution);

  bool mesh_updated;
  std::vector<BSCurve> plane_curves_;
  Point2DVector parameters_;
  TriMesh mesh_;
};
