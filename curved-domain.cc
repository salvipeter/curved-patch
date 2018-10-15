#include "curved-domain.hh"

#define ANSI_DECLARATORS
#define REAL double
#define VOID void
extern "C" {
#include "triangle.h"
}

#include "lsq-plane.hh"

CurvedDomain::CurvedDomain() {
}

CurvedDomain::~CurvedDomain() {
}

bool
CurvedDomain::update() {
  return false;
}

const Point2DVector &
CurvedDomain::parameters(size_t resolution) const {
  return parameters_;
}

TriMesh
CurvedDomain::meshTopology(size_t resolution) const {
  return {};
}
