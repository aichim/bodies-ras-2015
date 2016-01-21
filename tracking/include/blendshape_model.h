#pragma once

#include "mesh.h"


namespace af
{

template <class Scalar = double>
class BlendshapeModel
{
public:
  typedef boost::shared_ptr<BlendshapeModel<Scalar> > Ptr;
  typedef boost::shared_ptr<const BlendshapeModel<Scalar> > ConstPtr;

  BlendshapeModel ();

  bool
  load (const std::string &folder);

  bool
  save (const std::string &folder) const;

  bool
  loadBinary (const std::string &filename);

  bool
  saveBinary (const std::string &filename) const;

  void
  deepCopy (af::BlendshapeModel<Scalar> &copy) const;

  bool
  evaluate (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &coeffs,
            Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &points) const;

  void
  convertToTriMesh ();

  typename af::Mesh<Scalar>::Ptr neutral_mesh_;

  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> neutral_;
  Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> exps_;

  std::vector<std::string> exp_names_;

  int num_exps_;
  int num_pts_;

public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}

#include "impl/blendshape_model.hpp"

