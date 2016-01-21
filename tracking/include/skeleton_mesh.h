#pragma once

#include "skeleton.h"
#include "mesh.h"
#include "blendshape_model.h"

namespace af
{

class SkeletonMesh
{
public:
  af::Mesh<float>::Ptr mesh_, mesh_rest_;
  Skeleton skeleton_;
  Skeleton skeleton_rest_;

  SkeletonMesh ();

  void
  setRestMesh (af::Mesh<float>::ConstPtr mesh_rest)
  {
    mesh_rest_.reset (new af::Mesh<float> (*mesh_rest));
    mesh_.reset (new af::Mesh<float> (*mesh_rest));
  }

  void
  setRestSkeleton (Skeleton &skel)
  {
    skeleton_rest_ = skel;
    skeleton_ = skel;
  }

  void
  bakeGeometry ();

  void
  initIndices ();

  float
  IKWithPointToPlane (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                      pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                      const pcl::Correspondences &corresps,
                      const float weight_reg_zero = 5e-1);  ///1. for 10k corresps

  float
  IKWithPointToPoint (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                      const pcl::Correspondences &corresps,
                      const float weight_reg_zero = 5e-1);  ///1. for 10k corresps

  float
  IKCustom (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
            pcl::PointCloud<pcl::Normal>::ConstPtr normals,
            const pcl::Correspondences &corresps,
            const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &contour_normals,
            const std::vector<Eigen::Vector3f> &joint_pos,
            const Eigen::VectorXf &vars_frame_prev,
            const af::PosePCA &pose_pca,
            Eigen::VectorXf &vars_result,
            const float weight_p2plane = 1.,
            const float weight_contour = 0.5,
            const float weight_reg_zero = 5e-1,
            const float weight_close_prev = 1e-1,
            const float weight_nite = 1e-2,
            const float weight_pca_proj = 1.,
            const float weight_pca_dev = 5e-2);

  float
  IKCustomFixed (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                 const pcl::Correspondences &corresps,
                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &contour_normals,
                 const std::vector<Eigen::Vector3f> &joint_pos,
                 const Eigen::VectorXf &vars_frame_prev,
                 const af::PosePCA &pose_pca,
                 Eigen::VectorXf &vars_result,
                 const float weight_p2plane = 1.,
                 const float weight_contour = 0.5,
                 const float weight_reg_zero = 5e-1,
                 const float weight_close_prev = 1e-1,
                 const float weight_nite = 1e-2,
                 const float weight_pca_proj = 1.,
                 const float weight_pca_dev = 5e-2);

  float
  IKCustomFelix (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                 const pcl::Correspondences &corresps,
                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                 const std::vector<Eigen::Vector3f> &nite_joint_pos,
                 const Eigen::VectorXf &vars_frame_prev,
                 Eigen::VectorXf &vars_result,
                 const float weight_p2plane = 1.,
                 const float weight_contour = 0.5,
                 const float weight_reg_zero = 5e-1,
                 const float weight_close_prev = 1e-1,
                 const float weight_nite = 1e-2);


  float
  accumulateModelingConstraints (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                 const pcl::Correspondences &corresps,
                                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                                 af::BlendshapeModel<float>::ConstPtr bs_model,
                                 Eigen::VectorXf &var_bs,
                                 Eigen::MatrixXf &M_lhs, Eigen::VectorXf &M_rhs,
                                 const int num_frames = 1,
                                 const float weight_p2plane = 1., const float weight_contour = 2.);

  static float
  solveModelingProblem (const Eigen::MatrixXf &M_lhs, const Eigen::VectorXf &M_rhs,
                        Eigen::VectorXf &var_bs,
                        const float weight_bs_reg = 1e-1);

  static float
  solveModelingProblemGS (const Eigen::MatrixXf &M_lhs, const Eigen::VectorXf &M_rhs,
                          Eigen::VectorXf &var_bs,
                          const float weight_bs_reg_l1 = 1e-1);


private:
  std::vector<Skeleton::Node*> name_bindings_;
  std::vector<std::vector<std::pair<int, float> > > vertex_to_joints_weights_;
};

}
