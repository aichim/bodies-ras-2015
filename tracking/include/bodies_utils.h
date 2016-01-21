#pragma once

#include "common.h"
#include "Fbx_IO.h"


namespace af
{

struct SkeletonState
{
  enum TrackingState {NOT_TRACKED, INFERRED, TRACKED};
  std::vector<Eigen::Vector3f> joint_pos;
  std::vector<Eigen::Vector2f> joint_uv;
  std::vector<TrackingState> tracking_state;
  std::vector<float> quality_flag;
  std::vector<Eigen::Vector2i> link;
  std::vector<Eigen::Quaternionf> quat;

  SkeletonState (std::string &filename);

  bool
  load (std::string &filename);

  int
  numJoints () const { return (joint_pos.size ()); }
};


bool
loadBodyTrackingRGBDDataset (std::string &folder,
                             std::vector<SkeletonState*> &skeleton_frames,
                             std::vector<pcl::PointCloud<pcl::PointXYZRGBA>::Ptr> &point_clouds,
                             af::CameraParams &camera_params_rgb, af::CameraParams &camera_params_depth);

struct PosePCA
{
  PosePCA () {}
  PosePCA (const Eigen::VectorXf &_mean,
           const Eigen::VectorXf &_std_devs,
           const Eigen::MatrixXf &_modes)
    : mean (_mean), std_devs (_std_devs), modes (_modes)
  {}

  Eigen::VectorXf mean;
  Eigen::VectorXf std_devs;
  Eigen::MatrixXf modes;

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

bool
loadPCA (std::string path,
         PosePCA &pca);

bool
savePCA (std::string path,
         const PosePCA &pca);


void
clampAngles (float &angles);


class SkeletonMesh;
bool
saveTrackingAndModelingResults (const std::string &filename,
                                const af::SkeletonMesh &skeleton_mesh,
                                const std::vector<Eigen::VectorXf> &tracking_results);

bool
loadTrackingAndModelingResults (const std::string &filename,
                                af::SkeletonMesh &skeleton_mesh,
                                std::vector<Eigen::VectorXf> &tracking_results);

}
