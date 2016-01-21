#pragma once

#include "mesh.h"
#include "blendshape_model.h"
#include <pcl/common/common.h>
#include <pcl/console/parse.h>
#include <pcl/console/print.h>
#include <pcl/console/time.h>

#include <pcl/registration/registration.h>
#include <opencv2/opencv.hpp>


namespace af
{
typedef Eigen::Matrix<Eigen::Vector3d, Eigen::Dynamic, Eigen::Dynamic> MatrixXX3d;
typedef Eigen::Matrix<Eigen::Vector4d, Eigen::Dynamic, Eigen::Dynamic> MatrixXX4d;

struct CameraParams
{
  std::string image_path;
  int height, width;
  float fx, fy, cx, cy;
  float skew, k1, k2, k3, p1, p2;
  Eigen::Affine3f pose;

  void
  print () const;

  inline Eigen::Matrix3f
  getK () const
  {
    Eigen::Matrix3f K (Eigen::Matrix3f::Identity ());
    K (0, 0) = fx;
    K (1, 1) = fy;
    K (0, 2) = cx;
    K (1, 2) = cy;
    return (K);
  }

};


void
clamp (double &val, const double min, const double max);


void
getFilesFromDirectory (const std::string path_dir,
                       const std::string extension,
                       std::vector<std::string> &files);

void
getFilesFromDirectory (const std::string path_dir,
                       std::vector<std::string> &files);

void
getDirectoriesFromDirectory (const std::string path_dir,
                             std::vector<std::string> &dirs);

std::string
getBaseName (const std::string &filename);

std::string
getExtension (const std::string &filename);

/**
 * @brief isMatrixHealthy checks if the matrix contains non-finite values (inf or nan)
 * @param matrix
 * @return true if matrix contains only finite values, false otherwise
 */
bool
isMatrixHealthy (const Eigen::MatrixXd &matrix);


bool
readIndicesFile (const std::string &filename,
                 std::vector<size_t> &indices);

bool
writeIndicesFile (const std::string &filename,
                  const std::vector<size_t> &indices);

af::Mesh<float>::Ptr
generateCylinder (const double radius, const double length,
                  const int num_circle_samples = 50, const int num_length_samples = 50);

af::Mesh<float>::Ptr
generateCylinder (Eigen::Vector3f &orig, Eigen::Vector3f &tip,
                  const int num_circle_samples = 50, const int num_length_samples = 50,
                  const double radius = 0.01);

af::Mesh<float>::Ptr
generateSphere (const Eigen::Vector3f &center = Eigen::Vector3f::Zero (),
                const float radius = 0.03,
                const int num_samples = 20);


void
displayBlendshapeWeights (const Eigen::VectorXf &values,
                          cv::Mat &img,
                          const int bar_height = 100,
                          const int bar_width = 10);

void
saveTextForWordCloud (const std::string &filename,
                      const af::BlendshapeModel<float>::ConstPtr bs_model,
                      const Eigen::VectorXf &bs_weights);

bool
saveMatrix (const std::string filename, const Eigen::MatrixXf &matrix);

bool
loadMatrix (const std::string filename, Eigen::MatrixXf &matrix);
}

