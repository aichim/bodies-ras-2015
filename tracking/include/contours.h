#ifndef CONTOURS_H_
#define CONTOURS_H_

#include "common.h"

#ifndef __APPLE__
  #include <GL/glew.h>
#endif

#ifdef __APPLE__
  #include <GLUT/glut.h>
#else
  #include <GL/glut.h>
#endif


namespace af
{
class ContourFBO
{
public:
  ContourFBO ()
  {}

  void
  setMesh (af::Mesh<float>::ConstPtr mesh,
           const std::vector<size_t> &valid_vertex_indices = std::vector<size_t> ())
  {
    mesh_ = mesh;
    valid_vertex_indices_ = valid_vertex_indices;
  }

  void
  setCameraParams (const CameraParams &params);

  void
  draw ();

  void
  getContourPoints (std::vector<int> &contour_indices);

  void
  getRenderImage (cv::Mat &render);

  void
  getDepthMap (cv::Mat &depth_map);

  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

protected:
  af::Mesh<float>::ConstPtr mesh_;
  std::vector<size_t> valid_vertex_indices_;


  Eigen::Matrix4f proj_matrix_;
  Eigen::Matrix4f model_matrix_;

  static GLuint fbo_id_, rb_id_, tex_id_, depth_id_;
  static int win_id_;
  int width_, height_;

  float z_near_, z_far_;

  af::CameraParams camera_params_;
};
}


#endif // CONTOURS_H_
