#ifndef VISUALIZATION_UTILS_H_
#define VISUALIZATION_UTILS_H_

#ifndef __APPLE__
  #include <GL/glew.h>
#endif

#include "common.h"

#include <QGLViewer/qglviewer.h>
#include <QtOpenGL/QGLWidget>
#include <QMutex>

namespace af
{
  struct CameraParams;

class MeshVisualizationWidget : public QGLViewer
{
public:
  MeshVisualizationWidget (std::string title = "", QWidget *parent = NULL);

  void
  setMesh (const af::Mesh<float> &new_mesh);

  void
  getMesh (af::Mesh<float> &mesh);

  void
  setPointCloud (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud);

  void
  addCorrespondences (const pcl::Correspondences &corresps, const Eigen::Vector3i &color)
  {
    mesh_mutex_.lock ();
    corresps_vec_.push_back (std::make_pair (corresps, color));
    mesh_mutex_.unlock ();
  }

  void
  clearCorrespondences ()
  {
    mesh_mutex_.lock ();
    corresps_vec_.clear ();
    mesh_mutex_.unlock ();
  }


  void
  setMarkers (const std::vector<Eigen::Vector3f> &markers)
  {
    mesh_mutex_.lock ();
    markers_ = markers;
    mesh_mutex_.unlock ();

    updateGL ();
  }


  void
  setTexture (cv::Mat &texture);

  bool
  enableTexture ()
  {
    mesh_mutex_.lock ();
    enable_texture_ = !enable_texture_;
    mesh_mutex_.unlock ();

    updateGL ();
    return (enable_texture_);
  }

  void
  setEnableTexture (bool enable_texture)
  {
    mesh_mutex_.lock ();
    enable_texture_ = enable_texture;
    mesh_mutex_.unlock ();

    updateGL ();
  }

protected:
  virtual void
  draw ();

  virtual void
  init ();


  QSharedPointer<af::Mesh<float> > mesh_;
  pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud_;

  std::vector<std::pair<pcl::Correspondences, Eigen::Vector3i> > corresps_vec_;
  std::vector<Eigen::Vector3f> markers_;


  QMutex mesh_mutex_;
  bool enable_texture_;
  cv::Mat *texture_;
  GLuint tex_id_;

  bool mesh_changed_;
};
}

#endif // VISUALIZATION_UTILS_H_
