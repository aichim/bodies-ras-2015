#include "visualization_utils.h"

#include <QGLViewer/manipulatedCameraFrame.h>
#include <QGLViewer/manipulatedFrame.h>


using namespace std;

af::MeshVisualizationWidget::MeshVisualizationWidget (std::string title,
                                                      QWidget *parent)
  : QGLViewer (parent)
  , enable_texture_ (false)
  , mesh_changed_ (true)
  , tex_id_ (0)
{
  mesh_ = QSharedPointer<af::Mesh<float> > (new af::Mesh<float> ());
  texture_ = new cv::Mat ();

  setWindowTitle (QString (title.c_str ()));
}


void
af::MeshVisualizationWidget::init ()
{
#ifndef __APPLE__
  glewInit ();
#endif

  camera ()->setSceneRadius (5.f);
  qglviewer::Vec pos (0., 0., 0.0);
  camera ()->setPosition (pos);
  camera ()->lookAt (qglviewer::Vec (0., 0., 2.));
  camera ()->frame ()->setSpinningSensitivity (100.f);
  setBackgroundColor (QColor (255, 255, 255));
}


void
af::MeshVisualizationWidget::setTexture (cv::Mat &texture)
{
//  glCheckError ();
  glEnable (GL_TEXTURE_2D);

  cv::flip (texture, *texture_, 0);

  if (tex_id_ == 0)
  {
    glGenTextures (1, &tex_id_);
    PCL_INFO ("Gen texture was called\n");
  }
  glBindTexture (GL_TEXTURE_2D, tex_id_);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
  glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
  glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, texture_->rows, texture_->cols,
                0, GL_BGR, GL_UNSIGNED_BYTE, texture_->ptr ());

//  glCheckError ();
}

void
af::MeshVisualizationWidget::draw ()
{
//  glCheckError ();
  mesh_mutex_.lock ();

  glViewport (0, 0, width (), height ());
  camera ()->loadProjectionMatrix ();
  camera ()->loadModelViewMatrix ();


  glClearColor (0.6, 0.6, 0.6, 0.);
  glClearDepth (1.);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);

  glEnable (GL_LIGHT0);
  glEnable (GL_LIGHTING);

  const GLfloat pos[4] = {0., 0., 100., 1.0};
  glLightfv (GL_LIGHT0, GL_POSITION, pos);
  const GLfloat black[4] = {0., 0., 0., 1.};
  const GLfloat white[4] = {1., 1., 1., 1.};
  const GLfloat gray[4] = {0.7, 0.7, 0.7, 1.};
  glLightfv (GL_LIGHT0, GL_DIFFUSE, black);
  glLightfv (GL_LIGHT0, GL_AMBIENT, gray);
  glEnable (GL_COLOR_MATERIAL);
  glMaterialfv (GL_FRONT_AND_BACK, GL_DIFFUSE, white);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT, white);

  /// TODO set the materials

  glShadeModel (GL_SMOOTH);



  /// Draw the face
  if (enable_texture_ && tex_id_ != 0)
    glBindTexture (GL_TEXTURE_2D, tex_id_);
  else
    glBindTexture (GL_TEXTURE_2D, 0);


  glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
  glLineWidth (0.5);
  glColor3ub (255, 255, 255);

//glCheckError ();
  glBegin (GL_QUADS);
  for (int quad_i = 0; quad_i < mesh_->quad_vertices_.cols (); ++quad_i)
  {
    for (size_t v_i = 0; v_i < 4; ++v_i)
    {
      glNormal3f (mesh_->normals_ (0, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->normals_ (1, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->normals_ (2, mesh_->quad_vertices_ (v_i, quad_i)));
//      glTexCoord2f (mesh_->tex_coords_ (0, mesh_->quad_tex_coords_ (v_i, quad_i)),
//                    mesh_->tex_coords_ (1, mesh_->quad_tex_coords_ (v_i, quad_i)));
      glVertex3f (mesh_->vertices_ (0, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (1, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (2, mesh_->quad_vertices_ (v_i, quad_i)));
    }
  }
  glEnd ();


  glBegin (GL_TRIANGLES);
  for (int tri_i = 0; tri_i < mesh_->tri_vertices_.cols (); ++tri_i)
  {
    for (size_t v_i = 0; v_i < 3; ++v_i)
    {
      glNormal3f (mesh_->normals_ (0, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->normals_ (1, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->normals_ (2, mesh_->tri_vertices_ (v_i, tri_i)));
//      glTexCoord2f (mesh_->tex_coords_ (0, mesh_->tri_tex_coords_ (v_i, tri_i)),
//                    mesh_->tex_coords_ (1, mesh_->tri_tex_coords_ (v_i, tri_i)));
      glVertex3f (mesh_->vertices_ (0, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (1, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (2, mesh_->tri_vertices_ (v_i, tri_i)));
    }
  }
  glEnd ();

/*
  /// Draw tiny spheres for the vertices
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  GLUquadricObj *quadric = gluNewQuadric ();
  gluQuadricDrawStyle (quadric, GLU_FILL);

  for (size_t v_i = 0; v_i < mesh_->vertices_.cols (); ++v_i)
  {
    glPushMatrix();
    glTranslated (mesh_->vertices_ (0, v_i),
                  mesh_->vertices_ (1, v_i),
                  mesh_->vertices_ (2, v_i));
    gluSphere (quadric, .01, 6, 6);
    glPopMatrix();
  }
  gluDeleteQuadric (quadric);
*/

  glBindTexture (GL_TEXTURE_2D, 0);
//  glCheckError ();

  /// Draw the point cloud
  glDisable (GL_TEXTURE_2D);
  glDisable (GL_BLEND);
  glDisable (GL_COLOR_MATERIAL);
  glPointSize (3.);
  glDisable (GL_LIGHT0);
  glDisable (GL_LIGHTING);
  glShadeModel (GL_FLAT);

  glBegin (GL_POINTS);
  if (cloud_)
    for (size_t i = 0; i < cloud_->size (); ++i)
      if (std::isfinite (cloud_->at (i).x))
      {
        glColor3ub (cloud_->at (i).r, cloud_->at (i).g, cloud_->at (i).b);
//        glColor3ub (255, 0, 0);
        glVertex3f (cloud_->at (i).x, cloud_->at (i).y, cloud_->at (i).z);
      }
  glEnd ();


  /// Draw the correspondences
  glLineWidth (0.9);
  glBegin (GL_LINES);
  for (auto c_it = corresps_vec_.begin (); c_it != corresps_vec_.end (); ++c_it)
  {
    glColor3ub (c_it->second (0), c_it->second (1), c_it->second (2));
    for (size_t c_i = 0; c_i < c_it->first.size (); ++c_i)
    {
      glVertex3f (mesh_->vertices_ (0, c_it->first[c_i].index_query),
                  mesh_->vertices_ (1, c_it->first[c_i].index_query),
                  mesh_->vertices_ (2, c_it->first[c_i].index_query));
      glVertex3f (cloud_->at (c_it->first[c_i].index_match).x,
                  cloud_->at (c_it->first[c_i].index_match).y,
                  cloud_->at (c_it->first[c_i].index_match).z);
    }
  }
  glEnd ();



  /// Draw the markers
  glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
  GLUquadricObj *quadric = gluNewQuadric ();
  gluQuadricDrawStyle (quadric, GLU_FILL);
  glColor3ub (255, 255, 50);
  for (size_t i = 0; i < markers_.size (); ++i)
  {
    glPushMatrix();
    glTranslated (markers_[i] (0), markers_[i] (1), markers_[i] (2));
    gluSphere (quadric, .02, 10, 10);
    glPopMatrix();

    qglviewer::Vec screenPos = camera()->projectedCoordinatesOf (qglviewer::Vec (markers_[i] (0), markers_[i] (1), markers_[i] (2)));
    char str[32]; sprintf (str, "%ld", i);
    drawText((int)screenPos[0], (int)screenPos[1], str);
  }

  glColor3ub (255, 255, 255);
  mesh_mutex_.unlock ();
}


void
af::MeshVisualizationWidget::setMesh (const af::Mesh<float> &new_mesh)
{
  mesh_mutex_.lock ();
  mesh_changed_ = true;
  *mesh_ = new_mesh;
  mesh_mutex_.unlock ();

  updateGL ();
}


void
af::MeshVisualizationWidget::getMesh (af::Mesh<float> &mesh)
{
  mesh_mutex_.lock ();
  mesh = *mesh_;
  mesh_mutex_.unlock ();
}


void
af::MeshVisualizationWidget::setPointCloud (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud)
{
  mesh_mutex_.lock ();
  cloud_ = cloud;
  mesh_mutex_.unlock ();

  updateGL ();
}
