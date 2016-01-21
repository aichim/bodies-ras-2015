#include "contours.h"

void
af::ContourFBO::setCameraParams (const CameraParams &params)
{
  proj_matrix_ = Eigen::Matrix4f::Identity ();

  z_near_ = 0.1;
  z_far_ = 50.;

  proj_matrix_ (0, 0) = 2. * params.fx / static_cast<float> (params.width);
  proj_matrix_ (0, 2) = (static_cast<float> (params.width) - 2. * params.cx) / static_cast<float> (params.width);
  proj_matrix_ (1, 1) = 2. * params.fy / static_cast<float> (params.height);
  proj_matrix_ (1, 2) = (-static_cast<float> (params.height) + 2. * params.cy) / static_cast<float> (params.height);
  proj_matrix_ (2, 2) = (-z_far_ - z_near_) / (z_far_ - z_near_);
  proj_matrix_ (2, 3) = -2. * z_far_ * z_near_ / (z_far_ - z_near_);
  proj_matrix_ (3, 2) = -1.;
  proj_matrix_ (3, 3) = 0;

  model_matrix_ = params.pose.inverse ().matrix ().cast<float> ();
  model_matrix_.row (1) *= -1.;
  model_matrix_.row (2) *= -1.;

  width_ = params.width;
  height_ = params.height;

  camera_params_ = params;
}


/*
void
af::ContourFBO::draw ()
{
  // glCheckError ();

  /// Do the actual drawing
  glDisable (GL_TEXTURE_2D);
  glDisable (GL_BLEND);
  glDisable (GL_COLOR_MATERIAL);
  glClearColor (0., 0., 0., 0.);
  glClearDepth (1.);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  glDisable (GL_LIGHT0);
  glDisable (GL_LIGHTING);

  // glCheckError ();

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  glMultMatrixf (proj_matrix_.data ());

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  glMultMatrixf (model_matrix_.data ());

  // glCheckError ();
  glShadeModel (GL_FLAT);

  // glCheckError ();

  // if (!wireframe_)
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
  // else
    // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

  glBegin (GL_QUADS);
  for (int quad_i = 0; quad_i < mesh_->quad_vertices_.cols (); ++quad_i)
  {
    for (size_t v_i = 0; v_i < 4; ++v_i)
    {
      glVertex3f (mesh_->vertices_ (0, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (1, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (2, mesh_->quad_vertices_ (v_i, quad_i)));
    }
  }
  glEnd ();

  // glCheckError ();

  glBegin (GL_TRIANGLES);
  for (int tri_i = 0; tri_i < mesh_->tri_vertices_.cols (); ++tri_i)
  {
    unsigned char r, g, b;
    size_t temp = tri_i + 1; /// to avoid confusion with the black background
    r = temp / (256 * 256);
    g = (temp % (256 * 256)) / 256;
    b = temp % 256;

    glColor3ub (r, g, b);
    for (size_t v_i = 0; v_i < 3; ++v_i)
    {
      glVertex3f (mesh_->vertices_ (0, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (1, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (2, mesh_->tri_vertices_ (v_i, tri_i)));
    }
  }
  glEnd ();

  // glCheckError ();

  glColor3ub (255, 255, 255);
}*/


void
af::ContourFBO::draw ()
{
  // glCheckError ();

  /// Do the actual drawing
  glDisable (GL_TEXTURE_2D);
//  glEnable (GL_BLEND);
  glEnable (GL_COLOR_MATERIAL);
  glClearColor (0., 0., 0., 0.);
  glClearDepth (1.);
  glClear (GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glEnable(GL_DEPTH_TEST);
  const GLfloat pos[4] = {0., 20., -20., 1.0};
  glLightfv (GL_LIGHT0, GL_POSITION, pos);
  const GLfloat black[4] = {0., 0., 0., 1.};
  const GLfloat white[4] = {1., 1., 1., 1.};
  const GLfloat darkgrey[4] = {0.1, 0.1, 0.1, 1.};
  const GLfloat grey[4] = {0.5, 0.5, 0.5, 1.};


  glEnable (GL_LIGHT0);
  glEnable (GL_LIGHTING);

  glDisable (GL_TEXTURE_2D);
  const GLfloat material_color[4] = {0.55, 0.75, 0.97, 1.};
  glColor3fv (material_color);
  glLightfv (GL_LIGHT0, GL_SPECULAR, black);
  glLightfv (GL_LIGHT0, GL_AMBIENT, darkgrey);
  glLightfv (GL_LIGHT0, GL_DIFFUSE, grey);



  // glCheckError ();

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();
  glMultMatrixf (proj_matrix_.data ());

  glMatrixMode (GL_MODELVIEW);
  glLoadIdentity ();
  glMultMatrixf (model_matrix_.data ());

  // glCheckError ();
  glShadeModel (GL_SMOOTH);

  // glCheckError ();

  // if (!wireframe_)
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
  // else
    // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

  glBegin (GL_QUADS);
  for (int quad_i = 0; quad_i < mesh_->quad_vertices_.cols (); ++quad_i)
  {
    for (size_t v_i = 0; v_i < 4; ++v_i)
    {
      glNormal3f (mesh_->normals_ (0, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->normals_ (1, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->normals_ (2, mesh_->quad_vertices_ (v_i, quad_i)));
      glVertex3f (mesh_->vertices_ (0, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (1, mesh_->quad_vertices_ (v_i, quad_i)),
                  mesh_->vertices_ (2, mesh_->quad_vertices_ (v_i, quad_i)));
    }
  }
  glEnd ();

  // glCheckError ();

  glBegin (GL_TRIANGLES);
  for (int tri_i = 0; tri_i < mesh_->tri_vertices_.cols (); ++tri_i)
  {
    for (size_t v_i = 0; v_i < 3; ++v_i)
    {
      glNormal3f (mesh_->normals_ (0, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->normals_ (1, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->normals_ (2, mesh_->tri_vertices_ (v_i, tri_i)));
      glVertex3f (mesh_->vertices_ (0, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (1, mesh_->tri_vertices_ (v_i, tri_i)),
                  mesh_->vertices_ (2, mesh_->tri_vertices_ (v_i, tri_i)));
    }
  }
  glEnd ();

  // glCheckError ();
}


GLuint af::ContourFBO::fbo_id_ = 0;
GLuint af::ContourFBO::rb_id_ = 0;
GLuint af::ContourFBO::tex_id_ = 0;
GLuint af::ContourFBO:: depth_id_ = 0;
int af::ContourFBO::win_id_ = 0;


void
display ()
{}


void
af::ContourFBO::getRenderImage (cv::Mat &render)
{
  /// Get current OpenGL state
  int fbo_id_prev, rb_id_prev;
  glGetIntegerv (GL_FRAMEBUFFER_BINDING, &fbo_id_prev);
  glGetIntegerv (GL_RENDERBUFFER_BINDING, &rb_id_prev);


  /// Check if GLUT was initialized already (by ourselves or by VTK)
  /// Don't do it a second time
  if (glutGetWindow () == 0)
  {
    PCL_ERROR ("Creating glut window\n");
    int argc = 0;
    char** argv = nullptr;
    glutInit (&argc, argv);
    glutInitWindowSize (320, 320);
    win_id_ = glutCreateWindow ("glut dummy window render_with_texture");
#ifndef __APPLE__
    glewInit ();
#endif
    glutDisplayFunc (display);
  }

  win_id_ = glutGetWindow ();
  glutSetWindow (win_id_);
  // glCheckError ();

  if (rb_id_ == 0)
  {
    glGenFramebuffers (1, &fbo_id_);
    glGenRenderbuffers (1, &rb_id_);
    glGenTextures (1, &tex_id_);
    glBindFramebuffer (GL_FRAMEBUFFER, fbo_id_);

    glBindTexture (GL_TEXTURE_2D, tex_id_);
    glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, width_, height_, 0, GL_RGB, GL_UNSIGNED_BYTE, 0);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

    glFramebufferTexture2D (GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, tex_id_, 0);

    glBindRenderbuffer (GL_RENDERBUFFER, rb_id_);
    glRenderbufferStorage (GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width_, height_);

    glFramebufferRenderbuffer (GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_id_);

    GLenum status = glCheckFramebufferStatus (GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE)
      PCL_ERROR ("Error with our framebuffer.\n");
  }

  // glCheckError ();

  glViewport (0, 0, width_, height_);

  // glCheckError ();

  draw ();

  // glCheckError ();

  /// Get the data from the texture into a cv::Mat
  glReadBuffer (GL_COLOR_ATTACHMENT0);

  // glCheckError ();

  render = cv::Mat::zeros (height_, width_, CV_8UC3);
  glPixelStorei (GL_PACK_ALIGNMENT, (render.step & 3) ? 1 : 4);
  glPixelStorei (GL_PACK_ROW_LENGTH, render.step / render.elemSize ());
  glReadPixels (0, 0, width_, height_, GL_BGR, GL_UNSIGNED_BYTE, render.data);
  cv::flip (render, render, 0);

  // glCheckError ();

  glDeleteTextures (1, &tex_id_);
  glDeleteRenderbuffers (1, &rb_id_);
  glBindFramebuffer (GL_FRAMEBUFFER, 0);
  glDeleteFramebuffers (1, &fbo_id_);
  tex_id_ = rb_id_ = fbo_id_ = 0;
  // glCheckError ();


  /// Put back the old OpenGL state
  glBindFramebuffer (GL_FRAMEBUFFER, fbo_id_prev);
  glBindRenderbuffer (GL_RENDERBUFFER, rb_id_prev);
}



void
af::ContourFBO::getDepthMap (cv::Mat &depth_map)
{
  /// Get current OpenGL state
  int fbo_id_prev, rb_id_prev;
  glGetIntegerv (GL_FRAMEBUFFER_BINDING, &fbo_id_prev);
  glGetIntegerv (GL_RENDERBUFFER_BINDING, &rb_id_prev);

  /// Check if GLUT was initialized already (by ourselves or by VTK)
  /// Don't do it a second time
  if (glutGetWindow () == 0)
  {
    printf ("Creating glut window\n");
    int argc = 0;
    char** argv = nullptr;
    glutInit (&argc, argv);
    glutInitWindowSize (320, 320);
    win_id_ = glutCreateWindow ("glut dummy window render_with_texture");
#ifndef __APPLE__
    glewInit ();
#endif
    glutDisplayFunc (display);
  }

  win_id_ = glutGetWindow ();
  glutSetWindow (win_id_);
//  glCheckError ();

  if (rb_id_ == 0)
  {
    glGenFramebuffers (1, &fbo_id_);
    glGenRenderbuffers (1, &rb_id_);
    glGenTextures (1, &depth_id_);
    glBindFramebuffer (GL_FRAMEBUFFER, fbo_id_);

    glBindTexture (GL_TEXTURE_2D, depth_id_);
    glTexImage2D (GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT24, width_, height_, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

    glFramebufferTexture2D (GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, depth_id_, 0);

    glBindRenderbuffer (GL_RENDERBUFFER, rb_id_);
    glRenderbufferStorage (GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width_, height_);
    glFramebufferRenderbuffer (GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_id_);

    glDrawBuffer(GL_NONE);
    glReadBuffer(GL_NONE);

    GLenum status = glCheckFramebufferStatus (GL_FRAMEBUFFER);
    if (status != GL_FRAMEBUFFER_COMPLETE)
      printf ("Error with our framebuffer.\n");
  }

//  glCheckError ();

  glViewport (0, 0, width_, height_);

//  glCheckError ();

  draw ();

//  glCheckError ();

  /// Get the data from the texture into a cv::Mat
  glReadBuffer (GL_DEPTH_ATTACHMENT);

//  glCheckError ();

  depth_map = cv::Mat::zeros (height_, width_, CV_32FC1);
  glPixelStorei (GL_PACK_ALIGNMENT, (depth_map.step & 3) ? 1 : 4);
  glPixelStorei (GL_PACK_ROW_LENGTH, depth_map.step / depth_map.elemSize ());
  glReadPixels (0, 0, width_, height_, GL_DEPTH_COMPONENT, GL_FLOAT, depth_map.data);
  cv::flip (depth_map, depth_map, 0);

//  glCheckError ();

  glDeleteTextures (1, &depth_id_);
  glDeleteRenderbuffers (1, &rb_id_);
  glBindFramebuffer (GL_FRAMEBUFFER, 0);
  glDeleteFramebuffers (1, &fbo_id_);
  tex_id_ = rb_id_ = fbo_id_ = 0;
//  glCheckError ();


  /// Put back the old OpenGL state
  glBindFramebuffer (GL_FRAMEBUFFER, fbo_id_prev);
  glBindRenderbuffer (GL_RENDERBUFFER, rb_id_prev);
}


void
af::ContourFBO::getContourPoints (std::vector<int> &contour_indices)
{
  cv::Mat depthmap;
  getDepthMap (depthmap);
//  cv::imshow ("depth", depthmap);

  /// Convert the depth map to meters
  depthmap = 2.0 * z_near_ * z_far_ / (z_far_ + z_near_ - (2.0 * depthmap - 1.0) * (z_far_ - z_near_));
//  cv::imshow ("depth meters", depthmap);


  cv::Mat contour_depth = cv::Mat::zeros (height_, width_, CV_8UC1);
  float d;
  for (int y = 0; y < depthmap.rows; ++y)
    for (int x = 0; x < depthmap.cols; ++x)
    {
      d = depthmap.at<float> (cv::Point (x, y));
      if (d > z_far_ * 0.9)
        continue;

      bool stop = false;
      for (int dx = -1; dx <= 1 && !stop; ++dx)
        for (int dy = -1; dy <= 1 && !stop; ++dy)
          if (fabs (d - depthmap.at<float> (cv::Point (x + dx, y + dy))) > 0.1)
          {
            contour_depth.at<unsigned char> (cv::Point (x, y)) = 255;
            stop = true;
          }
    }
  cv::dilate (contour_depth, contour_depth, cv::Mat ());

  Eigen::Matrix3f K = camera_params_.getK ().cast<float> ();
  cv::Mat image_sparse = cv::Mat::zeros (height_, width_, CV_8UC3);
  for (size_t i = 0; i < valid_vertex_indices_.size (); ++i)
  {
    Eigen::Vector3f p_proj = K * (camera_params_.pose.inverse ().cast<float> () * mesh_->vertices_.col (valid_vertex_indices_[i]));
    p_proj /= p_proj (2);

    if (contour_depth.at<unsigned char> (cv::Point (p_proj (0), p_proj (1))) != 0)
    {
      cv::circle (image_sparse, cv::Point (p_proj (0), p_proj (1)), 2, cv::Scalar (0, 0, 255));
      contour_indices.push_back (valid_vertex_indices_[i]);
    }
  }

//  cv::imshow ("contour_depth", contour_depth);
//  cv::imshow ("image contour sparse from depth", image_sparse);
//  cv::waitKey ();
}

