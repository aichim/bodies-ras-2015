#include "bodies_utils.h"
#include <pcl/visualization/pcl_visualizer.h>
#include "skeleton_mesh.h"
#include <pcl/filters/morphological_filter.h>


//#define DEBUG
//#define BIWI


void
erodePointcloud (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud_in,
                 pcl::PointCloud<pcl::PointXYZRGBA> &cloud_out,
                 const int window_size = 1)
{
  cloud_out = *cloud_in;
  pcl::PointXYZRGBA p_nan;
  p_nan.x = p_nan.y = p_nan.z = std::numeric_limits<float>::quiet_NaN ();
  for (size_t v = window_size; v < cloud_in->height - window_size; ++v)
    for (size_t u = window_size; u < cloud_in->width - window_size; ++u)
    {
      int count_valid = 0;
      for (int du = -window_size; du <= window_size; ++du)
        for (int dv = -window_size; dv <= window_size; ++dv)
          if (pcl_isfinite ((*cloud_in)[(v + dv) * 640 + u + du].z))
            count_valid ++;

      if (count_valid != (2 * window_size + 1) * (2 * window_size + 1))
        cloud_out[v * 640 + u] = p_nan;
    }

}

bool
af::loadBodyTrackingRGBDDataset (std::string &folder,
                                 std::vector<SkeletonState*> &skeleton_frames,
                                 std::vector<pcl::PointCloud<pcl::PointXYZRGBA>::Ptr> &clouds,
                                 af::CameraParams &camera_rgb, af::CameraParams &camera_depth)
{
  std::vector<std::string> files;
  af::getFilesFromDirectory (folder, files);

#ifdef BIWI
  int num_frames = files.size () / 5;
#else
  int num_frames = (files.size ()) / 4;
#endif
  std::vector<cv::Mat*> depth_imgs, color_imgs, user_imgs;

  for (size_t f_i = 0; f_i < files.size (); ++f_i)
  {
    std::string base = af::getBaseName (files[f_i]);
//    std::cout << "reading file: " << base << std::endl;
    if (base.find ("_depth") != std::string::npos)
    {
      /// Depth map
      cv::Mat *depth_img = new cv::Mat ();
      *depth_img = cv::imread (files[f_i], -1);
      depth_imgs.push_back (depth_img);
    }
    else if (base.find ("_rgb") != std::string::npos)
    {
      /// Color image
      cv::Mat *color_img = new cv::Mat ();
      *color_img = cv::imread (files[f_i]);
      color_imgs.push_back (color_img);
    }
    else if (base.find ("_skel") != std::string::npos)
    {
      /// Skeleton information file
      skeleton_frames.push_back (new SkeletonState (files[f_i]));
    }
    else if (base.find ("_userMap") != std::string::npos)
    {
      /// User map image
      cv::Mat *user_img = new cv::Mat ();
      *user_img = cv::imread (files[f_i]);
      cv::cvtColor (*user_img, *user_img, cv::COLOR_BGR2GRAY);
//      *user_img *= 100;
      user_imgs.push_back (user_img);
    }
    else
    {
//      std::cerr << "skipping file: " << files[f_i] << std::endl;
    }
  }

  PCL_ERROR ("--> Read in data: depth maps %zu, color images %zu, user images %zu, skeletons %zu.\n",
             depth_imgs.size (), color_imgs.size (), user_imgs.size (), skeleton_frames.size ());

  /// Hardcode the camera parameters
  camera_rgb.fx = 538.9; //525.;
  camera_rgb.fy = 539.2; //525.;
  camera_rgb.cx = 321.2; //319.5;
  camera_rgb.cy = 244.3;  //239.5;
#ifdef BIWI
  camera_rgb.fx *= 2.;
  camera_rgb.fy *= 2.;
  camera_rgb.cx *= 2.;
  camera_rgb.cy *= 2.;
#endif
  camera_rgb.width = color_imgs.front ()->cols;
  camera_rgb.height = color_imgs.front ()->rows;
  Eigen::Matrix4f pose (Eigen::Matrix4f::Identity ());
  pose.block<3, 1> (0, 3) = Eigen::Vector3f (0.025, 0., 0.);
  camera_rgb.pose = Eigen::Affine3f (pose);

  camera_depth.fx = 563.2; //575.8;
  camera_depth.fy = 562.2; //575.8;
  camera_depth.cx = 320.49; //319.5;
  camera_depth.cy = 244.7; //239.5;
  camera_depth.width = depth_imgs.front ()->cols;
  camera_depth.height = depth_imgs.front ()->rows;
  camera_depth.pose = Eigen::Affine3f::Identity ();

  /// HACK - seems to work for hardware-registered OpenNI frames
  camera_rgb = camera_depth;



  pcl::PointXYZRGBA p_nan;
  p_nan.x = p_nan.y = p_nan.z = std::numeric_limits<float>::quiet_NaN ();

  /// Compose the point clouds
  clouds.resize (num_frames);
  for (size_t frame_i = 0; frame_i < num_frames; ++frame_i)
  {
    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud (new pcl::PointCloud<pcl::PointXYZRGBA> ());
    cloud->width = depth_imgs[frame_i]->cols;
    cloud->height = depth_imgs[frame_i]->rows;
    cloud->points.resize (cloud->width * cloud->height, p_nan);
    cloud->is_dense = false;
    for (size_t v = 0; v < depth_imgs[frame_i]->rows; ++v)
      for (size_t u = 0; u < depth_imgs[frame_i]->cols; ++u)
      {
        /// TODO unsigned short for the Matteo datasets
        if (user_imgs[frame_i]->at<unsigned char> (cv::Point (u, v)) == 0)
        {
//          cloud->push_back (p_nan);
          continue;
        }

        double z = 1e-3 * static_cast<double> (depth_imgs[frame_i]->at<unsigned short> (cv::Point (u, v)));

        pcl::PointXYZRGBA p;
        p.getVector3fMap () = Eigen::Vector3f (z / camera_depth.fx * (static_cast<float> (u) - camera_depth.cx),
                                               z / camera_depth.fy * (static_cast<float> (camera_depth.height - 1) - static_cast<float> (v) - camera_depth.cy),
//                                               z / camera_depth.fy * (static_cast<double> (v) - camera_depth.cy),
                                               z);

        if (pcl_isfinite (p.z) && fabs (p.z) > 1e-3)
        {
          Eigen::Vector3f p_proj = camera_depth.getK () * (p.getVector3fMap ());
          p_proj /= p_proj (2);

          /// TODO bilinear interpolation here
          cv::Vec3b color = color_imgs[frame_i]->at<cv::Vec3b> (cv::Point (p_proj (0), static_cast<float> (camera_rgb.height - 1) - p_proj (1)));
//          cv::Vec3b color = color_imgs[frame_i]->at<cv::Vec3b> (cv::Point (p_proj (0), p_proj (1)));
          p.r = color [2];
          p.g = color [1];
          p.b = color [0];
          p.a = 255;
        }

        int index = (camera_depth.height - 1 - v) * cloud->width + u;
        (*cloud)[index] = p;
//        cloud->push_back (p);
//        PCL_ERROR ("DEBUG shit: index (%zu %zu height %d) -> %d, clouds size by push_back %zu\n",
//                   u, v, cloud->height, index, cloud->size ());
      }
    clouds[frame_i] = cloud;

    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_eroded (new pcl::PointCloud<pcl::PointXYZRGBA> ());
    for (size_t iter_i = 0; iter_i < 2; ++iter_i)
    {
      erodePointcloud (clouds[frame_i], *cloud_eroded);
      *clouds[frame_i] = *cloud_eroded;
    }
  }

#ifdef DEBUG
  /// Visualize the data
  pcl::visualization::PCLVisualizer vis, vis_erode;
  pcl::visualization::PointCloudColorHandlerRGBAField<pcl::PointXYZRGBA> rgba (clouds[0]);
  vis.addPointCloud<pcl::PointXYZRGBA> (clouds[0], rgba, "cloud");
  for (size_t frame_i = 0; frame_i < num_frames; ++frame_i)
  {
    PCL_ERROR ("Showing frame %zu / %d\n", frame_i, num_frames);
    cv::imshow ("depth", *depth_imgs[frame_i]);
    cv::imshow ("color", *color_imgs[frame_i]);
    cv::imshow ("user map", *user_imgs[frame_i]);

    cv::Mat img_skeleton = cv::Mat::zeros (color_imgs[frame_i]->rows, color_imgs[frame_i]->cols, CV_8UC3);
    for (size_t j_i = 0; j_i < skeleton_frames[frame_i]->numJoints (); ++j_i)
    {
      cv::circle (img_skeleton, cv::Point (skeleton_frames[frame_i]->joint_uv[j_i] (0),
                                           skeleton_frames[frame_i]->joint_uv[j_i] (1)),
                                           2, cv::Scalar (0, 0, 255), -1);
      char str_num[8];
      sprintf (str_num, "%zu", j_i);
      cv::putText (img_skeleton, str_num, cv::Point (skeleton_frames[frame_i]->joint_uv[j_i] (0),
                                                     skeleton_frames[frame_i]->joint_uv[j_i] (1)),
                                                     cv::FONT_HERSHEY_SIMPLEX, 1., cv::Scalar (0, 0, 255));


      Eigen::Vector3f p_proj = camera_rgb.getK () * skeleton_frames[frame_i]->joint_pos[j_i];
      p_proj /= p_proj (2);
      cv::circle (img_skeleton, cv::Point (p_proj (0), p_proj (1)),
                  2, cv::Scalar (255, 0, 0), -1);
    }
    cv::imshow ("skeleton uv", img_skeleton);

    pcl::visualization::PointCloudColorHandlerRGBAField<pcl::PointXYZRGBA> rgba (clouds[frame_i]);
    vis.updatePointCloud<pcl::PointXYZRGBA> (clouds[frame_i], rgba, "cloud");


    pcl::visualization::PointCloudColorHandlerRGBAField<pcl::PointXYZRGBA> rgba_eroded (cloud_eroded);
    if (!vis_erode.updatePointCloud<pcl::PointXYZRGBA> (cloud_eroded, rgba_eroded, "cloud_eroded"))
        vis_erode.addPointCloud<pcl::PointXYZRGBA> (cloud_eroded, rgba_eroded, "cloud_eroded");


    vis.spin ();
    cv::waitKey (1);
  }
#endif


  return (true);
}


af::SkeletonState::SkeletonState (std::string &filename)
{
  load (filename);
}

/*
bool
af::SkeletonState::load (std::string &filename)
{
  std::ifstream file (filename);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening skeleton file %s for reading.\n", filename.c_str ());
    return (false);
  }

  std::string str;
  while (file.good ())
  {
    /// Joint 3d position
    Eigen::Vector3d pos;
    getline (file, str, ',');
    std::stringstream ss (str);
    ss >> pos[0];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> pos[1];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> pos[2];
    this->joint_pos.push_back (pos);

    /// Confidence
    getline (file, str, '\n');
    ss.clear ();
    ss.str (str);
    double conf;
    ss >> conf;
  }

  return (true);
}
*/

bool
af::SkeletonState::load (std::string &filename)
{
  std::ifstream file (filename);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening skeleton file %s for reading.\n", filename.c_str ());
    return (false);
  }

  std::string str;
  while (file.good ())
  {
    /// Person id
    getline (file, str, ',');
    if (!file.good ())
      continue;
    std::stringstream ss (str);
    double person_id;
    ss >> person_id;

//    if (person_id < 0.5)
//      continue;

    /// Joint 3d position
    Eigen::Vector3f pos;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> pos[0];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> pos[1];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> pos[2];
    /// TODO uncomment this for my own datasets
    pos /= 1000.;
    this->joint_pos.push_back (pos);

    /// Joint 2d projection
    Eigen::Vector2f uv;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> uv[0];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> uv[1];
    /// HACK for the BiWi datasets
#ifdef BIWI
    uv *= 2.;
#endif

    this->joint_uv.push_back (uv);

    /// Joint tracking state
    int ts;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> ts;
    this->tracking_state.push_back (static_cast<SkeletonState::TrackingState> (ts));

    /// Quality flag
    float quality;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> quality;
    this->quality_flag.push_back (quality);

    /// Orientation indices
    Eigen::Vector2i link;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> link[0];
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> link[1];
    this->link.push_back (link);

    /// Quaternion
    Eigen::Quaternionf quat;
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> quat.x ();
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> quat.y ();
    getline (file, str, ',');
    ss.clear ();
    ss.str (str);
    ss >> quat.z ();
    getline (file, str, '\n');
    ss.clear ();
    ss.str (str);
    ss >> quat.w ();
    this->quat.push_back (quat);

/*
    if (this->tracking_state.back () == SkeletonState::NOT_TRACKED)
    {
      this->joint_pos.pop_back ();
      this->joint_uv.pop_back ();
      this->quality_flag.pop_back ();
      this->link.pop_back ();
      this->quat.pop_back ();
    }*/
  }

#ifndef BIWI

#endif

//  printf ("frame has %d joint positions\n", joint_pos.size ());
  return (true);
}


bool
af::loadPCA (std::string path,
             PosePCA &pca)
{
  std::ifstream file (path.c_str (), std::ios_base::binary | std::ios_base::in);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening PCA file %s\n", path.c_str ());
    return (false);
  }

  int32_t file_id;
  file.read (reinterpret_cast<char*> (&file_id), sizeof (file_id));
  PCL_INFO ("PCA model file id: %ld\n", file_id);
  if (file_id != 1121)
  {
    PCL_ERROR ("Old PCA model found\n");
    return (false);
  }

  int32_t num_dims, num_modes;
  file.read (reinterpret_cast<char*> (&num_dims), sizeof (num_dims));
  file.read (reinterpret_cast<char*> (&num_modes), sizeof (num_modes));

  pca.mean.resize (num_dims, 1);
  pca.modes.resize (num_dims, num_modes);
  pca.std_devs.resize (num_modes, 1);

  file.read (reinterpret_cast<char*> (&pca.mean (0, 0)), sizeof (float) * num_dims);
  file.read (reinterpret_cast<char*> (&pca.modes (0, 0)), sizeof (float) * num_dims * num_modes);
  file.read (reinterpret_cast<char*> (&pca.std_devs (0, 0)), sizeof (float) * num_modes);

  file.close ();

  return (true);
}


bool
af::savePCA (std::string path,
             const PosePCA &pca)
{
  std::ofstream file (path.c_str (), std::ios_base::binary | std::ios_base::out);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening PCA file %s\n", path.c_str ());
    return (false);
  }

  const int32_t file_id = 1121;
  file.write (reinterpret_cast<const char*> (&file_id), sizeof (file_id));

  const int32_t num_dims = pca.mean.rows (),
                num_modes = pca.modes.cols ();
  printf ("saving pca pose model with %d dims and %d modes\n", num_dims, num_modes);
  file.write (reinterpret_cast<const char*> (&num_dims), sizeof (num_dims));
  file.write (reinterpret_cast<const char*> (&num_modes), sizeof (num_modes));

  file.write (reinterpret_cast<const char*> (&pca.mean (0, 0)), sizeof (float) * num_dims);
  file.write (reinterpret_cast<const char*> (&pca.modes (0, 0)), sizeof (float) * num_dims * num_modes);
  file.write (reinterpret_cast<const char*> (&pca.std_devs (0, 0)), sizeof (float) * num_modes);

  file.close ();

  return (true);
}


void
af::clampAngles (float &angles)
{
//  for (size_t i = 0; i < angles.rows (); ++i)
//  {
    while (angles < -M_PI)
      angles  += 2. * M_PI;
    while (angles  > M_PI)
      angles  -= 2. * M_PI;
//  }
}


bool
af::saveTrackingAndModelingResults (const std::string &filename,
                                    const af::SkeletonMesh &skeleton_mesh,
                                    const std::vector<Eigen::VectorXf> &tracking_results)
{
  std::ofstream file (filename);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error: could not open results file %s for writing.\n", filename.c_str ());
    return (false);
  }

  /// Save skeleton data
  file << skeleton_mesh.skeleton_rest_.nodes.size () << std::endl;
  for (size_t n_i = 0; n_i < skeleton_mesh.skeleton_rest_.nodes.size (); ++n_i)
  {
    for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
        file << skeleton_mesh.skeleton_rest_.nodes[n_i]->local_transformation (i, j) << " ";
    file << std::endl;
  }

  /// Save rest mesh data
  file << skeleton_mesh.mesh_rest_->vertices_.cols () << std::endl;
  for (size_t v_i = 0; v_i < skeleton_mesh.mesh_rest_->vertices_.cols (); ++v_i)
    file << skeleton_mesh.mesh_rest_->vertices_.col (v_i).transpose () << std::endl;

  /// Save the frame tracking data
  file << tracking_results.size () << std::endl;
  for (size_t f_i = 0; f_i < tracking_results.size (); ++f_i)
    file << tracking_results[f_i].transpose () << std::endl;

  file.close ();

  return (true);
}


bool
af::loadTrackingAndModelingResults (const std::string &filename,
                                    af::SkeletonMesh &skeleton_mesh,
                                    std::vector<Eigen::VectorXf> &tracking_results)
{
  std::ifstream file (filename);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error: could not open results file %s for reading.\n", filename.c_str ());
    return (false);
  }

  /// Load the skeleton data - assume the skeleton has been loaded - we only modify the local transformations
  int num_nodes;
  file >> num_nodes;
  if (num_nodes != skeleton_mesh.skeleton_rest_.nodes.size ())
  {
    PCL_ERROR ("Error: number of nodes differs: %d vs %zu.\n", num_nodes, skeleton_mesh.skeleton_rest_.nodes.size ());
    return (false);
  }
  for (size_t n_i = 0; n_i < num_nodes; ++n_i)
    for (size_t i = 0; i < 4; ++i)
      for (size_t j = 0; j < 4; ++j)
        file >> skeleton_mesh.skeleton_rest_.nodes[n_i]->local_transformation (i, j);

  /// Load the mesh data - assumes the topology of the mesh is the same
  int num_vertices;
  file >> num_vertices;
  if (num_vertices != skeleton_mesh.mesh_rest_->vertices_.cols ())
  {
    PCL_ERROR ("Error: number of vertices differs: %d vs %d.\n", num_vertices, skeleton_mesh.mesh_rest_->vertices_.cols ());
    return (false);
  }
  for (size_t v_i = 0; v_i < num_vertices; ++v_i)
    for (size_t i = 0; i < 3; ++i)
      file >> skeleton_mesh.mesh_rest_->vertices_ (i, v_i);

  /// Load the frame tracking data
  int num_frames;
  file >> num_frames;
  tracking_results.resize (num_frames, Eigen::VectorXf::Zero (3 * num_nodes + 3));
  for (size_t f_i = 0; f_i < num_frames; ++f_i)
    for (size_t i = 0; i < tracking_results.front ().rows (); ++i)
      file >> tracking_results[f_i] (i);


  file.close ();

  *skeleton_mesh.mesh_ = *skeleton_mesh.mesh_rest_;
  skeleton_mesh.skeleton_ = skeleton_mesh.skeleton_rest_;

  return (true);
}
