#include "common.h"
#include "bodies_utils.h"
#include "skeleton_mesh.h"

#include <pcl/PolygonMesh.h>
#include <pcl/visualization/pcl_visualizer.h>
#include <pcl/PCLPointCloud2.h>
#include <pcl/io/pcd_io.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/integral_image_normal.h>

#include "SkeletonIOFbx.h"

#include <boost/filesystem.hpp>

//#define ENABLE_VIS

af::CameraParams camera_params_rgb, camera_params_depth;


void
convertAFMeshToPCLPolygonMesh (af::Mesh<float>::ConstPtr mesh,
                               pcl::PolygonMesh &polygon_mesh,
                               pcl::PointCloud<pcl::PointXYZRGBA> &pc)
{
  pc.clear ();
  /// Convert the vertices
  for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
  {
    pcl::PointXYZRGBA p;
    p.getVector3fMap () = mesh->vertices_.col (v_i).cast<float> ();
    p.r = 0.; p.g = 255.; p.b = 0.;
    pc.push_back (p);
  }
  pc.width = 1;
  pc.height = pc.points.size ();
  pcl::PCLPointCloud2 pc2;
  pcl::toPCLPointCloud2 (pc, pc2);
  polygon_mesh.cloud = pc2;

  /// Convert the topology
  polygon_mesh.polygons.resize (mesh->tri_vertices_.cols ());
  for (size_t tri_i = 0; tri_i < mesh->tri_vertices_.cols (); ++tri_i)
  {
    pcl::Vertices verts;
    for (size_t i = 0; i < 3; ++i)
      verts.vertices.push_back (mesh->tri_vertices_ (i, tri_i));
    polygon_mesh.polygons[tri_i] = verts;
  }

  std::cerr << "polygonmesh has " << pc2.height * pc2.width << " vertices and " << polygon_mesh.polygons.size () << " triangles.\n";
}


void
findPointCloudContour (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                       std::vector<int> &indices_contour)
{
//  cv::Mat img_debug = cv::Mat::zeros (cloud->height, cloud->width, CV_8UC3);
//  img_debug.setTo (cv::Scalar (255, 255, 255));

  for (size_t v = 1; v < cloud->height - 1; ++v)
    for (size_t u = 1; u < cloud->width - 1; ++u)
      if (pcl_isfinite ((*cloud)[v * 640 + u].z))
      {
        int count_valid = 0;
        for (int du = -1; du <= 1; ++du)
          for (int dv = -1; dv <= 1; ++dv)
            if (pcl_isfinite ((*cloud)[(v + dv) * 640 + u + du].z) &&
                fabs ((*cloud)[(v + dv) * 640 + u + du].z - (*cloud)[v * 640 + u].z) < 0.03)
              count_valid ++;
        if (count_valid < 7)
        {
          indices_contour.push_back (v * 640 + u);
//          cv::circle (img_debug, cv::Point (u, v), 2, cv::Scalar (0, 0, 0));
        }
      }

//  cv::imshow ("img contour depths", img_debug);
}



void
kdtreeCorrespondencesOMP (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud_image,
                          pcl::PointCloud<pcl::Normal>::ConstPtr normals_image,
                          af::Mesh<float>::ConstPtr mesh,
                          const std::vector<bool> &mesh_good_indices,
                          pcl::Correspondences &corresps)
{
  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_mesh (new pcl::PointCloud<pcl::PointXYZRGBA> ());
  pcl::IndicesPtr indices (new std::vector<int> ());
  pcl::PointXYZRGBA p;
  for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
  {
    float angle = (-Eigen::Vector3f::UnitZ ()).dot (mesh->normals_.col (v_i));
    if (mesh_good_indices[v_i] &&
        angle >= cos (85. * M_PI / 180.))
      indices->push_back (v_i);

    p.getVector3fMap () = mesh->vertices_.col (v_i);
    cloud_mesh->push_back (p);
  }
  PCL_ERROR ("kdtreecorrespsomp size: %zu\n", indices->size ());

  pcl::search::KdTree<pcl::PointXYZRGBA> kdtree_mesh;
  kdtree_mesh.setInputCloud (cloud_mesh, indices);

  const size_t sampling_rate = 4;
#pragma omp parallel for
  for (size_t v = 0; v < cloud_image->height; v += sampling_rate)
    for (size_t u = 0; u < cloud_image->width; u += sampling_rate)
      if (pcl_isfinite ((*cloud_image)[v * 640 + u].z) &&
          pcl_isfinite ((*normals_image)[v * 640 + u].normal_x))
      {
        std::vector<int> nn_indices;
        std::vector<float> nn_dists;
        kdtree_mesh.nearestKSearch ((*cloud_image)[v * 640 + u], 1, nn_indices, nn_dists);
#pragma omp critical
        corresps.push_back (pcl::Correspondence (nn_indices.front (), v * 640 + u, nn_dists.front ()));
      }
}



void
contourNormals (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud_image,
                cv::Mat &img_dx, cv::Mat &img_dy,
                cv::Mat &img_debug)
{
  pcl::console::TicToc timer;
  timer.tic ();

  cv::Mat img_cloud = 255 * cv::Mat::ones (cloud_image->height, cloud_image->width, CV_8UC1);
//  img_debug = cv::Mat::zeros (cloud_image->height, cloud_image->width, CV_8UC3);
//  img_debug.setTo (cv::Scalar (255, 255, 255));
#pragma omp parallel for
  for (size_t v = 0; v < cloud_image->height; v += 1)
    for (size_t u = 0; u < cloud_image->width; u += 1)
      if (pcl_isfinite ((*cloud_image)[v * 640 + u].z))
      {
        img_cloud.at<unsigned char> (cv::Point (u, v)) = 0;
//        img_debug.at<cv::Vec3b> (cv::Point (u, v)) = cv::Vec3b (255, 255, 255);
//        cv::circle (img_debug, cv::Point (u, v), 2, cv::Scalar (0, 0, 0));
      }

  cv::Mat img_cloud_blur;
  cv::GaussianBlur (img_cloud, img_cloud_blur, cv::Size (-1, -1), 1.5);

//  cv::imshow ("img cloud blur", img_cloud_blur);


//  cv::Mat img_dx, img_dy;
  cv::Sobel (img_cloud_blur, img_dx, CV_32F, 1, 0);
  cv::Sobel (img_cloud_blur, img_dy, CV_32F, 0, 1);

  double duration = timer.toc ();
  PCL_ERROR ("###TIME### Contour normals: %f ms\n", duration);

  /*
  /// DEBUG Visualization
  cv::Mat img_dx_norm, img_dy_norm;
  cv::normalize (img_dx, img_dx_norm, 0., 255., cv::NORM_MINMAX);
  cv::normalize (img_dy, img_dy_norm, 0., 255., cv::NORM_MINMAX);

  cv::imshow ("body", img_cloud);
  cv::imshow ("body blur", img_cloud_blur);
  cv::imshow ("body dx", img_dx_norm / 255.);
  cv::imshow ("body dy", img_dy_norm / 255.);
  cv::waitKey ();
  */
}


#ifdef ENABLE_VIS
pcl::visualization::PCLVisualizer vis_aux, vis_mesh;
#endif
void
correspsContour (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud_image,
                 af::Mesh<float>::ConstPtr mesh,
                 const std::vector<bool> &mesh_good_indices,
                 pcl::Correspondences &corresps,
                 std::vector<Eigen::Vector3f> &contour_normals)
{
  Eigen::Matrix3f K = camera_params_depth.getK ();
  cv::Mat img_dx, img_dy, img_debug;
  contourNormals (cloud_image, img_dx, img_dy, img_debug);

  pcl::IndicesPtr indices_contour (new std::vector<int> ());
  findPointCloudContour (cloud_image, *indices_contour);
  PCL_ERROR ("Contour has %zu points.\n", indices_contour->size ());
  pcl::search::KdTree<pcl::PointXYZRGBA> kdtree_contour;
  kdtree_contour.setInputCloud (cloud_image, indices_contour);

  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_mesh_contour (new pcl::PointCloud<pcl::PointXYZRGBA> ());
#pragma omp parallel for
  for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
  {
    if (!mesh_good_indices[v_i])
      continue;
    float angle = (-Eigen::Vector3f::UnitZ ()).dot (mesh->normals_.col (v_i));
    if (angle < cos (90. * M_PI / 180.) ||
        angle > cos (70. * M_PI / 180.))
      continue;

    pcl::PointXYZRGBA p_q;
    p_q.getVector3fMap () = mesh->vertices_.col (v_i);
    cloud_mesh_contour->push_back (p_q);

    Eigen::Vector3f normal_mesh_3d = K * mesh->normals_.col (v_i);
    Eigen::Vector2f normal_mesh_2d (normal_mesh_3d (0), normal_mesh_3d (1));
    normal_mesh_2d.normalize ();

    /// Search within the radius, and choose the closest with a good angle
    std::vector<int> nn_indices;
    std::vector<float> nn_dists;

    kdtree_contour.radiusSearch (p_q, 0.1, nn_indices, nn_dists);
//    kdtree_contour.nearestKSearch (p_q, 50, nn_indices, nn_dists);
    int min_nn = -1;
    float min_dist = std::numeric_limits<float>::max ();
    Eigen::Vector2f min_normal;
    for (size_t n_i = 0; n_i < nn_indices.size (); ++n_i)
    {
      /// Get the normal
      int u = nn_indices[n_i] % 640;
      int v = nn_indices[n_i] / 640;
      Eigen::Vector2f normal_image (img_dx.at<float> (cv::Point (u, v)),
                                    img_dy.at<float> (cv::Point (u, v)));
      normal_image.normalize ();
      if (normal_mesh_2d.dot (normal_image) > 0. &&
          nn_dists[n_i] < min_dist)
      {
        min_nn = n_i;
        min_dist = nn_dists[n_i];
        min_normal = normal_image;
      }
    }

#pragma omp critical
{
    if (min_nn != -1)
    {
      corresps.push_back (pcl::Correspondence (v_i, nn_indices[min_nn], nn_dists[min_nn]));
      Eigen::Vector3f normal (min_normal (0), min_normal (1), 0.);
      contour_normals.push_back (normal);
    }
}

    /*
    /// Debug stuff
    cv::Mat img_mesh_debug = img_debug;
    img_mesh_debug.setTo (cv::Scalar (255, 255, 255));
    for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
    {
      if (!mesh_good_indices[v_i])
        continue;
      float angle = (-Eigen::Vector3f::UnitZ ()).dot (mesh->normals_.col (v_i));
      if (angle < cos (90. * M_PI / 180.) ||
          angle > cos (70. * M_PI / 180.))
        continue;

      Eigen::Vector3f p_proj = camera_params_depth.getK () * mesh->vertices_.col (v_i);
      p_proj /= p_proj (2);

      cv::circle (img_mesh_debug, cv::Point (p_proj (0), p_proj (1)), 2, cv::Scalar (0, 0, 0));
    }
    cv::imshow ("mesh contour", img_mesh_debug);
    cv::waitKey ();
*/

/*
    /// Get only the closest neighbor and reject it if bad
    kdtree_contour.nearestKSearch (p_q, 1, nn_indices, nn_dists);

    cloud_mesh_contour->push_back (p_q);

    if ((p_q.getVector3fMap () - (*cloud_image)[nn_indices.front ()].getVector3fMap ()).norm () > 0.05)
      continue;

    /// Get the normal
    int u = nn_indices.front () % 640;
    int v = nn_indices.front () / 640;
    Eigen::Vector2f normal_image (img_dx.at<float> (cv::Point (u, v)),
                                  img_dy.at<float> (cv::Point (u, v)));
    normal_image.normalize ();

    cv::line (img_debug, cv::Point (u, v), cv::Point (u + 10. * normal_image (0),
                                                      v + 10. * normal_image (1)),
              cv::Scalar (0, 0, 255), 1);

    Eigen::Vector3f normal_mesh_3d = K * mesh->normals_.col (v_i);
    Eigen::Vector2f normal_mesh_2d (normal_mesh_3d (0), normal_mesh_3d (1));
    normal_mesh_2d.normalize ();

    if (normal_mesh_2d.dot (normal_image) > 0.)
    {
      corresps.push_back (pcl::Correspondence (v_i, nn_indices.front (), nn_dists.front ()));

      Eigen::Vector3f normal (normal_image (0), normal_image (1), 0.);
      contour_normals.push_back (normal);
    }
    */

  }

//  cv::imshow ("normals", img_debug);
//  cv::waitKey ();


#ifdef ENABLE_VIS
  /// Debugging visualization
  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGBA> color_mesh (cloud_mesh_contour, 255., 0., 0.);
  if (!vis_aux.updatePointCloud<pcl::PointXYZRGBA> (cloud_mesh_contour, color_mesh, "mesh_cloud_contour"))
    vis_aux.addPointCloud<pcl::PointXYZRGBA> (cloud_mesh_contour, color_mesh, "mesh_cloud_contour");

  pcl::ExtractIndices<pcl::PointXYZRGBA> ei;
  ei.setInputCloud (cloud_image);
  ei.setIndices (indices_contour);
  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr cloud_image_contour (new pcl::PointCloud<pcl::PointXYZRGBA> ());
  ei.filter (*cloud_image_contour);

  pcl::visualization::PointCloudColorHandlerCustom<pcl::PointXYZRGBA> color_cloud (cloud_image_contour, 0., 255., 0.);
  if (!vis_aux.updatePointCloud<pcl::PointXYZRGBA> (cloud_image_contour, color_cloud, "image_cloud_contour"))
    vis_aux.addPointCloud<pcl::PointXYZRGBA> (cloud_image_contour, color_cloud, "image_cloud_contour");


  /// Create the pointclouds for the correspondence visualization
  pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_vis_mesh (new pcl::PointCloud<pcl::PointXYZ> ()),
      cloud_vis_image (new pcl::PointCloud<pcl::PointXYZ> ());
  pcl::Correspondences corresps_vis;
  for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
  {
    pcl::PointXYZ p;
    p.getVector3fMap () = mesh->vertices_.col (corresps[c_i].index_query);
    cloud_vis_mesh->push_back (p);

    p.getVector3fMap () = (*cloud_image)[corresps[c_i].index_match].getVector3fMap ();
    cloud_vis_image->push_back (p);
    corresps_vis.push_back (pcl::Correspondence (c_i, c_i, 0.));
  }

  if (!vis_aux.updatePointCloud<pcl::PointXYZ> (cloud_vis_mesh, "mesh_cloud_contour_corresps"))
    vis_aux.addPointCloud<pcl::PointXYZ> (cloud_vis_mesh, "mesh_cloud_contour_corresps");
  if (!vis_aux.updatePointCloud<pcl::PointXYZ> (cloud_vis_image, "image_cloud_contour_corresps"))
    vis_aux.addPointCloud<pcl::PointXYZ> (cloud_vis_image, "image_cloud_contour_corresps");


  if (!vis_aux.updateCorrespondences<pcl::PointXYZ> (cloud_vis_mesh, cloud_vis_image, corresps_vis, "corresps"))
    vis_aux.addCorrespondences<pcl::PointXYZ> (cloud_vis_mesh, cloud_vis_image, corresps_vis, "corresps");

  vis_aux.spin ();
#endif
}


float
trackingSuccessful (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                    af::Mesh<float>::ConstPtr mesh,
                    const af::CameraParams &camera_params)
{
//  Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> cloud_bool (Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic>::Zero (cloud->width, cloud->height));

  Eigen::Matrix3f K = camera_params.getK ();
  Eigen::Vector3f view_dir (-Eigen::Vector3f::UnitZ ());

  /// Project each visible vertex to the image
  int count_cloud_hits = 0,
      count_visible_vertices = 0;
  for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
  {
    if (view_dir.dot (mesh->normals_.col (v_i)) < cos (80. * M_PI / 180.))
      continue;

    Eigen::Vector3f vertex = mesh->vertices_.col (v_i);
    Eigen::Vector3f v_proj = K * vertex;
    v_proj /= v_proj (2);
    int u = static_cast<int> (v_proj (0));
    int v = static_cast<int> (v_proj (1));

    if (u < 0 || u >= 640 ||
        v < 0 || v >= 480)
      continue;

    count_visible_vertices ++;

    if (!pcl_isfinite ((*cloud)[v * 640 + u].z))
      continue;

    count_cloud_hits ++;
//    cloud_bool (u, v) = true;
  }

  float explained_visible_vertices = static_cast<float> (count_cloud_hits) / static_cast<float> (count_visible_vertices);
  PCL_ERROR ("Percentage of explained visible vertices %d / %d: %f\n",
             count_cloud_hits, count_visible_vertices, explained_visible_vertices);

  return (explained_visible_vertices);
}


int
main (int argc,
      char **argv)
{
  std::string dataset_dir = "";
  pcl::console::parse_argument (argc, argv, "-dataset_dir", dataset_dir);

  std::string fbx_path = "";
  pcl::console::parse_argument (argc, argv, "-fbx", fbx_path);

  std::string bs_model_dir = "";
  pcl::console::parse_argument (argc, argv, "-bs_dir", bs_model_dir);

  std::string indices_path = "";
  pcl::console::parse_argument (argc, argv, "-good_points", indices_path);

  std::string pca_path = "";
  pcl::console::parse_argument (argc, argv, "-pca", pca_path);

  std::string out_dir = "";
  pcl::console::parse_argument (argc, argv, "-out_dir", out_dir);




  /// Load the blendshape/PCA model
  af::BlendshapeModel<float>::Ptr bs_model (new af::BlendshapeModel<float> ());
  bs_model->load (bs_model_dir);
  bs_model->convertToTriMesh ();
  bs_model->neutral_mesh_->computeNormalsUsingConnectivity ();

  af::Mesh<float>::Ptr mesh (new af::Mesh<float> ());
  *mesh = *bs_model->neutral_mesh_;

  /// Read the good indices
  std::vector<size_t> mesh_good_indices;
  af::readIndicesFile (indices_path, mesh_good_indices);
  std::vector<bool> mesh_good_indices_bool (mesh->vertices_.cols (), false);
  for (size_t i = 0; i < mesh_good_indices.size (); ++i)
    mesh_good_indices_bool[mesh_good_indices[i]] = true;

  boost::filesystem::create_directory (out_dir);


  /// Load the PCA model
  af::PosePCA pose_pca;
  af::loadPCA (pca_path, pose_pca);


  /// Load the skeleton
  af::SkeletonMesh skeleton_mesh;
  skeleton_mesh.setRestMesh (mesh);
  af::Skeleton skeleton_rest;
  af::SkeletonIOFbx (skeleton_rest).read (fbx_path);
  skeleton_mesh.setRestSkeleton (skeleton_rest);

  /// Load the input RGB-D frames
  std::vector<af::SkeletonState*> skeleton_frames;
  std::vector<pcl::PointCloud<pcl::PointXYZRGBA>::Ptr> point_clouds;
  af::loadBodyTrackingRGBDDataset (dataset_dir, skeleton_frames, point_clouds,
                                   camera_params_rgb, camera_params_depth);

  /// Scale the skeleton using the NiTE features throughout the whole sequence
  af::SkeletonMesh skeleton_mesh_scaled;
  skeleton_mesh_scaled.setRestMesh (mesh);
  skeleton_mesh_scaled.setRestSkeleton (skeleton_rest);
  af::Skeleton::scaleNiteToSLSkeleton (skeleton_frames, skeleton_mesh_scaled);

  /// Debug
  skeleton_mesh_scaled.bakeGeometry ();
  af::Mesh<float>::writeMeshOBJ (out_dir + "/body_original.obj", *skeleton_mesh.mesh_);
  af::Mesh<float>::writeMeshOBJ (out_dir + "/body_scaled.obj", *skeleton_mesh_scaled.mesh_);

  skeleton_mesh.setRestMesh (skeleton_mesh_scaled.mesh_);
  skeleton_mesh.setRestSkeleton (skeleton_mesh_scaled.skeleton_);
  skeleton_mesh.initIndices ();

  af::Mesh<float>::Ptr mesh_skel_test = skeleton_mesh.skeleton_.generateMesh ();
  af::Mesh<float>::writeMeshOBJ (out_dir + "/skeleton_rest.obj", *mesh_skel_test);

  /// HACK - change the neutral of the blendshape model with the scaled mesh
  /// TODO - should do it in a smarter way
  *bs_model->neutral_mesh_ = *skeleton_mesh_scaled.mesh_;
  bs_model->neutral_ = Eigen::Map<Eigen::MatrixXf> (bs_model->neutral_mesh_->vertices_.data (),
                                                    bs_model->neutral_mesh_->vertices_.cols () * 3, 1);

#ifdef ENABLE_VIS
  pcl::visualization::PCLVisualizer vis;
  pcl::visualization::PointCloudColorHandlerRGBAField<pcl::PointXYZRGBA> rgba (point_clouds[0]);
  vis.addPointCloud<pcl::PointXYZRGBA> (point_clouds[0], rgba, "cloud");
  vis.setBackgroundColor (0.5, 0.5, 0.5);
#endif


  pcl::console::TicToc timer, timer_frame;
  double duration = 0.;


  Eigen::VectorXf var_bs (Eigen::VectorXf::Zero (bs_model->num_exps_));
  bool skip_nite_features = false;

  std::vector<Eigen::VectorXf> angles_frame_vec (skeleton_frames.size (), Eigen::VectorXf::Zero (3 * skeleton_mesh.skeleton_.nodes.size () + 3));
  std::vector<float> tracking_scores (skeleton_frames.size (), 0.f);
  float s_prev = 1.;
  float s_cur = 1.;
  float gamma = 0.96;
  for (size_t glob_it = 0; glob_it < 3; ++glob_it)
  {
    s_prev = s_cur;
    PCL_ERROR ("===========> Global iteration %zu\n", glob_it);

    double total_modeling_error = 0.;
    Eigen::MatrixXf M_lhs (Eigen::MatrixXf::Zero (bs_model->num_exps_, bs_model->num_exps_));
    Eigen::VectorXf M_rhs (Eigen::VectorXf::Zero (bs_model->num_exps_));
    Eigen::MatrixXf M_lhs_iter (Eigen::MatrixXf::Zero (bs_model->num_exps_, bs_model->num_exps_));
    Eigen::VectorXf M_rhs_iter (Eigen::VectorXf::Zero (bs_model->num_exps_));

    int count_modeling_frames = 0, count_tracking_frames = 0;

    skeleton_mesh.skeleton_ = skeleton_mesh.skeleton_rest_;

    for (size_t frame_i = 0; frame_i < skeleton_frames.size () - 3; frame_i += 5)
    {
      if (frame_i != 0)
        angles_frame_vec[frame_i] = angles_frame_vec[frame_i - 1];

      timer_frame.tic ();
      PCL_ERROR ("====== Frame %zu ========\n", frame_i);


      skeleton_mesh.skeleton_.computePathsDFS (skeleton_mesh.skeleton_.getRoot ());

      timer.tic ();
//      if (!skip_nite_features)
        skeleton_mesh.skeleton_.IKWithJointPositions (skeleton_frames[frame_i]->joint_pos);
      duration = timer.toc ();
      PCL_ERROR ("###TIME### IK with joint positions took %f ms\n", duration);
      timer.tic ();
      PCL_ERROR ("IS Pose close to neutral -> %d\n", skeleton_mesh.skeleton_.closeToNeutralPose ());
      if (skeleton_mesh.skeleton_.closeToNeutralPose ())
        skip_nite_features = true;

      skeleton_mesh.bakeGeometry ();
      duration = timer.toc ();
      PCL_ERROR ("###TIME### Bake geometry: %f\n", duration);

#ifdef ENABLE_VIS
      pcl::visualization::PointCloudColorHandlerRGBAField<pcl::PointXYZRGBA> rgba (point_clouds[frame_i]);
      vis.updatePointCloud<pcl::PointXYZRGBA> (point_clouds[frame_i], rgba, "cloud");

      pcl::PolygonMesh polygon_mesh;
      pcl::PointCloud<pcl::PointXYZRGBA>::Ptr mesh_cloud (new pcl::PointCloud<pcl::PointXYZRGBA> ());
      convertAFMeshToPCLPolygonMesh (skeleton_mesh.mesh_, polygon_mesh, *mesh_cloud);
      //    if (!vis.updatePolygonMesh (polygon_mesh, "mesh"))
      //      vis.addPolygonMesh (polygon_mesh, "mesh");
      if (!vis.updatePointCloud (mesh_cloud, "mesh_cloud"))
        vis.addPointCloud (mesh_cloud, "mesh_cloud");

      /// Add spheres for the joint positions for debugging purposes
      for (size_t i = 0; i < skeleton_frames[frame_i]->joint_pos.size (); ++i)
      {
        pcl::PointXYZ center;
        center.getVector3fMap () = skeleton_frames[frame_i]->joint_pos[i];
        double radius = 0.03;
        double r = 255.;
        double b = 0.;
        double g = 0.;

        char str[128];
        sprintf (str, "blobs_%zu", i);
        if (!vis.updateSphere (center, radius, r, g, b, str))
          vis.addSphere (center, radius, r, g, b, str);
      }

      for (size_t i = 0; i < skeleton_mesh.skeleton_.nodes.size (); ++i)
      {
        pcl::PointXYZ center;
        center.getVector3fMap () = (skeleton_mesh.skeleton_.nodes[i]->global_transformation * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        double radius = 0.03;
        double r = 0.;
        double b = 255.;
        double g = 0.;

        char str[128];
        sprintf (str, "tracked_%zu", i);
        if (!vis.updateSphere (center, radius, r, g, b, str))
          vis.addSphere (center, radius, r, g, b, str);
      }
#endif


      PCL_ERROR ("IK with SDK joint positions...\n");
//      vis.spin ();

      PCL_ERROR ("Computing normals with %zu vertices...\n", point_clouds[frame_i]->size ());
      timer.tic ();
      /// Compute the point cloud normals, make sure they are oriented properly
      /// TODO use the IntegralImages method - find good parameters for that.
      pcl::NormalEstimationOMP<pcl::PointXYZRGBA, pcl::Normal> normal_estimation;
//      pcl::IntegralImageNormalEstimation<pcl::PointXYZRGBA, pcl::Normal> normal_estimation;
      normal_estimation.setInputCloud (point_clouds[frame_i]);
      normal_estimation.setKSearch (10);
//      normal_estimation.setRectSize (5, 5);
//      normal_estimation.setNormalEstimationMethod (normal_estimation.AVERAGE_3D_GRADIENT);
//      normal_estimation.setMaxDepthChangeFactor(0.02f);
//      normal_estimation.setNormalSmoothingSize(10.0f);
      pcl::PointCloud<pcl::Normal>::Ptr normal_cloud (new pcl::PointCloud<pcl::Normal> ());
      normal_estimation.compute (*normal_cloud);
      duration = timer.toc ();
      PCL_ERROR ("normal cloud size: %d\n", normal_cloud->size ());

      PCL_ERROR ("###TIME### Computing normals took: %f\n", duration);

#ifdef ENABLE_VIS
      vis.removeShape ("normals");
      vis.addPointCloudNormals<pcl::PointXYZRGBA, pcl::Normal> (point_clouds[frame_i], normal_cloud, 10, 0.05, "normals");
#endif


      double p2plane_tracking_error_prev = std::numeric_limits<double>::max ();
      for (size_t reg_iter = 0; reg_iter < 50; ++reg_iter)
      {
        timer.tic ();

        PCL_ERROR ("--- Registration iteration %zu ---\n", reg_iter);
        pcl::console::TicToc timer_small;

        /// Compute the mesh normals
        skeleton_mesh.mesh_->computeNormalsUsingConnectivity ();

        /// Compute correspondences (from the mesh to the point cloud)
        pcl::Correspondences corresps, corresps_contour;

//        projectiveCorrespondences (point_clouds[frame_i], skeleton_mesh.mesh_, corresps);
        kdtreeCorrespondencesOMP (point_clouds[frame_i], normal_cloud, skeleton_mesh.mesh_, mesh_good_indices_bool, corresps);
        std::vector<Eigen::Vector3f> normals_contour;
        correspsContour (point_clouds[frame_i], skeleton_mesh.mesh_, mesh_good_indices_bool, corresps_contour, normals_contour);
        duration = timer_small.toc ();


        /// Filter the correspondences
        pcl::Correspondences corresps_filtered;
        for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
        {
          /// Threshold on the distance
          float dist = (skeleton_mesh.mesh_->vertices_.col (corresps[c_i].index_query) -
                         point_clouds[frame_i]->at (corresps[c_i].index_match).getVector3fMap ()).norm ();
          if (dist > 0.08)
            continue;

          /// Throw away normal incompatibilities
          float cos_angle = skeleton_mesh.mesh_->normals_.col (corresps[c_i].index_query).dot (
                            normal_cloud->at (corresps[c_i].index_match).getNormalVector3fMap ());
//          if (fabs (cos_angle) < cos (60. * M_PI / 180.))
//            continue;

          corresps_filtered.push_back (corresps[c_i]);
        }


        PCL_ERROR ("Correspondences: %zu, after filtering %zu\n", corresps.size (), corresps_filtered.size ());
        PCL_ERROR ("Correspondences contour: %zu\n", corresps_contour.size ());

        duration = timer.toc ();
        PCL_ERROR ("###TIME### Computing correspondences took %f ms\n", duration);

        /// Solve the IK problem
        timer.tic ();
//        corresps_filtered.resize (1);

//        double error = skeleton_mesh.IKWithPointToPlane (point_clouds[frame_i], normal_cloud, corresps_filtered, 5e-1);
//        double error = skeleton_mesh.IKWithPointToPoint (point_clouds[frame_i], corresps_filtered, 5e-1);

#ifdef ENABLE_VIS
        if (!vis.updateCorrespondences<pcl::PointXYZRGBA> (mesh_cloud, point_clouds[frame_i], corresps_filtered, "corresps"))
          vis.addCorrespondences<pcl::PointXYZRGBA> (mesh_cloud, point_clouds[frame_i], corresps_filtered, "corresps");
        vis.spin ();
#endif

        double error = skeleton_mesh.IKCustom (point_clouds[frame_i], normal_cloud,
                                               corresps_filtered, corresps_contour, normals_contour,
                                               skeleton_frames[frame_i]->joint_pos,
                                               angles_frame_vec[frame_i],
                                               pose_pca,
                                               angles_frame_vec[frame_i],
                                               1., 3., 0./*4e-1*/, 5e-2, 1., 1e-1, 1e-3);
        duration = timer.toc ();
        PCL_ERROR ("###TIME### IK with point2plane took %f ms.\n", duration);

#ifdef ENABLE_VIS
        /// Display the registered mesh
        skeleton_mesh.bakeGeometry ();
        convertAFMeshToPCLPolygonMesh (skeleton_mesh.mesh_, polygon_mesh, *mesh_cloud);
        //      if (!vis.updatePolygonMesh (polygon_mesh, "mesh"))
        //        vis.addPolygonMesh (polygon_mesh, "mesh");
        if (!vis.updatePointCloud (mesh_cloud, "mesh_cloud"))
          vis.addPointCloud (mesh_cloud, "mesh_cloud");
        vis.spin ();
#endif

        PCL_ERROR ("After refinement\n");
        PCL_ERROR ("---------> Registration error: %f\n", error);

        if (reg_iter > 1 &&
            (p2plane_tracking_error_prev - error) / p2plane_tracking_error_prev < 5e-2) //1e-2)
        {
          duration = timer_frame.toc ();
          PCL_ERROR ("###TIME### Frame processing: %f\n", duration);

#ifdef ENABLE_VIS
          if (!vis.updatePointCloud (mesh_cloud, "mesh_cloud"))
            vis.addPointCloud (mesh_cloud, "mesh_cloud");

          /// Display the correspondences
          if (!vis.updateCorrespondences<pcl::PointXYZRGBA> (mesh_cloud, point_clouds[frame_i], corresps_filtered, "corresps"))
            vis.addCorrespondences<pcl::PointXYZRGBA> (mesh_cloud, point_clouds[frame_i], corresps_filtered, "corresps");

          char name[64];
          sprintf (name, "iteration_%zu frame_%04zu", glob_it, frame_i);
          vis.setWindowName (name);
          vis.spin ();
#endif

          /// Check how good our tracking was
          float tracking_score = trackingSuccessful (point_clouds[frame_i], skeleton_mesh.mesh_, camera_params_depth);
          tracking_scores[frame_i] = tracking_score;

          count_tracking_frames ++;

          /// Accumulate the constraints for modeling
          if (tracking_score > 55. / 100.)
          {
            timer.tic ();
            M_lhs_iter.setZero ();
            M_rhs_iter.setZero ();
            total_modeling_error += skeleton_mesh.accumulateModelingConstraints (point_clouds[frame_i], normal_cloud, corresps_filtered,
                                                                                 corresps_contour, normals_contour,
                                                                                 bs_model,
                                                                                 var_bs,
                                                                                 M_lhs_iter, M_rhs_iter,
                                                                                 skeleton_frames.size (),
                                                                                 1, 3.);

            s_cur = gamma * s_prev + 1.;
            M_lhs = gamma * s_prev / s_cur * M_lhs + 1. / s_cur * M_lhs_iter;
            M_rhs = gamma * s_prev / s_cur * M_rhs + 1. / s_cur * M_rhs_iter;
            s_prev = s_cur;
            af::SkeletonMesh::solveModelingProblemGS (M_lhs, M_rhs, var_bs, 7e-6); //5e-3);
            bs_model->evaluate (var_bs, skeleton_mesh.mesh_rest_->vertices_);

            /// Print activated bs
            PCL_ERROR ("Activated modeling bs:\n");
            for (size_t i = 0; i < var_bs.rows (); ++i)
              if (fabs (var_bs (i)) > 0.05)
                PCL_ERROR ("  %s - %f\n", bs_model->exp_names_[i].c_str (), var_bs (i));

            count_modeling_frames ++;

            duration = timer.toc ();
            PCL_ERROR ("###TIME### Accumulating modeling constraints took %f ms.\n", duration);


#ifdef ENABLE_VIS
            pcl::PolygonMesh polygon_mesh_vis;
            pcl::PointCloud<pcl::PointXYZRGBA> pc_vis;
            convertAFMeshToPCLPolygonMesh (skeleton_mesh.mesh_rest_, polygon_mesh_vis, pc_vis);

            if (!vis_mesh.updatePolygonMesh (polygon_mesh_vis, "mesh updated"))
              vis_mesh.addPolygonMesh (polygon_mesh_vis, "mesh updated");

            cv::Mat img_bs_weights;
            af::displayBlendshapeWeights (var_bs, img_bs_weights);
            cv::imshow ("bs weights", img_bs_weights);
            cv::waitKey ();
#endif
          }


          /// Some debug stuff
          /// Save each frame skeleton and mesh
          char str[512];
          sprintf (str, "%s/frame_mesh_iter_%zu_frame_%04zu.obj", out_dir.c_str (), glob_it, frame_i);
          skeleton_mesh.mesh_->computeNormalsUsingConnectivity ();
          af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_);
          sprintf (str, "%s/frame_skeleton_iter_%zu_frame_%04zu.obj", out_dir.c_str (), glob_it, frame_i);
          af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.skeleton_.generateMesh ());

          sprintf (str, "%s/mesh_opt_%zu_%zu.obj", out_dir.c_str (), glob_it, frame_i);
          af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_rest_);

          cv::Mat img_bs_weights;
          af::displayBlendshapeWeights (var_bs, img_bs_weights);
          sprintf (str, "%s/bs_weights_%zu_%zu.png", out_dir.c_str (), glob_it, frame_i);
          cv::imwrite (str, img_bs_weights);


          break;
        }

        p2plane_tracking_error_prev = error;
      }
    }

    PCL_ERROR ("Total modeling error accumulated at iteration %zu: %f\n",
               glob_it, total_modeling_error);
    PCL_ERROR ("### Tracked frames %d, modeling frames %d\n", count_tracking_frames, count_modeling_frames);

    /// Print activated bs
    PCL_ERROR ("Activated modeling bs:\n");
    for (size_t i = 0; i < var_bs.rows (); ++i)
      if (fabs (var_bs (i)) > 0.05)
        PCL_ERROR ("  %s - %f\n", bs_model->exp_names_[i].c_str (), var_bs (i));


    char str[512];
    cv::Mat img_bs_weights;
    af::displayBlendshapeWeights (var_bs, img_bs_weights);
    sprintf (str, "%s/bs_weights_%zu.png", out_dir.c_str (), glob_it);
    cv::imwrite (str, img_bs_weights);

//    cv::imshow ("bs weights", img_bs_weights);
//    cv::waitKey ();



    /// Apply the bs weights on the scaled mesh.
    bs_model->evaluate (var_bs, skeleton_mesh.mesh_rest_->vertices_);
    sprintf (str, "%s/body_modeled_iteration_%zu.obj", out_dir.c_str (), glob_it);
    skeleton_mesh.mesh_rest_->computeNormalsUsingConnectivity ();
    af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_rest_);
    std::cerr << "resulting bs weights: " << var_bs.transpose () << std::endl;

    sprintf (str, "%s/bs_cloud_text_%zu.txt", out_dir.c_str (), glob_it);
    af::saveTextForWordCloud (str, bs_model, var_bs);
  }






  ///// Texture the body
  const int TEXTURE_SIZE = 2048;
  std::vector<af::Texel> texels;
  PCL_ERROR ("Compute texel geometry.\n");
  af::Mesh<float>::computeTexelGeometry (*skeleton_mesh.mesh_rest_,
                                         texels, TEXTURE_SIZE);
  cv::Mat texture_acc (cv::Mat::zeros (TEXTURE_SIZE, TEXTURE_SIZE, CV_32FC3));
  Eigen::MatrixXf texture_weights_acc (Eigen::MatrixXf::Zero (TEXTURE_SIZE, TEXTURE_SIZE));

  /// Go through each frame and accumulate texture colors
  for (size_t frame_i = 0; frame_i < skeleton_frames.size () - 3; frame_i += 5)
  {
    PCL_ERROR ("Texture processing frame %zu.\n", frame_i);

    /// Skip frame if the tracking score was quite bad.
    if (tracking_scores[frame_i] < 55. / 100.)
      continue;

    skeleton_mesh.skeleton_.applyPose (angles_frame_vec[frame_i]);
    skeleton_mesh.bakeGeometry ();
    skeleton_mesh.mesh_->computeNormalsUsingConnectivity ();

    /// Save the mesh for debugging purposes
//    char str[512];
//    sprintf (str, "mesh_frame_%04zu.obj", frame_i);
//    af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_);

    af::Mesh<float>::ConstPtr mesh_frame = skeleton_mesh.mesh_;

    PCL_ERROR ("Going through each texel\n");

    pcl::PointCloud<pcl::PointXYZ>::Ptr pcl_debug (new pcl::PointCloud<pcl::PointXYZ> ());
    cv::Mat img_debug (cv::Mat::zeros (480, 640, CV_8UC3));
    for (size_t v = 0; v < (*point_clouds[frame_i]).height; v += 1)
      for (size_t u = 0; u < (*point_clouds[frame_i]).width; u += 1)
      {
        int index = v * 640 + u;
        if (pcl_isfinite ((*point_clouds[frame_i])[index].z))
          img_debug.at<cv::Vec3b> (cv::Point (u, v)) = cv::Vec3b ((*point_clouds[frame_i])[index].b, (*point_clouds[frame_i])[index].g, (*point_clouds[frame_i])[index].r);
      }

    for (size_t t_i = 0; t_i < texels.size (); ++t_i)
    {
      Eigen::Vector3f n = texels[t_i].bary (0) * mesh_frame->normals_.col (mesh_frame->tri_vertices_ (0, texels[t_i].tri_id)) +
                          texels[t_i].bary (1) * mesh_frame->normals_.col (mesh_frame->tri_vertices_ (1, texels[t_i].tri_id)) +
                          texels[t_i].bary (2) * mesh_frame->normals_.col (mesh_frame->tri_vertices_ (2, texels[t_i].tri_id));
      n.normalize ();
      float cos_angle = (-Eigen::Vector3f::UnitZ ()).dot (n);

      if (cos_angle < cos (60. * M_PI / 180.))
        continue;

      /// Use the squared as the weight
      cos_angle *= cos_angle;


      Eigen::Vector3f v = texels[t_i].bary (0) * mesh_frame->vertices_.col (mesh_frame->tri_vertices_ (0, texels[t_i].tri_id)) +
                          texels[t_i].bary (1) * mesh_frame->vertices_.col (mesh_frame->tri_vertices_ (1, texels[t_i].tri_id)) +
                          texels[t_i].bary (2) * mesh_frame->vertices_.col (mesh_frame->tri_vertices_ (2, texels[t_i].tri_id));
      Eigen::Vector3f v_proj = camera_params_depth.getK () * v; //camera_params_rgb.getK () * camera_params_rgb.pose * v;
      v_proj /= v_proj (2);

      if (v_proj (0) < 0. || v_proj (1) < 0. ||
          v_proj (0) > 640. || v_proj (1) > 480.)
        continue;


      /// TODO bilinear interpolation on cv::Mats instead of PCL datastructures
      int index = static_cast<int> (v_proj (1)) * 640 + static_cast<int> (v_proj (0));
//      if (!pcl_isfinite ((*point_clouds[frame_i])[index].z))
//        continue;

      pcl::PointXYZ p;
      p.getVector3fMap () = v;
      pcl_debug->push_back (p);

      cv::Vec3f color_proj (cos_angle * static_cast<float> ((*point_clouds[frame_i])[index].b),
                            cos_angle * static_cast<float> ((*point_clouds[frame_i])[index].g),
                            cos_angle * static_cast<float> ((*point_clouds[frame_i])[index].r));
      texture_acc.at<cv::Vec3f> (cv::Point (texels[t_i].u, texels[t_i].v)) += color_proj;
      texture_weights_acc (texels[t_i].u, texels[t_i].v) += cos_angle;
      /// Debugging
      img_debug.at<cv::Vec3b> (cv::Point (v_proj (0), v_proj (1))) = cv::Vec3b (0, 0, 255); //cv::Vec3b ((*point_clouds[frame_i])[index].b, (*point_clouds[frame_i])[index].g, (*point_clouds[frame_i])[index].r);
    }

//    sprintf (str, "cloud_frame_%04zu.pcd", frame_i);
//    pcl::io::savePCDFileBinary (str, *pcl_debug);
//    cv::imshow ("img debug", img_debug);
//    cv::waitKey ();
  }

  /// Normalize the accumulated colors
  PCL_ERROR ("Normalizing the accumulated texture.\n");
  for (size_t t_i = 0; t_i < texels.size (); ++t_i)
  {
    cv::Vec3f col = texture_acc.at<cv::Vec3f> (cv::Point (texels[t_i].u, texels[t_i].v));
    col[0] /= texture_weights_acc (texels[t_i].u, texels[t_i].v);
    col[1] /= texture_weights_acc (texels[t_i].u, texels[t_i].v);
    col[2] /= texture_weights_acc (texels[t_i].u, texels[t_i].v);

    texture_acc.at<cv::Vec3f> (cv::Point (texels[t_i].u, texels[t_i].v)) = col;
  }

  PCL_ERROR ("Save the resulting texture.\n");
  cv::imwrite (out_dir + "/texture.png", texture_acc);


  return (0);
}
