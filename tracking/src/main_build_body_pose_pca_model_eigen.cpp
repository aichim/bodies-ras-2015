#include "common.h"
#include "bodies_utils.h"
#include "SkeletonIOFbx.h"

#include "skeleton_mesh.h"

#include <pcl/visualization/pcl_visualizer.h>




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
}

int
main (int argc,
      char **argv)
{
  std::string dataset_dir = "";
  pcl::console::parse_argument (argc, argv, "-dataset_dir", dataset_dir);

  std::string fbx_path = "";
  pcl::console::parse_argument (argc, argv, "-fbx", fbx_path);

  std::string obj_path = "";
  pcl::console::parse_argument (argc, argv, "-obj", obj_path);

  std::string out_path = "";
  pcl::console::parse_argument (argc, argv, "-out", out_path);


  std::vector<af::SkeletonState*> skeleton_frames;
  std::vector<pcl::PointCloud<pcl::PointXYZRGBA>::Ptr> point_clouds;
  af::CameraParams camera_params_rgb, camera_params_depth;
  af::loadBodyTrackingRGBDDataset (dataset_dir, skeleton_frames, point_clouds,
                                   camera_params_rgb, camera_params_depth);


  af::Skeleton skeleton_rest;
  af::SkeletonIOFbx (skeleton_rest).read (fbx_path);


  af::Mesh<float> mesh_aux;
  af::Mesh<float>::readMeshOBJ (obj_path, mesh_aux);
  af::Mesh<float>::Ptr mesh (new af::Mesh<float> ());
  af::Mesh<float>::convertToTriangleMesh (mesh_aux, *mesh);
  mesh->computeNormalsUsingConnectivity ();


  af::SkeletonMesh skeleton_mesh;
  skeleton_mesh.setRestMesh (mesh);
  skeleton_mesh.setRestSkeleton (skeleton_rest);

  af::SkeletonMesh skeleton_mesh_scaled;
  skeleton_mesh_scaled.setRestMesh (mesh);
  skeleton_mesh_scaled.setRestSkeleton (skeleton_rest);
  af::Skeleton::scaleNiteToSLSkeleton (skeleton_frames, skeleton_mesh_scaled);

  skeleton_mesh_scaled.bakeGeometry ();
  af::Mesh<float>::writeMeshOBJ ("body_original.obj", *skeleton_mesh.mesh_);
  af::Mesh<float>::writeMeshOBJ ("body_scaled.obj", *skeleton_mesh_scaled.mesh_);

  skeleton_mesh.setRestMesh (skeleton_mesh_scaled.mesh_);
  skeleton_mesh.setRestSkeleton (skeleton_mesh_scaled.skeleton_);
  skeleton_mesh.initIndices ();
  skeleton_mesh.bakeGeometry ();
  af::Mesh<float>::writeMeshOBJ ("body_scaled_copy.obj", *skeleton_mesh.mesh_);



//  skeleton_frames.resize (60);

  int count_num_valid_frames = 0;
  for (size_t f_i = 0; f_i < skeleton_frames.size (); ++f_i)
    if (skeleton_frames[f_i]->joint_pos.size () != 0)
      count_num_valid_frames ++;
  PCL_ERROR ("Valid frames: %d / %zu\n", count_num_valid_frames, skeleton_frames.size ());
  int num_joints = skeleton_rest.nodes.size ();

  skeleton_mesh.skeleton_.computePathsDFS (skeleton_mesh.skeleton_.getRoot ());


  /// TODO user interaction to choose the correctly aligned skeletons and skip the bad ones

  pcl::visualization::PCLVisualizer vis;
  pcl::PolygonMesh polygon_mesh;
  pcl::PointCloud<pcl::PointXYZRGBA>::Ptr mesh_cloud (new pcl::PointCloud<pcl::PointXYZRGBA> ());


  Eigen::MatrixXd data_matrix (Eigen::MatrixXd::Zero (count_num_valid_frames, 3 * (num_joints - 1)));
  int matrix_col = 0;
  for (size_t frame_i = 0; frame_i < skeleton_frames.size (); ++frame_i)
  {
    PCL_ERROR ("Frame %zu / %zu\n", frame_i, skeleton_frames.size ());
    if (skeleton_frames[frame_i]->joint_pos.size () == 0)
      continue;

    /// Compute the IK angles
    skeleton_mesh.skeleton_.IKWithJointPositions (skeleton_frames[frame_i]->joint_pos);
    skeleton_mesh.bakeGeometry ();

    /// Display the result and wait for the key press
    convertAFMeshToPCLPolygonMesh (skeleton_mesh.mesh_, polygon_mesh, *mesh_cloud);
    if (!vis.updatePolygonMesh (polygon_mesh, "mesh"))
      vis.addPolygonMesh (polygon_mesh, "mesh");
    if (!vis.updatePointCloud (point_clouds[frame_i], "pcl"))
      vis.addPointCloud (point_clouds[frame_i], "pcl");

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

//    vis.spin ();
    vis.spinOnce ();


    /// Ignore the rotation of the first node, as it's the global rotation of the body
    for (size_t n_i = 1; n_i < skeleton_mesh.skeleton_.nodes.size (); ++n_i)
    {
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (0));
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (1));
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (2));
      data_matrix.block <1, 3> (matrix_col, 3 * (n_i - 1)) = skeleton_mesh.skeleton_.nodes[n_i]->angles.transpose ().cast<double> ();
    }

    matrix_col ++;
  }

  PCL_ERROR ("Built the data matrix.\n");

  PCL_ERROR ("Sum of all the elements is not zero, need to de-mean: %f\n", data_matrix.sum ());

  /// Compute the mean for each column and subtract it
  Eigen::VectorXd mean (data_matrix.cols ());
  for (size_t d_i = 0; d_i < data_matrix.cols (); ++d_i)
  {
    mean (d_i) = data_matrix.col (d_i).sum () / static_cast<double> (data_matrix.rows ());
    data_matrix.col (d_i) -= Eigen::VectorXd::Ones (data_matrix.rows ()) * mean (d_i);
  }

  /// Check
  PCL_ERROR ("Sum of all the elements should be zero: %f\n", data_matrix.sum ());

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensolver;
  eigensolver.compute (data_matrix.transpose () * data_matrix / static_cast<double> (count_num_valid_frames - 1));
  Eigen::MatrixXd U = eigensolver.eigenvectors ();
  Eigen::VectorXd evs = eigensolver.eigenvalues ();

  PCL_INFO ("sizes: U (%zu, %zu), eigen values %zu %zu\n",
            U.rows (), U.cols (), evs.rows (), evs.cols ());
  PCL_INFO ("Eigenvalues: ");
  double ev_sum = 0.;
  for (int i = evs.rows () - 1; i >= 0; --i)
    ev_sum += evs (i);

  double ev_partial_sum = 0.;
  for (int i = evs.rows () - 1; i >= 0; --i)
  {
    ev_partial_sum += evs (i);
    PCL_INFO ("%zu - %f, it covers %.2f%% of the data explanation, and cummulated so far %.2f%%\n",
              evs.rows () - 1 - i, evs (i), evs (i) / ev_sum * 100., ev_partial_sum / ev_sum * 100.);
  }
  PCL_INFO ("\n");


  /// Create the modes
  const int n_modes = 10;
  Eigen::MatrixXd modes = U.block (0, U.cols () - n_modes, U.rows (), n_modes);
  std::cerr << "modes:\n" << modes << std::endl;
  std::cerr << "mean: " << mean.transpose () << std::endl;

  Eigen::VectorXd std_devs = evs.block (evs.rows () - n_modes, 0, n_modes, 1);

  std::cerr << "std devs: " << std_devs.transpose () << std::endl;


  /// Save things for debugging purposes
  /// Save the mean pose
  skeleton_mesh.skeleton_.applyAngles (mean);
  skeleton_mesh.bakeGeometry ();
  af::Mesh<float>::writeMeshOBJ ("body_mean_pose.obj", *skeleton_mesh.mesh_);

  /// For each mode, save the corresponding meshes with +-1, +-2, +-3 std dev
  char str[512];
  for (size_t mode_i = 0; mode_i < n_modes; ++mode_i)
  {
    for (int i = -5; i <= 5; i += 2)
    {
      Eigen::VectorXd angles = mean + modes.col (mode_i) * std_devs (mode_i) * static_cast<double> (i);
      skeleton_mesh.skeleton_.applyAngles (angles);
      skeleton_mesh.bakeGeometry ();
      sprintf (str, "body_mode_%zu_%d.obj", mode_i, i);
      af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_);
    }
  }


  af::PosePCA pca (mean.cast<float> (),
                   std_devs.cast<float> (),
                   modes.cast<float> ());
  af::savePCA (out_path, pca);


  /// debugging the I/O
  af::PosePCA pca2;
  af::loadPCA (out_path, pca2);
  std::cerr << pca2.modes <<std::endl;

  std::cerr << "mean diff: " << (pca.mean - pca2.mean).squaredNorm () << std::endl;
  std::cerr << "std dev diff: " << (pca.std_devs - pca2.std_devs).squaredNorm () << std::endl;
  std::cerr << "modes diff: " << (pca.modes - pca2.modes).squaredNorm () << std::endl;

  std::cerr << modes.transpose () * modes << std::endl;

  /// Test the PCA error with the training set
  double total_reconstruction_error = 0.;
  for (size_t f_i = 0; f_i < skeleton_frames.size (); ++f_i)
  {
    if (skeleton_frames[f_i]->joint_pos.size () == 0)
      continue;

    skeleton_mesh.skeleton_.IKWithJointPositions (skeleton_frames[f_i]->joint_pos);
    Eigen::VectorXf vars_abs (Eigen::VectorXf::Zero (num_joints * 3 + 3));
    for (size_t n_i = 0; n_i < skeleton_mesh.skeleton_.nodes.size (); ++n_i)
    {
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (0));
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (1));
      af::clampAngles (skeleton_mesh.skeleton_.nodes[n_i]->angles (2));

      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_mesh.skeleton_.nodes[n_i]->angles.cast<float> ();
    }
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_mesh.skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    Eigen::VectorXf vars_abs_unproj_diff = vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pca.mean -
                                           pca.modes * pca.modes.transpose () *
                                           (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pca.mean);
//    for (size_t i = 0; i < pca.std_devs.rows (); ++i)
//      vars_abs_unproj_diff (i) /= pca.std_devs (i);

    std::cerr << "vars abs: " << vars_abs.transpose () << "\nprojected: " << vars_abs_unproj_diff.transpose () << std::endl;

    PCL_ERROR ("pca residual frame %zu: %f\n", f_i, vars_abs_unproj_diff.squaredNorm ());


    total_reconstruction_error += vars_abs_unproj_diff.squaredNorm ();

    /// Display the result and wait for the key press
    skeleton_mesh.bakeGeometry ();
//    convertAFMeshToPCLPolygonMesh (skeleton_mesh.mesh_, polygon_mesh, *mesh_cloud);
//    if (!vis.updatePolygonMesh (polygon_mesh, "mesh"))
//      vis.addPolygonMesh (polygon_mesh, "mesh");
//    vis.spinOnce ();
//    vis.spin ();
  }
  total_reconstruction_error /= static_cast<double> (count_num_valid_frames);

  PCL_ERROR ("Total reconstruction error: %f\n", total_reconstruction_error);

  return (0);
}
