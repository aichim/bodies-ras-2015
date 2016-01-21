#include "common.h"
#include "bodies_utils.h"
#include "skeleton_mesh.h"
#include "SkeletonIOFbx.h"


void
smoothTracking (const std::vector<Eigen::VectorXf> &vec_in,
                std::vector<Eigen::VectorXf> &vec_out)
{
  vec_out = vec_in;
  for (size_t f_i = 1; f_i < vec_in.size () - 1; ++f_i)
  {
    vec_out[f_i] = (vec_in[f_i - 1] + 2. * vec_in[f_i] + vec_in[f_i + 1]) / 4.;
  }
}



int
main (int argc,
      char **argv)
{
  std::string fbx_path = "";
  pcl::console::parse_argument (argc, argv, "-fbx", fbx_path);

  std::string mesh_path = "";
  pcl::console::parse_argument (argc, argv, "-mesh", mesh_path);

  std::string results_path = "";
  pcl::console::parse_argument (argc, argv, "-results", results_path);

  std::string out_dir = "";
  pcl::console::parse_argument (argc, argv, "-out_dir", out_dir);


  af::Mesh<float>::Ptr mesh (new af::Mesh<float> ());
  af::Mesh<float>::readMeshOBJ (mesh_path, *mesh);

  af::SkeletonMesh skeleton_mesh;
  skeleton_mesh.setRestMesh (mesh);
  af::Skeleton skeleton_rest;
  af::SkeletonIOFbx (skeleton_rest).read (fbx_path);
  skeleton_mesh.setRestSkeleton (skeleton_rest);

  std::vector<Eigen::VectorXf> tracking_results;
  af::loadTrackingAndModelingResults (results_path, skeleton_mesh, tracking_results);

  /// Done reading

  std::vector<Eigen::VectorXf> tracking_results_smoothed;
//  smoothTracking (tracking_results, tracking_results_smoothed);
  tracking_results_smoothed = tracking_results;



  /// Save the resulting meshes
  boost::filesystem::create_directory (out_dir);
  char str[512];
  for (size_t f_i = 0; f_i < tracking_results_smoothed.size (); f_i += 5) //++f_i)
  {
    skeleton_mesh.skeleton_.applyPose (tracking_results_smoothed[f_i]);
    skeleton_mesh.bakeGeometry ();
    skeleton_mesh.mesh_->computeNormalsUsingConnectivity ();

    sprintf (str, "%s/mesh_frame_%04zu.obj", out_dir.c_str (), f_i);
    af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.mesh_);

    sprintf (str, "%s/skeleton_frame_%04zu.obj", out_dir.c_str (), f_i);
    af::Mesh<float>::writeMeshOBJ (str, *skeleton_mesh.skeleton_.generateMesh ());
  }



  return (0);
}
