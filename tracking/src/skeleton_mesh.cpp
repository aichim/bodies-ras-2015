#include "skeleton_mesh.h"

//#define DEBUG_JACOBIANS


af::SkeletonMesh::SkeletonMesh ()
{
}

void
af::SkeletonMesh::bakeGeometry ()
{
//  PCL_ERROR ("Baking in the geometry...\n");
  skeleton_rest_.computeGlobalTransformations ();
  skeleton_.computeGlobalTransformations ();

  /// Set the current mesh values to 0
  *mesh_ = *mesh_rest_;
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    for (size_t v_i = 0; v_i < skeleton_.nodes[n_i]->vertex_indices.size (); ++v_i)
      mesh_->vertices_.col (skeleton_.nodes[n_i]->vertex_indices[v_i]) = Eigen::Vector3f::Zero ();
//  PCL_ERROR ("The skeleton has %zu nodes.\n", skeleton_.nodes.size ());

  /// Aggregate the new mesh vertices
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
  {
    Eigen::Matrix4f transf_diff = skeleton_.nodes[n_i]->global_transformation * skeleton_rest_.nodes[n_i]->global_transformation.inverse ();
//    std::cerr << "Global translation of node " << skeleton_.nodes[n_i]->name << " is: " << skeleton_.nodes[n_i]->global_transformation.block<3, 1> (0, 3).transpose () << std::endl;
//    PCL_ERROR ("Node %zu %s has %zu vertices.\n", n_i, skeleton_.nodes[n_i]->name.c_str (), skeleton_.nodes[n_i]->vertex_indices.size ());
    for (size_t v_i = 0; v_i < skeleton_.nodes[n_i]->vertex_indices.size (); ++v_i)
    {
      Eigen::Vector3f pos = mesh_rest_->vertices_.col (skeleton_.nodes[n_i]->vertex_indices[v_i]);
      Eigen::Vector3f pos_tr = Eigen::Affine3f (transf_diff) * pos;
      mesh_->vertices_.col (skeleton_.nodes[n_i]->vertex_indices[v_i]) += skeleton_.nodes[n_i]->vertex_weights[v_i] * pos_tr;
    }
  }
}


/*
void
af::SkeletonMesh::IKWithPointToPlane (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                      pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                      const pcl::Correspondences &corresps)
{
  const double weight_damping = 1e1;
  int num_joints = skeleton_.nodes.size ();

  for (size_t iter_i = 0; iter_i < 25; ++iter_i)
  {
    bakeGeometry ();

    Eigen::MatrixXd J (Eigen::MatrixXd::Zero (corresps.size (), num_joints * 3 + 3));
    Eigen::VectorXd b (Eigen::VectorXd::Zero (corresps.size ()));

    /// Set up the constraints
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      Eigen::Vector3d point = cloud->at (corresps[c_i].index_match).getVector3fMap ().cast<double> ();
      Eigen::Vector3d normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ().cast<double> ();
      Eigen::Vector4d mesh_point_rest;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point_rest (3) = 1.;
      Eigen::Vector4d mesh_point;
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      mesh_point (3) = 1.;

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        double weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        Eigen::Matrix4d transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id; //skeleton_.map_name_to_num_[node_on_path->name];
          Eigen::Matrix4d T_prev, T_next;
          T_prev = node_on_path->global_transformation;
          T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

          Eigen::Matrix4d jac = T_prev * Skeleton::jacrot (T_next * transform_joint_rest_inv * mesh_point_rest);
          J.block<1, 3> (c_i, 3 * node_id) += weight_joint * normal.transpose () * jac.block<3, 3> (0, 0);
        }


        Eigen::Matrix4d transform_joint = node_end->global_transformation;
        Eigen::Matrix4d transf_diff = transform_joint * transform_joint_rest_inv;

        b (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3d::Identity ();
      }

//      b (c_i) = normal.transpose () * (mesh_point.block<3, 1> (0, 0) - point);
    }

    PCL_ERROR ("Registration error at inner iteration %zu: %f\n", iter_i, b.squaredNorm ());

    /// Solve the system
    Eigen::VectorXd var_params = (J.transpose () * J +
                                  weight_damping * Eigen::MatrixXd::Identity (3 * num_joints + 3, 3 * num_joints + 3)).ldlt ().solve (-J.transpose () * b);

    /// Put back the parameters
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    {
      int node_id = n_i; //skeleton_.map_name_to_num_[skeleton_.nodes[n_i]->name];
      Eigen::Matrix3d rot_inc = (Eigen::AngleAxisd (var_params (3 * node_id + 0), Eigen::Vector3d::UnitX ()) *
                                 Eigen::AngleAxisd (var_params (3 * node_id + 1), Eigen::Vector3d::UnitY ()) *
                                 Eigen::AngleAxisd (var_params (3 * node_id + 2), Eigen::Vector3d::UnitZ ())).matrix ();

      skeleton_.nodes[n_i]->local_transformation.block<3, 3> (0, 0) = skeleton_.nodes[n_i]->local_transformation.block<3, 3> (0, 0) * rot_inc;
    }
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += var_params.block<3, 1> (3 * num_joints, 0);
  }

  std::cerr << "HipCenter translation " << skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) << std::endl;
}
*/

float
af::SkeletonMesh::IKWithPointToPlane (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                      pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                      const pcl::Correspondences &corresps,
                                      const float weight_reg_zero)
{
  /*
  /// DEBUG - paths from root
  PCL_ERROR ("DEBUG PATHS FROM ROOT:\n");
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
  {
    PCL_ERROR ("   Path from root for node %s:\n", skeleton_.nodes[n_i]->name.c_str ());
    for (size_t n_j = 0; n_j < skeleton_.nodes[n_i]->path_from_root.size (); ++n_j)
      PCL_ERROR ("      %s\n", skeleton_.nodes[n_i]->path_from_root[n_j]->name.c_str ());
    PCL_ERROR ("\n");
  }
  */

  int num_joints = skeleton_.nodes.size ();

  /// Initialize some matrices
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1.;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b (Eigen::VectorXf::Zero (corresps.size ()));

  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));

  for (size_t iter_i = 0; iter_i < 150; ++iter_i)
  {
    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    Eigen::VectorXf vars_abs (num_joints * 3 + 3);
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    J.setZero ();
    b.setZero ();

    /// Set up the constraints
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

          J (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

    /// Solve the system
    lhs = J.transpose () * J +
          weight_reg_zero * J_reg_zero;
    rhs = J.transpose () * (-b) +
          weight_reg_zero * J_reg_zero * (-vars_abs);
    float total_error = b.squaredNorm () + (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ();

    if ((total_error_prev - total_error) / total_error_prev < 1e-3)
    {
      weight_damping *= 10.;
//      PCL_ERROR ("-> damping increased to ");
      if (weight_damping > 1e8)
        break;

      /// Take back the angles and the transformation
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars.block<3, 1> (3 * num_joints, 0);
    }
    else
    {
//      PCL_ERROR ("Registration error at inner iteration %zu: %f\n", iter_i, b.squaredNorm ());
//      PCL_ERROR ("Reg zero error: %f\n", (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ());
//      PCL_ERROR ("Total error: %f\n", total_error);

      weight_damping /= 2.;
//      PCL_ERROR ("-> damping decreased to ");

      total_error_prev = total_error;
      PCL_ERROR ("Pose error accepted: %f\n", b.squaredNorm ());
    }
    PCL_ERROR ("damping weight %f\n", weight_damping);
    PCL_ERROR ("Errors: p2plane %f, reg_zero %f\n",
               b.squaredNorm (), weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ());


    lhs += weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                       3 * num_joints + 3);

    timer_sec.tic ();
    vars = lhs.ldlt ().solve (rhs);
    duration = timer_sec.toc ();
//    PCL_ERROR ("###TIME### LDLT solve %f ms.\n", duration);

    /// Put back the angles and the transformation
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      skeleton_.nodes[n_i]->angles += vars.block<3, 1> (3 * n_i, 0);
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars.block<3, 1> (3 * num_joints, 0);
  }

  return (total_error_prev / static_cast<double> (corresps.size ()));
}


float
af::SkeletonMesh::IKWithPointToPoint (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                      const pcl::Correspondences &corresps,
                                      const float weight_reg_zero)
{
  int num_joints = skeleton_.nodes.size ();

  /// Initialize some matrices
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1.;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J (Eigen::MatrixXf::Zero (3 * corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b (Eigen::VectorXf::Zero (3 * corresps.size ()));

  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));

  for (size_t iter_i = 0; iter_i < 150; ++iter_i)
  {
    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    Eigen::VectorXf vars_abs (num_joints * 3 + 3);
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);
//    std::cerr << vars_abs.transpose () << std::endl;

    J.setZero ();
    b.setZero ();

    /// Set up the constraints
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

          J.block<3, 1> (3 * c_i, 3 * node_id + 0) += weight_joint * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J.block<3, 1> (3 * c_i, 3 * node_id + 1) += weight_joint * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J.block<3, 1> (3 * c_i, 3 * node_id + 2) += weight_joint * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b.block<3, 1> (3 * c_i, 0) += weight_joint * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J.block<3, 3> (3 * c_i, 3 * num_joints) += weight_joint * Eigen::Matrix3f::Identity ();
      }
    }

    /// Solve the system
    lhs = J.transpose () * J +
          weight_reg_zero * J_reg_zero;
    rhs = J.transpose () * (-b) +
          weight_reg_zero * J_reg_zero * (-vars_abs);
    float total_error = b.squaredNorm () + (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ();

    if ((total_error_prev - total_error) / total_error_prev < 1e-3)
    {
      weight_damping *= 10.;
//      PCL_ERROR ("-> damping increased to ");
      if (weight_damping > 1e8)
        break;

      /// Take back the angles and the transformation
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars.block<3, 1> (3 * num_joints, 0);
    }
    else
    {
//      PCL_ERROR ("Registration error at inner iteration %zu: %f\n", iter_i, b.squaredNorm ());
//      PCL_ERROR ("Reg zero error: %f\n", (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ());
//      PCL_ERROR ("Total error: %f\n", total_error);

      weight_damping /= 2.;
//      PCL_ERROR ("-> damping decreased to ");

      total_error_prev = total_error;
      PCL_ERROR ("Pose error accepted: %f\n", b.squaredNorm ());
    }
    PCL_ERROR ("damping weight %f\n", weight_damping);
    PCL_ERROR ("Errors: p2plane %f, reg_zero %f\n",
               b.squaredNorm (), weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ());


    lhs += weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                       3 * num_joints + 3);

    timer_sec.tic ();
    vars = lhs.ldlt ().solve (rhs);
    duration = timer_sec.toc ();
//    PCL_ERROR ("###TIME### LDLT solve %f ms.\n", duration);

    /// Put back the angles and the transformation
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      skeleton_.nodes[n_i]->angles += vars.block<3, 1> (3 * n_i, 0);
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars.block<3, 1> (3 * num_joints, 0);
  }

  return (total_error_prev / static_cast<double> (corresps.size ()));
}

void
af::SkeletonMesh::initIndices ()
{
  /// Connect the vertices to the joints
  vertex_to_joints_weights_ = std::vector<std::vector<std::pair<int, float> > > (mesh_rest_->vertices_.cols ());

  for (size_t n_i = 0; n_i < skeleton_rest_.nodes.size (); ++n_i)
  {
    for (size_t v_i = 0; v_i < skeleton_rest_.nodes[n_i]->vertex_indices.size (); ++v_i)
      vertex_to_joints_weights_[skeleton_rest_.nodes[n_i]->vertex_indices[v_i]].push_back (std::make_pair (n_i,
                                                                                                           skeleton_rest_.nodes[n_i]->vertex_weights[v_i]));
  }
}


/*
float
af::SkeletonMesh::IKCustom (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                            pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                            const pcl::Correspondences &corresps,
                            const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                            const std::vector<Eigen::Vector3f> &nite_joint_pos,
                            const Eigen::VectorXf &vars_frame_prev,
                            const af::PosePCA &pose_pca,
                            Eigen::VectorXf &vars_result,
                            const float weight_p2plane,
                            const float weight_contour,
                            const float weight_reg_zero,
                            const float weight_close_prev,
                            const float weight_nite,
                            const float weight_pca_proj,
                            const float weight_pca_dev)
{
  /// Initialize some matrices
  int num_joints = skeleton_.nodes.size ();
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1.;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J_p2plane (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_p2plane (Eigen::VectorXf::Zero (corresps.size ()));

  Eigen::MatrixXf J_nite (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_contour (Eigen::VectorXf::Zero (corresps_contour.size ()));

  Eigen::MatrixXf J_nite (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
  Eigen::VectorXf b_nite (Eigen::VectorXf::Zero (nite_joint_pos.size () * 3));

  Eigen::MatrixXf MMt = pose_pca.modes * pose_pca.modes.transpose ();
  Eigen::MatrixXf J_pca_proj (Eigen::MatrixXf::Zero (3 * (num_joints - 1), 3 * num_joints + 3));
  J_pca_proj.block (0, 3, 3 * (num_joints - 1), 3 * (num_joints - 1)) = Eigen::MatrixXf::Identity (3 * (num_joints - 1), 3 * (num_joints - 1)) - MMt;
  Eigen::VectorXf b_pca_proj (Eigen::VectorXf::Zero (3 * (num_joints - 1)));

  Eigen::MatrixXf J_pca_dev (Eigen::MatrixXf::Zero (pose_pca.std_devs.rows (), 3 * num_joints + 3));
  J_pca_dev.block (0, 3, pose_pca.std_devs.rows (), 3 * (num_joints - 1)) = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose ();
  Eigen::VectorXf b_pca_dev (Eigen::VectorXf::Zero (pose_pca.std_devs.rows ()));


  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));

  Eigen::VectorXf vars_abs (Eigen::VectorXf::Zero (num_joints * 3 + 3));

  for (size_t iter_i = 0; iter_i < 150; ++iter_i)
  {
    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    J_p2plane.setZero ();
    b_p2plane.setZero ();
    J_nite.setZero ();
    b_contour.setZero ();

    timer_sec.tic ();
    /// 1. Set up the point2plane constraints
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

          J_p2plane (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_p2plane (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_p2plane.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }
    duration = timer_sec.toc ();
    PCL_ERROR ("###TIME### Building the p2plane system: %f ms.\n", duration);

    timer_sec.tic ();
    /// 2. Set up the contour constraints
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

          J_nite (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_nite (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_nite (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_contour (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_nite.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }
    duration = timer_sec.toc ();
    PCL_ERROR ("###TIME### Building the contour system: %f ms.\n", duration);


    timer_sec.tic ();
    /// 3. Set up the constraints for the nite joint positions
    /// Set up the constraints
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        af::Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 0) += T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 1) += T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 2) += T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b_nite.block<3, 1> (3 * c_i, 0) = node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i];

      J_nite.block<3, 3> (3 * c_i, 3 * num_joints) += Eigen::Matrix3f::Identity ();
    }
    duration = timer_sec.toc ();
    PCL_ERROR ("###TIME### Building the nite system: %f ms.\n", duration);


    timer_sec.tic ();
    /// 4. PCA regularization for the joint angles
    b_pca_proj = (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean) -
                  MMt * (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean);
    b_pca_dev = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose () *
                (vars_abs.block (3, 0, 3 * (num_joints - 1), 1)  - pose_pca.mean);
    duration = timer_sec.toc ();
    PCL_ERROR ("###TIME### Building the pca system: %f ms.\n", duration);


    timer_sec.tic ();
    /// Solve the system
    lhs = weight_p2plane * J_p2plane.transpose () * J_p2plane +
          weight_contour * J_nite.transpose () * J_nite +
          weight_nite * J_nite.transpose () * J_nite +
          weight_reg_zero * J_reg_zero +
          weight_close_prev * J_reg_zero +
          weight_pca_proj * J_pca_proj.transpose () * J_pca_proj +
          weight_pca_dev * J_pca_dev.transpose () * J_pca_dev;
    rhs = weight_p2plane * J_p2plane.transpose () * (-b_p2plane) +
          weight_contour * J_nite.transpose () * (-b_contour) +
          weight_nite * J_nite.transpose () * (-b_nite) +
          weight_reg_zero * J_reg_zero * (-vars_abs) +
          weight_close_prev * J_reg_zero * (-vars_abs + vars_frame_prev) +
          weight_pca_proj * J_pca_proj.transpose () * (-b_pca_proj) +
          weight_pca_dev * J_pca_dev.transpose () * (-b_pca_dev);
    float total_error = weight_p2plane * b_p2plane.squaredNorm () +
                        weight_contour * b_contour.squaredNorm () +
                        weight_nite * b_nite.squaredNorm () +
                        weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm () +
                        weight_close_prev * (J_reg_zero * (vars_frame_prev - vars_abs)).squaredNorm () +
                        weight_pca_proj * b_pca_proj.squaredNorm () +
                        weight_pca_dev * b_pca_dev.squaredNorm ();
    duration = timer_sec.toc ();
    PCL_ERROR ("###TIME### Building the lhs and rhs: %f ms.\n", duration);


    if ((total_error_prev - total_error) / total_error_prev < 5e-2)
    {
      weight_damping *= 10.;
//      PCL_ERROR ("-> damping increased to ");

      /// Take back the angles and the transformation
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars.block<3, 1> (3 * num_joints, 0);

      if (weight_damping > 1e5)
        break;
    }
    else
    {
//      PCL_ERROR ("Registration error at inner iteration %zu: %f\n", iter_i, b.squaredNorm ());
//      PCL_ERROR ("Reg zero error: %f\n", (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ());
//      PCL_ERROR ("Total error: %f\n", total_error);

      weight_damping /= 2.;
//      PCL_ERROR ("-> damping decreased to ");

      total_error_prev = total_error;
//      PCL_ERROR ("Pose error accepted: %f\n", b.squaredNorm ());
    }
    PCL_ERROR ("damping weight %f\n", weight_damping);
    PCL_ERROR ("Errors: p2plane %f, contour %f, nite: %f, reg_zero %f, reg_pca_proj %f, reg_pca_dev %f\n",
               weight_p2plane * b_p2plane.squaredNorm (),
               weight_contour * b_contour.squaredNorm (),
               weight_nite * b_nite.squaredNorm (),
               weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm (),
               weight_pca_proj * b_pca_proj.squaredNorm (),
               weight_pca_dev * b_pca_dev.squaredNorm ());
    PCL_ERROR ("-- Total error: %f\n", total_error);


    lhs += weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                       3 * num_joints + 3);

    timer_sec.tic ();
    vars = lhs.ldlt ().solve (rhs);
    duration = timer_sec.toc ();
//    PCL_ERROR ("###TIME### LDLT solve %f ms.\n", duration);

    /// Put back the angles and the transformation
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    {
      skeleton_.nodes[n_i]->angles += vars.block<3, 1> (3 * n_i, 0);
      af::clampAngles (skeleton_.nodes[n_i]->angles (0));
      af::clampAngles (skeleton_.nodes[n_i]->angles (1));
      af::clampAngles (skeleton_.nodes[n_i]->angles (2));
    }
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars.block<3, 1> (3 * num_joints, 0);
  }

  vars_result = vars_abs;
  return (total_error_prev / static_cast<double> (corresps.size ()));
}
*/



int total_ikcustom_iterations = 0;
int count_ikcustom_calls = 0;
float
af::SkeletonMesh::IKCustom (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                            pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                            const pcl::Correspondences &corresps,
                            const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                            const std::vector<Eigen::Vector3f> &nite_joint_pos,
                            const Eigen::VectorXf &vars_frame_prev,
                            const af::PosePCA &pose_pca,
                            Eigen::VectorXf &vars_result,
                            const float weight_p2plane,
                            const float weight_contour,
                            const float weight_reg_zero,
                            const float weight_close_prev,
                            const float weight_nite,
                            const float weight_pca_proj,
                            const float weight_pca_dev)
{
  /// Initialize some matrices
  int num_joints = skeleton_.nodes.size ();
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1e-12;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J_reg_close_prev (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));

  Eigen::MatrixXf J_p2plane (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_p2plane (Eigen::VectorXf::Zero (corresps.size ()));

  Eigen::MatrixXf J_contour (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_contour (Eigen::VectorXf::Zero (corresps_contour.size ()));

  Eigen::MatrixXf J_nite (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
  Eigen::VectorXf b_nite (Eigen::VectorXf::Zero (nite_joint_pos.size () * 3));

  Eigen::MatrixXf MMt = pose_pca.modes * pose_pca.modes.transpose ();
  Eigen::MatrixXf J_pca_proj (Eigen::MatrixXf::Zero (3 * (num_joints - 1), 3 * num_joints + 3));
  J_pca_proj.block (0, 3, 3 * (num_joints - 1), 3 * (num_joints - 1)) = Eigen::MatrixXf::Identity (3 * (num_joints - 1), 3 * (num_joints - 1)) - MMt;
  Eigen::VectorXf b_pca_proj (Eigen::VectorXf::Zero (3 * (num_joints - 1)));

  Eigen::MatrixXf J_pca_dev (Eigen::MatrixXf::Zero (pose_pca.std_devs.rows (), 3 * num_joints + 3));
  J_pca_dev.block (0, 3, pose_pca.std_devs.rows (), 3 * (num_joints - 1)) = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose ();
  Eigen::VectorXf b_pca_dev (Eigen::VectorXf::Zero (pose_pca.std_devs.rows ()));


  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));
  Eigen::VectorXf vars_inc (Eigen::VectorXf::Zero (3 * num_joints + 3));
  Eigen::VectorXf vars_abs (Eigen::VectorXf::Zero (3 * num_joints + 3));


  double error_p2plane, error_contour, error_reg_zero, error_close_prev, error_features, error_pca_proj, error_pca_dev;
  double total_error_normalized;
  size_t iter_i = 0;
  for (; iter_i < 150; ++iter_i)
  {
    error_p2plane = error_contour = error_reg_zero = error_close_prev = error_features = error_pca_proj = error_pca_dev = 0.;
    total_error_normalized = 0.;
    vars_inc.setZero ();
    if (iter_i != 0)
      vars_inc = (lhs + weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                                    3 * num_joints + 3)).ldlt ().solve (rhs);

    /// Update all the nodes temporarily
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      skeleton_.nodes[n_i]->angles += vars_inc.block<3, 1> (3 * n_i, 0);
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars_inc.block<3, 1> (3 * num_joints, 0);
    skeleton_.computeGlobalTransformations ();
    bakeGeometry ();
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    /// Compute the error
    double total_error = 0.;
    /// 1. Point2plane error
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_p2plane += weight_p2plane * val * val;
    }
    total_error += error_p2plane;

    /// 2. Contour error
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_contour += weight_contour * val * val;
    }
    total_error += error_contour;

    /// 3. Features Error
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      error_features += weight_nite * (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]).squaredNorm ();
    }
    total_error += error_features;

    /// 4. PCA regularization error
    error_pca_proj = weight_pca_proj * ((vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean) -
                                        MMt * (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean)).squaredNorm ();
    total_error += error_pca_proj;

    error_pca_dev = weight_pca_dev * (pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose () *
                                      (vars_abs.block (3, 0, 3 * (num_joints - 1), 1)  - pose_pca.mean)).squaredNorm ();
    total_error += error_pca_dev;

    /// 5. Reg zero and 6. Smoothness errors
    error_reg_zero = weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ();
    error_close_prev = weight_close_prev * (J_reg_close_prev * (-vars_abs + vars_frame_prev)).squaredNorm ();
    total_error += error_reg_zero + error_close_prev;

    /// Check if the error is good
    if ((total_error_prev - total_error) / total_error_prev < 1e-4) //5e-3) //5e-4) //5e-3) //1e-5)
    {
      /// Error bad

      /// Take out the last incremental updates
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars_inc.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars_inc.block<3, 1> (3 * num_joints, 0);

      /// Increase the damping
      weight_damping *= 100.;

      /// Stop if we are already damping too much
      if (weight_damping > 1e3)
        break;

      /// Go back to solving the problem with a different damping weight
      continue;
    }
    else
    {
      /// Error good
      weight_damping /= 5.;
      if (weight_damping < 1e-12)
        weight_damping = 1e-12;
      total_error_prev = total_error;

      total_error_normalized = error_p2plane / static_cast<double> (corresps.size ()) +
                               error_contour / static_cast<double> (corresps_contour.size ()) +
                               error_features + error_reg_zero + error_close_prev + error_pca_proj + error_pca_dev;
      /// Continue with computing new jacobians
    }

    PCL_ERROR ("Iteration: %zu - error %f,   damping %f\n", iter_i, total_error, weight_damping);
    PCL_ERROR ("   Errors: p2plane %f, contour %f, features %f, reg_zero %f, close_prev %f, pca_proj %f, pca_dev %f.\n",
               error_p2plane, error_contour, error_features, error_reg_zero, error_close_prev, error_pca_proj, error_pca_dev);



    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
//    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    J_p2plane.setZero ();
    b_p2plane.setZero ();
    J_contour.setZero ();
    b_contour.setZero ();
    J_nite.setZero ();
    b_nite.setZero ();

    /// 1. Set up the point2plane constraints
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_p2plane (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_p2plane (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_p2plane.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_p2plane_numdiff (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
    const double epsilon = 1e-3;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }

      for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
        J_p2plane_numdiff (c_i, vec_i) -= val / epsilon;
    }

    Eigen::VectorXf pose_current = skeleton_.getCurrentPose ();
    for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (pos_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
      {
        int v_i = corresps[c_i].index_query;
        mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
        mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
        const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
        const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

        /// Go through all of the joints that influence this vertex
        double val = 0.;
        for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
        {
          int j_index = vertex_to_joints_weights_[v_i][j_i].first;
          float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
          Skeleton::Node* node_end = skeleton_.nodes [j_index];
          transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

          transf_diff = node_end->global_transformation * transform_joint_rest_inv;
          val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
        }
        J_p2plane_numdiff (c_i, pos_i) += val / epsilon;
      }
    }
    PCL_ERROR ("DEBUG Jacobian p2plane diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_p2plane - J_p2plane_numdiff).squaredNorm (), (J_p2plane - J_p2plane_numdiff).squaredNorm () / J_p2plane_numdiff.squaredNorm (),
               J_p2plane.squaredNorm (), J_p2plane_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();
#endif // DEBUG_JACOBIANS

    /// 2. Set up the contour constraints
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_contour (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_contour (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_contour.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_contour_numdiff (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }

      for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
        J_contour_numdiff (c_i, vec_i) -= val / epsilon;
    }

    pose_current = skeleton_.getCurrentPose ();
    for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (vec_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
      {
        int v_i = corresps_contour[c_i].index_query;
        mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
        mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
        const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
        const Eigen::Vector3f &normal = normals_contour[c_i];

        /// Go through all of the joints that influence this vertex
        double val = 0.;
        for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
        {
          int j_index = vertex_to_joints_weights_[v_i][j_i].first;
          float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
          Skeleton::Node* node_end = skeleton_.nodes [j_index];
          transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

          transf_diff = node_end->global_transformation * transform_joint_rest_inv;
          val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
        }
        J_contour_numdiff (c_i, vec_i) += val / epsilon;
      }
    }
    PCL_ERROR ("DEBUG Jacobian contour diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_contour - J_contour_numdiff).squaredNorm (), (J_contour - J_contour_numdiff).squaredNorm () / J_contour_numdiff.squaredNorm (),
               J_contour.squaredNorm (), J_contour_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();
#endif // DEBUG_JACOBIANS

    /// 3. Set up the constraints for the nite joint positions
    /// Set up the constraints
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        af::Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 0) += T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 1) += T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 2) += T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b_nite.block<3, 1> (3 * c_i, 0) += node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i];

      J_nite.block<3, 3> (3 * c_i, 3 * num_joints) += Eigen::Matrix3f::Identity ();
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_nite_numdiff (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
        J_nite_numdiff.block<3, 1> (3 * c_i, pos_i) -= (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]) / epsilon;
    }


    pose_current = skeleton_.getCurrentPose ();
    for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (pos_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
      {
        /// Get the end effector for this constraint
        af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
        Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
        J_nite_numdiff.block<3, 1> (3 * c_i, pos_i) += (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]) / epsilon;
      }
    }

    PCL_ERROR ("DEBUG Jacobian nite diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_nite - J_nite_numdiff).squaredNorm (), (J_nite - J_nite_numdiff).squaredNorm () / J_nite_numdiff.squaredNorm (),
               J_nite.squaredNorm (), J_nite_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();

#endif // DEBUG_JACOBIANS

    /// 4. PCA regularization for the joint angles
    b_pca_proj = (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean) -
                  MMt * (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean);
    b_pca_dev = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose () *
                (vars_abs.block (3, 0, 3 * (num_joints - 1), 1)  - pose_pca.mean);

    /// Solve the system
    lhs = weight_p2plane * J_p2plane.transpose () * J_p2plane +
          weight_contour * J_contour.transpose () * J_contour +
          weight_nite * J_nite.transpose () * J_nite +
          weight_reg_zero * J_reg_zero +
          weight_close_prev * J_reg_close_prev +
          weight_pca_proj * J_pca_proj.transpose () * J_pca_proj +
          weight_pca_dev * J_pca_dev.transpose () * J_pca_dev;
    rhs = weight_p2plane * J_p2plane.transpose () * (-b_p2plane) +
          weight_contour * J_contour.transpose () * (-b_contour) +
          weight_nite * J_nite.transpose () * (-b_nite) +
          weight_reg_zero * J_reg_zero * (-vars_abs) +
          weight_close_prev * J_reg_close_prev * (-vars_abs + vars_frame_prev) +
          weight_pca_proj * J_pca_proj.transpose () * (-b_pca_proj) +
          weight_pca_dev * J_pca_dev.transpose () * (-b_pca_dev);
//    float total_error_new = weight_p2plane * b_p2plane.squaredNorm () +
//                            weight_contour * b_contour.squaredNorm () +
//                            weight_nite * b_nite.squaredNorm () +
//                            weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm () +
//                            weight_close_prev * (J_reg_zero * (vars_frame_prev - vars_abs)).squaredNorm () +
//                            weight_pca_proj * b_pca_proj.squaredNorm () +
//                            weight_pca_dev * b_pca_dev.squaredNorm ();

//    PCL_ERROR ("Sanity check %f vs %f\n",
//               total_error_prev, total_error_new);

    /// TODO clamp angles to nice values - see the old optimization
  }

  vars_result = Eigen::VectorXf (num_joints * 3 + 3);
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    vars_result.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
  vars_result.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

  total_ikcustom_iterations += iter_i;
  count_ikcustom_calls ++;
  PCL_ERROR ("AVG IKCUSTOM ITERATIONS: %f\n", static_cast<double> (total_ikcustom_iterations) / static_cast<double> (count_ikcustom_calls));

  return (total_error_prev / static_cast<double> (corresps.size ()));
}


float
af::SkeletonMesh::IKCustomFixed (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                 const pcl::Correspondences &corresps,
                                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                                 const std::vector<Eigen::Vector3f> &nite_joint_pos,
                                 const Eigen::VectorXf &vars_frame_prev,
                                 const af::PosePCA &pose_pca,
                                 Eigen::VectorXf &vars_result,
                                 const float weight_p2plane,
                                 const float weight_contour,
                                 const float weight_reg_zero,
                                 const float weight_close_prev,
                                 const float weight_nite,
                                 const float weight_pca_proj,
                                 const float weight_pca_dev)
{
  /// Initialize some matrices
  int num_joints = skeleton_.nodes.size ();
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1e-9;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J_reg_close_prev (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));

  Eigen::MatrixXf J_p2plane (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_p2plane (Eigen::VectorXf::Zero (corresps.size ()));

  Eigen::MatrixXf J_contour (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_contour (Eigen::VectorXf::Zero (corresps_contour.size ()));

  Eigen::MatrixXf J_nite (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
  Eigen::VectorXf b_nite (Eigen::VectorXf::Zero (nite_joint_pos.size () * 3));

  Eigen::MatrixXf MMt = pose_pca.modes * pose_pca.modes.transpose ();
  Eigen::MatrixXf J_pca_proj (Eigen::MatrixXf::Zero (3 * (num_joints - 1), 3 * num_joints + 3));
  J_pca_proj.block (0, 3, 3 * (num_joints - 1), 3 * (num_joints - 1)) = Eigen::MatrixXf::Identity (3 * (num_joints - 1), 3 * (num_joints - 1)) - MMt;
  Eigen::VectorXf b_pca_proj (Eigen::VectorXf::Zero (3 * (num_joints - 1)));

  Eigen::MatrixXf J_pca_dev (Eigen::MatrixXf::Zero (pose_pca.std_devs.rows (), 3 * num_joints + 3));
  J_pca_dev.block (0, 3, pose_pca.std_devs.rows (), 3 * (num_joints - 1)) = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose ();
  Eigen::VectorXf b_pca_dev (Eigen::VectorXf::Zero (pose_pca.std_devs.rows ()));


  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));
  Eigen::VectorXf vars_inc (Eigen::VectorXf::Zero (3 * num_joints + 3));


  double error_p2plane, error_contour, error_reg_zero, error_close_prev, error_features, error_pca_proj, error_pca_dev;
  double total_error_normalized;
  for (size_t iter_i = 0; iter_i < 150; ++iter_i)
  {
    error_p2plane = error_contour = error_reg_zero = error_close_prev = error_features = error_pca_proj = error_pca_dev = 0.;
    total_error_normalized = 0.;
    vars_inc.setZero ();
    if (iter_i != 0)
      vars_inc = (lhs + weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                                    3 * num_joints + 3)).ldlt ().solve (rhs);

    /// Update all the nodes temporarily
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      skeleton_.nodes[n_i]->angles += vars_inc.block<3, 1> (3 * n_i, 0);
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars_inc.block<3, 1> (3 * num_joints, 0);
    skeleton_.computeGlobalTransformations ();
    bakeGeometry ();
    Eigen::VectorXf vars_abs (num_joints * 3 + 3);
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    /// Compute the error
    double total_error = 0.;
    /// 1. Point2plane error
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_p2plane += weight_p2plane * val * val;
    }
    total_error += error_p2plane;

    /// 2. Contour error
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_contour += weight_contour * val * val;
    }
    total_error += error_contour;

    /// 3. Features Error
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      error_features += weight_nite * (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]).squaredNorm ();
    }
    total_error += error_features;

    /// 4. PCA regularization error
    error_pca_proj = weight_pca_proj * ((vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean) -
                                        MMt * (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean)).squaredNorm ();
    total_error += error_pca_proj;

    error_pca_dev = weight_pca_dev * (pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose () *
                                      (vars_abs.block (3, 0, 3 * (num_joints - 1), 1)  - pose_pca.mean)).squaredNorm ();
    total_error += error_pca_dev;

    /// 5. Reg zero and 6. Smoothness errors
    error_reg_zero = weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ();
    error_close_prev = weight_close_prev * (J_reg_close_prev * (-vars_abs + vars_frame_prev)).squaredNorm ();
    total_error += error_reg_zero + error_close_prev;

    /// Check if the error is good
    if ((total_error_prev - total_error) / total_error_prev < 1e-4)
    {
      /// Error bad

      /// Take out the last incremental updates
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars_inc.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars_inc.block<3, 1> (3 * num_joints, 0);

      /// Increase the damping
      weight_damping *= 10.;

      /// Stop if we are already damping too much
      if (weight_damping > 1e5)
        break;

      /// Go back to solving the problem with a different damping weight
      continue;
    }
    else
    {
      /// Error good
      weight_damping /= 2.;
      if (weight_damping < 1e-9)
        weight_damping = 1e-9;
      total_error_prev = total_error;

      total_error_normalized = error_p2plane / static_cast<double> (corresps.size ()) +
                               error_contour / static_cast<double> (corresps_contour.size ()) +
                               error_features + error_reg_zero + error_close_prev + error_pca_proj + error_pca_dev;
      /// Continue with computing new jacobians
    }

    PCL_ERROR ("Iteration: %zu - error %f,   damping %f\n", iter_i, total_error, weight_damping);
    PCL_ERROR ("   Errors: p2plane %f, contour %f, features %f, reg_zero %f, close_prev %f, pca_proj %f, pca_dev %f.\n",
               error_p2plane, error_contour, error_features, error_reg_zero, error_close_prev, error_pca_proj, error_pca_dev);



    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
//    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    J_p2plane.setZero ();
    b_p2plane.setZero ();
    J_contour.setZero ();
    b_contour.setZero ();
    J_nite.setZero ();
    b_nite.setZero ();

    /// 1. Set up the point2plane constraints
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_p2plane (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_p2plane (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_p2plane.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_p2plane_numdiff (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
    const double epsilon = 1e-3;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }

      for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
        J_p2plane_numdiff (c_i, vec_i) -= val / epsilon;
    }

    Eigen::VectorXf pose_current = skeleton_.getCurrentPose ();
    for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (pos_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
      {
        int v_i = corresps[c_i].index_query;
        mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
        mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
        const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
        const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

        /// Go through all of the joints that influence this vertex
        double val = 0.;
        for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
        {
          int j_index = vertex_to_joints_weights_[v_i][j_i].first;
          float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
          Skeleton::Node* node_end = skeleton_.nodes [j_index];
          transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

          transf_diff = node_end->global_transformation * transform_joint_rest_inv;
          val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
        }
        J_p2plane_numdiff (c_i, pos_i) += val / epsilon;
      }
    }
    PCL_ERROR ("DEBUG Jacobian p2plane diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_p2plane - J_p2plane_numdiff).squaredNorm (), (J_p2plane - J_p2plane_numdiff).squaredNorm () / J_p2plane_numdiff.squaredNorm (),
               J_p2plane.squaredNorm (), J_p2plane_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();
#endif // DEBUG_JACOBIANS

    /// 2. Set up the contour constraints
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_contour (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_contour (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_contour.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_contour_numdiff (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }

      for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
        J_contour_numdiff (c_i, vec_i) -= val / epsilon;
    }

    pose_current = skeleton_.getCurrentPose ();
    for (size_t vec_i = 0; vec_i < 3 * num_joints + 3; ++vec_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (vec_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
      {
        int v_i = corresps_contour[c_i].index_query;
        mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
        mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
        const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
        const Eigen::Vector3f &normal = normals_contour[c_i];

        /// Go through all of the joints that influence this vertex
        double val = 0.;
        for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
        {
          int j_index = vertex_to_joints_weights_[v_i][j_i].first;
          float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
          Skeleton::Node* node_end = skeleton_.nodes [j_index];
          transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

          transf_diff = node_end->global_transformation * transform_joint_rest_inv;
          val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
        }
        J_contour_numdiff (c_i, vec_i) += val / epsilon;
      }
    }
    PCL_ERROR ("DEBUG Jacobian contour diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_contour - J_contour_numdiff).squaredNorm (), (J_contour - J_contour_numdiff).squaredNorm () / J_contour_numdiff.squaredNorm (),
               J_contour.squaredNorm (), J_contour_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();
#endif // DEBUG_JACOBIANS

    /// 3. Set up the constraints for the nite joint positions
    /// Set up the constraints
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        af::Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 0) += T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 1) += T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 2) += T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b_nite.block<3, 1> (3 * c_i, 0) += node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i];

      J_nite.block<3, 3> (3 * c_i, 3 * num_joints) += Eigen::Matrix3f::Identity ();
    }

#ifdef DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_nite_numdiff (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
        J_nite_numdiff.block<3, 1> (3 * c_i, pos_i) -= (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]) / epsilon;
    }


    pose_current = skeleton_.getCurrentPose ();
    for (size_t pos_i = 0; pos_i < 3 * num_joints + 3; ++pos_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (pos_i) += epsilon;
      skeleton_.applyPose (pose_numdiff);
      bakeGeometry ();

      for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
      {
        /// Get the end effector for this constraint
        af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
        Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
        J_nite_numdiff.block<3, 1> (3 * c_i, pos_i) += (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]) / epsilon;
      }
    }

    PCL_ERROR ("DEBUG Jacobian nite diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J_nite - J_nite_numdiff).squaredNorm (), (J_nite - J_nite_numdiff).squaredNorm () / J_nite_numdiff.squaredNorm (),
               J_nite.squaredNorm (), J_nite_numdiff.squaredNorm ());

    skeleton_.applyPose (pose_current);
    bakeGeometry ();

#endif // DEBUG_JACOBIANS

    /// 4. PCA regularization for the joint angles
    b_pca_proj = (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean) -
                  MMt * (vars_abs.block (3, 0, 3 * (num_joints - 1), 1) - pose_pca.mean);
    b_pca_dev = pose_pca.std_devs.asDiagonal ().inverse () * pose_pca.modes.transpose () *
                (vars_abs.block (3, 0, 3 * (num_joints - 1), 1)  - pose_pca.mean);

    /// Solve the system
    lhs = weight_p2plane * J_p2plane.transpose () * J_p2plane +
          weight_contour * J_contour.transpose () * J_contour +
          weight_nite * J_nite.transpose () * J_nite +
          weight_reg_zero * J_reg_zero +
          weight_close_prev * J_reg_close_prev +
          weight_pca_proj * J_pca_proj.transpose () * J_pca_proj +
          weight_pca_dev * J_pca_dev.transpose () * J_pca_dev;
    rhs = weight_p2plane * J_p2plane.transpose () * (-b_p2plane) +
          weight_contour * J_contour.transpose () * (-b_contour) +
          weight_nite * J_nite.transpose () * (-b_nite) +
          weight_reg_zero * J_reg_zero * (-vars_abs) +
          weight_close_prev * J_reg_close_prev * (-vars_abs + vars_frame_prev) +
          weight_pca_proj * J_pca_proj.transpose () * (-b_pca_proj) +
          weight_pca_dev * J_pca_dev.transpose () * (-b_pca_dev);
//    float total_error_new = weight_p2plane * b_p2plane.squaredNorm () +
//                            weight_contour * b_contour.squaredNorm () +
//                            weight_nite * b_nite.squaredNorm () +
//                            weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm () +
//                            weight_close_prev * (J_reg_zero * (vars_frame_prev - vars_abs)).squaredNorm () +
//                            weight_pca_proj * b_pca_proj.squaredNorm () +
//                            weight_pca_dev * b_pca_dev.squaredNorm ();

//    PCL_ERROR ("Sanity check %f vs %f\n",
//               total_error_prev, total_error_new);

    /// TODO clamp angles to nice values - see the old optimization
  }

  vars_result = Eigen::VectorXf (num_joints * 3 + 3);
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    vars_result.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
  vars_result.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

  return (total_error_prev / static_cast<double> (corresps.size ()));
}






float
af::SkeletonMesh::IKCustomFelix (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                 const pcl::Correspondences &corresps,
                                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                                 const std::vector<Eigen::Vector3f> &nite_joint_pos,
                                 const Eigen::VectorXf &vars_frame_prev,
                                 Eigen::VectorXf &vars_result,
                                 const float weight_p2plane,
                                 const float weight_contour,
                                 const float weight_reg_zero,
                                 const float weight_close_prev,
                                 const float weight_nite)
{
  /// Initialize some matrices
  int num_joints = skeleton_.nodes.size ();
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0.f, 0.f, 0.f,
          0.f, 0.f, -1.f,
          0.f, 1.f, 0.f;
  Rydy << 0.f, 0.f, 1.f,
          0.f, 0.f, 0.f,
          -1.f, 0.f, 0.f;
  Rzdz << 0.f,-1.f, 0.f,
          1.f, 0.f, 0.f,
          0.f, 0.f, 0.f;

  pcl::console::TicToc timer_iteration, timer_sec;
  double duration = 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  float weight_damping = 1e-9;
  float total_error_prev = std::numeric_limits<float>::max ();

  /// Initialize storage for the precomputation of the rot_mat_dx
  std::vector<Eigen::Matrix3f> rot_mat_dx_vec (num_joints),
                               rot_mat_dy_vec (num_joints),
                               rot_mat_dz_vec (num_joints);

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf J_reg_close_prev (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));

  Eigen::MatrixXf J_p2plane (Eigen::MatrixXf::Zero (corresps.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_p2plane (Eigen::VectorXf::Zero (corresps.size ()));

  Eigen::MatrixXf J_contour (Eigen::MatrixXf::Zero (corresps_contour.size (), num_joints * 3 + 3));
  Eigen::VectorXf b_contour (Eigen::VectorXf::Zero (corresps_contour.size ()));

  Eigen::MatrixXf J_nite (Eigen::MatrixXf::Zero (nite_joint_pos.size () * 3, num_joints * 3 + 3));
  Eigen::VectorXf b_nite (Eigen::VectorXf::Zero (nite_joint_pos.size () * 3));

  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (num_joints * 3 + 3, num_joints * 3 + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (num_joints * 3 + 3));
  Eigen::VectorXf vars_inc (Eigen::VectorXf::Zero (3 * num_joints + 3));


  double error_p2plane, error_contour, error_reg_zero, error_close_prev, error_features, error_pca_proj, error_pca_dev;
  double total_error_normalized;
  for (size_t iter_i = 0; iter_i < 150; ++iter_i)
  {
    error_p2plane = error_contour = error_reg_zero = error_close_prev = error_features = error_pca_proj = error_pca_dev = 0.;
    total_error_normalized = 0.;
    vars_inc.setZero ();
    if (iter_i != 0)
      vars_inc = (lhs + weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                                    3 * num_joints + 3)).ldlt ().solve (rhs);

    /// Update all the nodes temporarily
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      skeleton_.nodes[n_i]->angles += vars_inc.block<3, 1> (3 * n_i, 0);
    skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars_inc.block<3, 1> (3 * num_joints, 0);
    skeleton_.computeGlobalTransformations ();
    bakeGeometry ();
    Eigen::VectorXf vars_abs (num_joints * 3 + 3);
    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    /// Compute the error
    double total_error = 0.;
    /// 1. Point2plane error
    Eigen::Vector4f mesh_point_rest, mesh_point;
    mesh_point_rest (3) = mesh_point (3) = 1.f;
    Eigen::Matrix4f transform_joint_rest_inv, transf_diff;
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_p2plane += weight_p2plane * val * val;
    }
    total_error += error_p2plane;

    /// 2. Contour error
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      double val = 0.;
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        val += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
      }
      error_contour += weight_contour * val * val;
    }
    total_error += error_contour;

    /// 3. Features Error
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      error_features += weight_nite * (node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i]).squaredNorm ();
    }
    total_error += error_features;

    /// 5. Reg zero and 6. Smoothness errors
    error_reg_zero = weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ();
    error_close_prev = weight_close_prev * (J_reg_close_prev * (-vars_abs + vars_frame_prev)).squaredNorm ();
    total_error += error_reg_zero + error_close_prev;

    /// Check if the error is good
    if ((total_error_prev - total_error) / total_error_prev < 1e-4)
    {
      /// Error bad

      /// Take out the last incremental updates
      for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
        skeleton_.nodes[n_i]->angles -= vars_inc.block<3, 1> (3 * n_i, 0);
      skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars_inc.block<3, 1> (3 * num_joints, 0);

      /// Increase the damping
      weight_damping *= 10.;

      /// Stop if we are already damping too much
      if (weight_damping > 1e5)
        break;

      /// Go back to solving the problem with a different damping weight
      continue;
    }
    else
    {
      /// Error good
      weight_damping /= 2.;
      if (weight_damping < 1e-9)
        weight_damping = 1e-9;
      total_error_prev = total_error;

      total_error_normalized = error_p2plane / static_cast<double> (corresps.size ()) +
                               error_contour / static_cast<double> (corresps_contour.size ()) +
                               error_features + error_reg_zero + error_close_prev + error_pca_proj + error_pca_dev;
      /// Continue with computing new jacobians
    }

    PCL_ERROR ("Iteration: %zu - error %f,   damping %f\n", iter_i, total_error, weight_damping);
    PCL_ERROR ("   Errors: p2plane %f, contour %f, features %f, reg_zero %f, close_prev %f, pca_proj %f, pca_dev %f.\n",
               error_p2plane, error_contour, error_features, error_reg_zero, error_close_prev, error_pca_proj, error_pca_dev);



    /// Precompute the rot_mat_d* for all the joints
    for (size_t n_i = 0; n_i < num_joints; ++n_i)
    {
      af::Skeleton::Node *node = skeleton_.nodes[n_i];
      Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                        Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                        Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().transpose ();
      rot_mat_dx_vec[n_i] = angles_rot_inv *
                            Rxdx * Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dy_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Rydy * Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
      rot_mat_dz_vec[n_i] = angles_rot_inv *
                            Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                            Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                            Rzdz * Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
    }

    duration = timer_iteration.toc ();
//    PCL_ERROR ("###TIME### Gauss-newton iteration took %f ms.\n", duration);
    timer_iteration.tic ();
    bakeGeometry ();

    for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    J_p2plane.setZero ();
    b_p2plane.setZero ();
    J_contour.setZero ();
    b_contour.setZero ();
    J_nite.setZero ();
    b_nite.setZero ();

    /// 1. Set up the point2plane constraints
    for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
    {
      int v_i = corresps[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_p2plane (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_p2plane (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_p2plane (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_p2plane.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

    /// 2. Set up the contour constraints
    for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
    {
      int v_i = corresps_contour[c_i].index_query;
      mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);
      mesh_point.block<3, 1> (0, 0) = mesh_->vertices_.col (v_i);
      const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
      const Eigen::Vector3f &normal = normals_contour[c_i];

      /// Go through all of the joints that influence this vertex
      for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
      {
        int j_index = vertex_to_joints_weights_[v_i][j_i].first;
        float weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
        Skeleton::Node* node_end = skeleton_.nodes [j_index];
        transform_joint_rest_inv = skeleton_rest_.getNode (node_end->name)->global_transformation.inverse ();

        /// Go through the kinematic chain up to this node and set the jacobians
        for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
        {
          Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

          int node_id = node_on_path->id;
          const Eigen::Matrix4f &T_prev = node_on_path->global_transformation;
          Eigen::Matrix4f T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation * transform_joint_rest_inv;

          J_contour (c_i, 3 * node_id + 0) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 1) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
          J_contour (c_i, 3 * node_id + 2) += weight_joint * normal.transpose () * T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * mesh_point_rest).block<3, 1> (0, 0);
        }

        transf_diff = node_end->global_transformation * transform_joint_rest_inv;
        b_contour (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);

        /// Add the jacobian for the translation of the body (HipCenter)
        J_contour.block<1, 3> (c_i, 3 * num_joints) += weight_joint * normal.transpose () * Eigen::Matrix3f::Identity ();
      }
    }

    /// 3. Set up the constraints for the nite joint positions
    /// Set up the constraints
    for (size_t c_i = 0; c_i < nite_joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      af::Skeleton::Node* node_end = skeleton_.getNode (skeleton_.map_num_to_name_ [c_i]);

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        af::Skeleton::Node *node_on_path = node_end->path_from_root[n_i];

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 0) += T_prev.block<3, 3> (0, 0) * rot_mat_dx_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 1) += T_prev.block<3, 3> (0, 0) * rot_mat_dy_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J_nite.block<3, 1> (3 * c_i, 3 * node_id + 2) += T_prev.block<3, 3> (0, 0) * rot_mat_dz_vec[node_id] * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b_nite.block<3, 1> (3 * c_i, 0) += node_end_position.block<3, 1> (0, 0) - nite_joint_pos[c_i];

      J_nite.block<3, 3> (3 * c_i, 3 * num_joints) += Eigen::Matrix3f::Identity ();
    }

    /// Solve the system
    lhs = weight_p2plane * J_p2plane.transpose () * J_p2plane +
          weight_contour * J_contour.transpose () * J_contour +
          weight_nite * J_nite.transpose () * J_nite +
          weight_reg_zero * J_reg_zero +
          weight_close_prev * J_reg_close_prev;
    rhs = weight_p2plane * J_p2plane.transpose () * (-b_p2plane) +
          weight_contour * J_contour.transpose () * (-b_contour) +
          weight_nite * J_nite.transpose () * (-b_nite) +
          weight_reg_zero * J_reg_zero * (-vars_abs) +
          weight_close_prev * J_reg_close_prev * (-vars_abs + vars_frame_prev);
  }

  vars_result = Eigen::VectorXf (num_joints * 3 + 3);
  for (size_t n_i = 0; n_i < skeleton_.nodes.size (); ++n_i)
    vars_result.block<3, 1> (3 * n_i, 0) = skeleton_.nodes[n_i]->angles;
  vars_result.block<3, 1> (3 * num_joints, 0) = skeleton_.getNode ("hip")->local_transformation.block<3, 1> (0, 3);

  return (total_error_prev / static_cast<double> (corresps.size ()));
}







////
///  MODELING
////
float
af::SkeletonMesh::accumulateModelingConstraints (pcl::PointCloud<pcl::PointXYZRGBA>::ConstPtr cloud,
                                                 pcl::PointCloud<pcl::Normal>::ConstPtr normals,
                                                 const pcl::Correspondences &corresps,
                                                 const pcl::Correspondences &corresps_contour, const std::vector<Eigen::Vector3f> &normals_contour,
                                                 af::BlendshapeModel<float>::ConstPtr bs_model,
                                                 Eigen::VectorXf &var_bs,
                                                 Eigen::MatrixXf &M_lhs, Eigen::VectorXf &M_rhs,
                                                 const int num_frames,
                                                 const float weight_p2plane, const float weight_contour)
{
  bakeGeometry ();

  Eigen::MatrixXf J_p2plane (Eigen::MatrixXf::Zero (corresps.size (), bs_model->num_exps_));
  Eigen::VectorXf b_p2plane (Eigen::VectorXf::Zero (corresps.size ()));
  Eigen::MatrixXf J_contour (Eigen::MatrixXf::Zero (corresps_contour.size (), bs_model->num_exps_));
  Eigen::VectorXf b_contour (Eigen::VectorXf::Zero (corresps_contour.size ()));

  /// Set up the constraints for the point2plane correspondences
  Eigen::Vector4f mesh_point_rest;
  mesh_point_rest (3) = 1.;
  Eigen::Matrix4f transf_diff;
  for (size_t c_i = 0; c_i < corresps.size (); ++c_i)
  {
    int v_i = corresps[c_i].index_query;
    const Eigen::Vector3f &point = cloud->at (corresps[c_i].index_match).getVector3fMap ();
    const Eigen::Vector3f &normal = normals->at (corresps[c_i].index_match).getNormalVector3fMap ();
    mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);

    /// Go through all of the joints that influence this vertex
    for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
    {
      int j_index = vertex_to_joints_weights_[v_i][j_i].first;
      double weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
      transf_diff = skeleton_.nodes[j_index]->global_transformation *
                    skeleton_rest_.nodes[j_index]->global_transformation.inverse ();

      J_p2plane.row (c_i) += weight_joint * normal.transpose () * transf_diff.block<3, 3> (0, 0) * bs_model->exps_.block (3 * v_i, 0, 3, bs_model->num_exps_);
      b_p2plane (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
    }
  }

  /// Set up the constraints for the contour correspondences
  for (size_t c_i = 0; c_i < corresps_contour.size (); ++c_i)
  {
    int v_i = corresps_contour[c_i].index_query;
    const Eigen::Vector3f &point = cloud->at (corresps_contour[c_i].index_match).getVector3fMap ();
    const Eigen::Vector3f &normal = normals_contour[c_i];
    mesh_point_rest.block<3, 1> (0, 0) = mesh_rest_->vertices_.col (v_i);

    /// Go through all of the joints that influence this vertex
    for (size_t j_i = 0; j_i < vertex_to_joints_weights_[v_i].size (); ++j_i)
    {
      int j_index = vertex_to_joints_weights_[v_i][j_i].first;
      double weight_joint = vertex_to_joints_weights_[v_i][j_i].second;
      transf_diff = skeleton_.nodes[j_index]->global_transformation *
                    skeleton_rest_.nodes[j_index]->global_transformation.inverse ();

      J_contour.row (c_i) += weight_joint * normal.transpose () * transf_diff.block<3, 3> (0, 0) * bs_model->exps_.block (3 * v_i, 0, 3, bs_model->num_exps_);
      b_contour (c_i) += weight_joint * normal.transpose () * (transf_diff.block<3, 4> (0, 0) * mesh_point_rest - point);
    }
  }

  const double num_frames_div = 1. / static_cast<double> (num_frames);
  /// Put the constraints into the big accumulated linear system
  M_lhs += num_frames_div * (weight_p2plane * J_p2plane.transpose () * J_p2plane +
                             weight_contour * J_contour.transpose () * J_contour);
//  M_rhs += num_frames_div * (weight_p2plane * J_p2plane.transpose () * (-b_p2plane) +
//                             weight_contour * J_contour.transpose () * (-b_contour));

  /// for GS - different linear system - in absolute values, not inc
  M_rhs += num_frames_div * (weight_p2plane * J_p2plane.transpose () * (J_p2plane * var_bs - b_p2plane) +
                             weight_contour * J_contour.transpose () * (J_contour * var_bs - b_contour));

  return (num_frames_div * (weight_p2plane * b_p2plane.squaredNorm () +
                            weight_contour * b_contour.squaredNorm ()));
}


float
af::SkeletonMesh::solveModelingProblem (const Eigen::MatrixXf &M_lhs, const Eigen::VectorXf &M_rhs,
                                        Eigen::VectorXf &var_bs,
                                        const float weight_bs_reg)
{
  /// Add the blendshape regularization
  /// Build up a new linear system - TODO check how inefficient this might be
  Eigen::MatrixXf A = M_lhs + weight_bs_reg * Eigen::MatrixXf::Identity (var_bs.rows (), var_bs.rows ());
  Eigen::VectorXf b = M_rhs + weight_bs_reg * (-var_bs);
  PCL_ERROR ("Blendshape regularization weight: %f\n", weight_bs_reg * var_bs.squaredNorm ());

  Eigen::VectorXf var_bs_inc = A.ldlt ().solve (b);

  var_bs += var_bs_inc;
}




float
af::SkeletonMesh::solveModelingProblemGS (const Eigen::MatrixXf &M_lhs, const Eigen::VectorXf &M_rhs,
                                          Eigen::VectorXf &var_bs,
                                          const float weight_bs_reg_l1)
{
  /// The Gauss-Seidel solving
  const int max_iter_gs = 1000;
  for (size_t iter_gs_i = 0; iter_gs_i < max_iter_gs; ++iter_gs_i)
  {
    for (size_t bs_i = 0; bs_i < var_bs.rows (); ++bs_i)
    {
      float v  = (M_lhs.row (bs_i).dot (var_bs) - M_lhs.coeff (bs_i, bs_i) * var_bs (bs_i) - M_rhs (bs_i));

      if(v > weight_bs_reg_l1)
        var_bs (bs_i) = (weight_bs_reg_l1 - v) / M_lhs.coeff (bs_i, bs_i);
      else if(v < -weight_bs_reg_l1)
        var_bs (bs_i) = (-weight_bs_reg_l1 - v) / M_lhs.coeff (bs_i, bs_i);
      else var_bs (bs_i) = 0.0f;

      /// Clamp it, clamp it!
      if (var_bs (bs_i) < 0.f)
        var_bs (bs_i) = 0.f;
      if (var_bs (bs_i) > 1.f)
        var_bs (bs_i) = 1.f;
    }
  }

  return 0.0f;
}



