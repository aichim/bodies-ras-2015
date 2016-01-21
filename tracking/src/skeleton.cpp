#include "skeleton.h"
#include "skeleton_mesh.h"

af::Skeleton::Node::Node (const af::Skeleton::Node &copy)
{
  parent = NULL;

  name = copy.name;
  id = copy.id;
  angles = copy.angles;

  local_transformation = copy.local_transformation;
  global_transformation = copy.global_transformation;
  angles = copy.angles;


  min_scaling = copy.min_scaling;
  max_scaling = copy.max_scaling;
  min_rotation = copy.min_rotation;
  max_rotation = copy.max_rotation;
  min_translation = copy.min_translation;
  max_translation = copy.max_translation;

  vertex_indices = copy.vertex_indices;
  vertex_weights = copy.vertex_weights;
}


af::Skeleton::Node::Node (std::string _name,
                          Eigen::Matrix4f &_local_transformation,
                          Eigen::Vector3f &_min_scaling,
                          Eigen::Vector3f &_max_scaling,
                          Eigen::Vector3f &_min_rotation,
                          Eigen::Vector3f &_max_rotation,
                          Eigen::Vector3f &_min_translation,
                          Eigen::Vector3f &_max_translation)
{
  parent = NULL;
  name = _name;

  local_transformation = _local_transformation;
  angles = Eigen::Vector3f::Zero ();
  global_transformation = Eigen::Matrix4f::Identity ();

  min_scaling = _min_scaling;
  max_scaling = _max_scaling;
  min_rotation = _min_rotation;
  max_rotation = _max_rotation;
  min_translation = _min_translation;
  max_translation = _max_translation;
}


Eigen::Matrix3f af::Skeleton::Rxdx, af::Skeleton::Rydy, af::Skeleton::Rzdz;
af::Skeleton::Skeleton ()
{
  bindJointNames ();

  Rxdx << 0., 0., 0.,
          0., 0., -1.,
          0., 1., 0.;
  Rydy << 0., 0., 1.,
          0., 0., 0.,
          -1., 0., 0.;
  Rzdz << 0.,-1., 0.,
          1., 0., 0.,
          0., 0., 0.;
}


af::Skeleton&
af::Skeleton::operator =(const af::Skeleton &copy)
{
  /// Copy the nodes without the pointers
  nodes.clear ();
  nodes.reserve (copy.nodes.size ());
  name_to_node.clear ();
  for (size_t i = 0; i < copy.nodes.size (); ++i)
  {
    Node *new_node = new Node (*copy.nodes[i]);
    nodes.push_back (new_node);
    name_to_node.insert( std::make_pair (new_node->name,new_node) );
  }

  /// Now add in the pointers
  for (size_t i = 0; i < copy.nodes.size (); ++i)
  {
    std::vector<Node*> children = copy.nodes[i]->getAllChildren ();

    for (size_t j = 0; j < children.size (); ++j)
    {
      nodes[i]->insertChild (getNode (children[j]->name));
      getNode (children[j]->name)->parent = nodes[i];
    }
  }

  return (*this);
}


af::Skeleton::Skeleton (const Skeleton &copy)
{
  /// Copy the nodes without the pointers
  nodes.reserve (copy.nodes.size ());
  for (size_t i = 0; i < copy.nodes.size (); ++i)
  {
    Node *new_node = new Node (*copy.nodes[i]);
    nodes.push_back (new_node);
    name_to_node.insert( std::make_pair (new_node->name,new_node) );
  }

  /// Now add in the pointers
  for (size_t i = 0; i < copy.nodes.size (); ++i)
  {
    std::vector<Node*> children = copy.nodes[i]->getAllChildren ();

    for (size_t j = 0; j < children.size (); ++j)
    {
      nodes[i]->insertChild (getNode (children[j]->name));
      getNode (children[j]->name)->parent = nodes[i];
    }
  }

  bindJointNames ();
}


void
af::Skeleton::insertNode (af::Skeleton::Node* new_node,
                          std::string parent_name)
{
    const std::string& name = new_node->getNodeName();
    if(name.empty()) assert(NODE_MUST_HAVE_A_NAME);
    // std::cout << "Inserting node '" << name.c_str() << "'" << std::endl;
    if(parent_name.empty() && nodes.size()!=0) assert(ROOT_MUST_BE_INSERTED_FIRST);
    bool success = name_to_node.insert( std::make_pair (name,new_node) ).second;
    if(!success) assert(ATTEMPTED_INSERT_DUPLICATE_NODE); /// No duplicates allowed

    /// attach the node
    nodes.push_back(new_node);
    new_node->id = nodes.size () - 1;

    /// Create parent-child relationship
    if(!parent_name.empty()){
        new_node->parent = getNode(parent_name);
        new_node->parent->insertChild(new_node);
    }
}


void
af::Skeleton::computeGlobalTransformations ()
{
  Node *node = getRoot ();

  Eigen::Matrix3f angles_rot = (Eigen::AngleAxisf (node->angles (0), Eigen::Vector3f::UnitX ()) *
                                Eigen::AngleAxisf (node->angles (1), Eigen::Vector3f::UnitY ()) *
                                Eigen::AngleAxisf (node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ();
  Eigen::Matrix4f aux_mat = Eigen::Matrix4f::Identity ();
  aux_mat.block<3, 3> (0, 0) = angles_rot;
  node->global_transformation = node->local_transformation * aux_mat;

  std::deque<Node*> queue;
  queue.push_back (node);

  /// BFS on the tree and accumulate the transformations downwards
  while (!queue.empty ())
  {
    node = queue.front ();
    queue.pop_front ();

    for (size_t i = 0; i < node->getAllChildren ().size (); ++i)
    {
      Node *child_node = node->getAllChildren ()[i];
      queue.push_back (child_node);

//      child_node->global_transformation = node->global_transformation * child_node->local_transformation;
      angles_rot = (Eigen::AngleAxisf (child_node->angles (0), Eigen::Vector3f::UnitX ()) *
                    Eigen::AngleAxisf (child_node->angles (1), Eigen::Vector3f::UnitY ()) *
                    Eigen::AngleAxisf (child_node->angles (2), Eigen::Vector3f::UnitZ ())).matrix ();
      aux_mat.block<3, 3> (0, 0) = angles_rot;
      child_node->global_transformation = node->global_transformation * child_node->local_transformation * aux_mat;
    }
  }
}


void
af::Skeleton::bindJointNames ()
{
//  map_num_to_name_.insert (std::make_pair (0, "HipCenter"));
//  map_num_to_name_.insert (std::make_pair (1, "Spine"));
//  map_num_to_name_.insert (std::make_pair (2, "ShoulderCenter"));
//  map_num_to_name_.insert (std::make_pair (3, "Head"));
//  map_num_to_name_.insert (std::make_pair (4, "ShoulderRight"));
//  map_num_to_name_.insert (std::make_pair (5, "ElbowRight"));
//  map_num_to_name_.insert (std::make_pair (6, "WristRight"));
//  map_num_to_name_.insert (std::make_pair (7, "HandRight"));
//  map_num_to_name_.insert (std::make_pair (8, "ShoulderLeft"));
//  map_num_to_name_.insert (std::make_pair (9, "ElbowLeft"));
//  map_num_to_name_.insert (std::make_pair (10, "WristLeft"));
//  map_num_to_name_.insert (std::make_pair (11, "HandLeft"));
//  map_num_to_name_.insert (std::make_pair (12, "HipRight"));
//  map_num_to_name_.insert (std::make_pair (13, "KneeRight"));
//  map_num_to_name_.insert (std::make_pair (14, "AnkleRight"));
//  map_num_to_name_.insert (std::make_pair (15, "FootRight"));
//  map_num_to_name_.insert (std::make_pair (16, "HipLeft"));
//  map_num_to_name_.insert (std::make_pair (17, "KneeLeft"));
//  map_num_to_name_.insert (std::make_pair (18, "AnkleLeft"));
//  map_num_to_name_.insert (std::make_pair (19, "FootLeft"));

/*
  map_num_to_name_.insert (std::make_pair (0, "Spine"));
  map_num_to_name_.insert (std::make_pair (1, "ShoulderCenter"));
  map_num_to_name_.insert (std::make_pair (8, "Head"));
  map_num_to_name_.insert (std::make_pair (2, "ShoulderLeft"));
  map_num_to_name_.insert (std::make_pair (4, "ElbowLeft"));
  map_num_to_name_.insert (std::make_pair (6, "WristLeft"));
  map_num_to_name_.insert (std::make_pair (3, "ShoulderRight"));
  map_num_to_name_.insert (std::make_pair (5, "ElbowRight"));
  map_num_to_name_.insert (std::make_pair (7, "WristRight"));
  map_num_to_name_.insert (std::make_pair (9, "HipLeft"));
  map_num_to_name_.insert (std::make_pair (11, "KneeLeft"));
  map_num_to_name_.insert (std::make_pair (13, "AnkleLeft"));
  map_num_to_name_.insert (std::make_pair (10, "HipRight"));
  map_num_to_name_.insert (std::make_pair (12, "KneeRight"));
  map_num_to_name_.insert (std::make_pair (14, "AnkleRight"));
*/

/*
  /// Naming convention for the datasets we collected
  /// USED FOR THE RAS PAPER
  map_num_to_name_.insert (std::make_pair (0, "chest"));
  map_num_to_name_.insert (std::make_pair (1, "neck"));
  map_num_to_name_.insert (std::make_pair (8, "head"));
  map_num_to_name_.insert (std::make_pair (2, "lShldr"));
  map_num_to_name_.insert (std::make_pair (4, "lForeArm"));
  map_num_to_name_.insert (std::make_pair (6, "lHand"));
  map_num_to_name_.insert (std::make_pair (3, "rShldr"));
  map_num_to_name_.insert (std::make_pair (5, "rForeArm"));
  map_num_to_name_.insert (std::make_pair (7, "rHand"));
  map_num_to_name_.insert (std::make_pair (9, "lThigh"));
  map_num_to_name_.insert (std::make_pair (11, "lShin"));
  map_num_to_name_.insert (std::make_pair (13, "lFoot"));
  map_num_to_name_.insert (std::make_pair (10, "rThigh"));
  map_num_to_name_.insert (std::make_pair (12, "rShin"));
  map_num_to_name_.insert (std::make_pair (14, "rFoot"));
*/

  /*
  /// For generating the training set from the TrailersPark BVH files
  map_num_to_name_.insert (std::make_pair (0, "chest"));
  map_num_to_name_.insert (std::make_pair (1, "lFoot"));
  map_num_to_name_.insert (std::make_pair (2, "rFoot"));
  map_num_to_name_.insert (std::make_pair (3, "lForeArm"));
  map_num_to_name_.insert (std::make_pair (4, "rForeArm"));
  map_num_to_name_.insert (std::make_pair (5, "lHand"));
  map_num_to_name_.insert (std::make_pair (6, "rHand"));
  map_num_to_name_.insert (std::make_pair (7, "head"));
  map_num_to_name_.insert (std::make_pair (8, "hip"));
  map_num_to_name_.insert (std::make_pair (9, "lShin"));
  map_num_to_name_.insert (std::make_pair (10, "rShin"));
  map_num_to_name_.insert (std::make_pair (11, "lThigh"));
  map_num_to_name_.insert (std::make_pair (12, "rThigh"));
  map_num_to_name_.insert (std::make_pair (13, "abdomen"));
*/

#if 1
  /// EVAL dataset mapping
  map_num_to_name_.insert (std::make_pair (0, "head"));
  map_num_to_name_.insert (std::make_pair (1, "neck"));
  map_num_to_name_.insert (std::make_pair (2, "rShldr"));
  map_num_to_name_.insert (std::make_pair (3, "rForeArm"));
  map_num_to_name_.insert (std::make_pair (4, "rHand"));
  map_num_to_name_.insert (std::make_pair (5, "lShldr"));
  map_num_to_name_.insert (std::make_pair (6, "lForeArm"));
  map_num_to_name_.insert (std::make_pair (7, "lHand"));
  map_num_to_name_.insert (std::make_pair (8, "rShin"));
  map_num_to_name_.insert (std::make_pair (9, "rFoot"));
  map_num_to_name_.insert (std::make_pair (10, "lShin"));
  map_num_to_name_.insert (std::make_pair (11, "lFoot"));
#else
  /// Felix's convention
  map_num_to_name_.insert (std::make_pair (0, "abdomen"));
  map_num_to_name_.insert (std::make_pair (1, "chest"));
  map_num_to_name_.insert (std::make_pair (2, "head"));
  map_num_to_name_.insert (std::make_pair (3, "hip"));
  map_num_to_name_.insert (std::make_pair (4, "lFoot"));
  map_num_to_name_.insert (std::make_pair (5, "lForeArm"));
  map_num_to_name_.insert (std::make_pair (6, "lHand"));
  map_num_to_name_.insert (std::make_pair (7, "lShin"));
  map_num_to_name_.insert (std::make_pair (8, "lShldr"));
  map_num_to_name_.insert (std::make_pair (9, "lThigh"));
  map_num_to_name_.insert (std::make_pair (10, "neck"));
  map_num_to_name_.insert (std::make_pair (11, "rFoot"));
  map_num_to_name_.insert (std::make_pair (12, "rForeArm"));
  map_num_to_name_.insert (std::make_pair (13, "rHand"));
  map_num_to_name_.insert (std::make_pair (14, "rShin"));
  map_num_to_name_.insert (std::make_pair (15, "rShldr"));
  map_num_to_name_.insert (std::make_pair (16, "rThigh"));
#endif

/*
  /// Naming convention for Matteo's datasets
  map_num_to_name_.insert (std::make_pair (0, "chest"));
  map_num_to_name_.insert (std::make_pair (1, "chest"));
  map_num_to_name_.insert (std::make_pair (2, "neck"));
  map_num_to_name_.insert (std::make_pair (3, "head"));
  map_num_to_name_.insert (std::make_pair (4, "rShldr"));
  map_num_to_name_.insert (std::make_pair (5, "rForeArm"));
  map_num_to_name_.insert (std::make_pair (6, "rHand"));
  map_num_to_name_.insert (std::make_pair (7, "rHand"));
  map_num_to_name_.insert (std::make_pair (8, "lShldr"));
  map_num_to_name_.insert (std::make_pair (9, "lForeArm"));
  map_num_to_name_.insert (std::make_pair (10, "lHand"));
  map_num_to_name_.insert (std::make_pair (11, "lHand"));
  map_num_to_name_.insert (std::make_pair (12, "rThigh"));
  map_num_to_name_.insert (std::make_pair (13, "rShin"));
  map_num_to_name_.insert (std::make_pair (14, "rFoot"));
  map_num_to_name_.insert (std::make_pair (15, "rFoot"));
  map_num_to_name_.insert (std::make_pair (16, "lThigh"));
  map_num_to_name_.insert (std::make_pair (17, "lShin"));
  map_num_to_name_.insert (std::make_pair (18, "lFoot"));
  map_num_to_name_.insert (std::make_pair (19, "lFoot"));
*/


  for (std::map<int, std::string>::const_iterator m_it = map_num_to_name_.begin ();
       m_it != map_num_to_name_.end (); ++m_it)
    map_name_to_num_.insert (std::make_pair (m_it->second, m_it->first));
}


void
af::Skeleton::clearPaths ()
{
  for (size_t i = 0; i < nodes.size (); ++i)
    nodes[i]->path_from_root.clear ();
}


void
af::Skeleton::computePathsDFS (Node* node)
{
  node->path_from_root.push_back (node);
  for (size_t c_i = 0; c_i < node->children.size (); ++c_i)
  {
    Node *child_node = node->children[c_i];
    child_node->path_from_root = node->path_from_root;

    computePathsDFS (child_node);
  }
}



Eigen::Matrix4f
af::Skeleton::jacrot (const Eigen::Vector4d &p)
{
  Eigen::Matrix4f mat;
  mat <<  0., p (2), -p (1), 0.,
      -p (2), 0., p (0), 0.,
      p (1), -p (0), 0., 0.,
      0., 0., 0., 1.;
  return (mat);
}


#if 0
void
af::Skeleton::IKWithJointPositions (std::vector<Eigen::Vector3f> &joint_pos)
{
  int num_joints = nodes.size ();

  /// TODO add regularization to keep the angles as close to the rest pose as possible

  std::vector<int> indices = {0, 9, 10};
  double error = 0.;
  for (size_t iter_i = 0; iter_i < 1000; ++iter_i)
  {
    computeGlobalTransformations ();

    Eigen::MatrixXf J (Eigen::MatrixXf::Zero (joint_pos.size () * 3, num_joints * 3 + 3));
    Eigen::VectorXf b (Eigen::VectorXf::Zero (joint_pos.size () * 3));

    /// Set up the constraints
    for (size_t i_i = 0; i_i < indices.size (); ++i_i)
    {
      size_t c_i = indices[i_i];
      /// Get the end effector for this constraint
      Node* node_end = getNode (map_num_to_name_ [c_i]);
//      std::cerr << "\nSetting up the jacobian for node " << node_end->name << " : ";

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        Node *node_on_path = node_end->path_from_root[n_i];
//        std::cerr << node_on_path->name << " -> ";

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        Eigen::Matrix4f jac = T_prev * jacrot (T_next * Eigen::Vector4f (0., 0., 0., 1.));
        J.block<3, 3> (3 * c_i, 3 * node_id) = jac.block<3, 3> (0, 0);
      }

      Eigen::Vector4d node_end_position = node_end->global_transformation * Eigen::Vector4d (0., 0., 0., 1.);
      b.block<3, 1> (3 * c_i, 0) = node_end_position.block<3, 1> (0, 0) - joint_pos[c_i];

      J.block<3, 3> (3 * c_i, 3 * num_joints) = Eigen::Matrix3f::Identity ();
    }
    error = b.squaredNorm ();
    if (iter_i == 0)
      PCL_ERROR ("IK with Joint positions optimization, initial error: %f ", error);

    PCL_ERROR ("IK error - iteration %zu -> %f\n", iter_i, b.squaredNorm ());

    /// Solve the system
    const double weight_damping = 1.;
    Eigen::VectorXf angles_inc = (J.transpose () * J +
                                  weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                                              3 * num_joints + 3)).ldlt ().solve (-J.transpose () * b);
//    std::cerr << "angles_inc: " << angles_inc.transpose () << std::endl;

    /// Put back the incremental angles
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
    {
      int node_id = n_i; //map_name_to_num_[nodes[n_i]->name];
      Eigen::Matrix3f rot_inc = (Eigen::AngleAxisf (angles_inc (3 * node_id + 0), Eigen::Vector3f::UnitX ()) *
                                 Eigen::AngleAxisf (angles_inc (3 * node_id + 1), Eigen::Vector3f::UnitY ()) *
                                 Eigen::AngleAxisf (angles_inc (3 * node_id + 2), Eigen::Vector3f::UnitZ ())).matrix ();

      nodes[n_i]->local_transformation.block<3, 3> (0, 0) = nodes[n_i]->local_transformation.block<3, 3> (0, 0) * rot_inc;
    }
    getNode ("hip")->local_transformation.block<3, 1> (0, 3) += angles_inc.block<3, 1> (3 * num_joints, 0);
  }

  PCL_ERROR (", final: %f\n", error);
}

#else
/*
void
af::Skeleton::IKWithJointPositions (const std::vector<Eigen::Vector3f> &joint_pos)
{
  int num_joints = nodes.size ();

  /// Initialize some matrices
  Eigen::Matrix3f Rxdx, Rydy, Rzdz;
  Rxdx << 0., 0., 0.,
          0., 0., -1.,
          0., 1., 0.;
  Rydy << 0., 0., 1.,
          0., 0., 0.,
          -1., 0., 0.;
  Rzdz << 0.,-1., 0.,
          1., 0., 0.,
          0., 0., 0.;

  Eigen::VectorXf vars (Eigen::VectorXf::Zero (3 * num_joints + 3));
  double weight_damping = 1e-3;
  double total_error_prev = std::numeric_limits<double>::max ();

  /// Initialize the variables
  for (size_t iter_i = 0; iter_i < 1000; ++iter_i)
  {
    Eigen::VectorXf vars_abs (num_joints * 3 + 3);
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = getNode ("hip")->local_transformation.block<3, 1> (0, 3);

//    std::cerr << "vars abs: " << vars_abs.transpose () << std::endl;

    computeGlobalTransformations ();

    Eigen::MatrixXf J (Eigen::MatrixXf::Zero (joint_pos.size () * 3, num_joints * 3 + 3));
    Eigen::VectorXf b (Eigen::VectorXf::Zero (joint_pos.size () * 3));

    /// Set up the constraints
    for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      Node* node_end = getNode (map_num_to_name_ [c_i]);
//      std::cerr << "\nSetting up the jacobian for node " << node_end->name << " : ";

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        Node *node_on_path = node_end->path_from_root[n_i];
//        std::cerr << node_on_path->name << " -> ";

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4f T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()) *
                                          Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()) *
                                          Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().inverse ();
        Eigen::Matrix3f rot_mat_dx = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () * Rxdx *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
        Eigen::Matrix3f rot_mat_dy = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () * Rydy *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
        Eigen::Matrix3f rot_mat_dz = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix () * Rzdz;

        J.block<3, 1> (3 * c_i, 3 * node_id + 0) = T_prev.block<3, 3> (0, 0) * rot_mat_dx * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 1) = T_prev.block<3, 3> (0, 0) * rot_mat_dy * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 2) = T_prev.block<3, 3> (0, 0) * rot_mat_dz * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b.block<3, 1> (3 * c_i, 0) = node_end_position.block<3, 1> (0, 0) - joint_pos[c_i];

      J.block<3, 3> (3 * c_i, 3 * num_joints) = Eigen::Matrix3f::Identity ();
    }

//    PCL_ERROR ("IK error - iteration %zu -> %f\n", iter_i, b.squaredNorm ());

    /// Solve the system
    const double weight_reg_zero = 1e-1; //1e-5;
    Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
    J_reg_zero.block<3, 3> (0, 0).setZero ();
    J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();
    Eigen::MatrixXf lhs = J.transpose () * J +
                          weight_reg_zero * J_reg_zero;
    Eigen::VectorXf rhs = J.transpose () * (-b) +
                          weight_reg_zero * J_reg_zero * (-vars_abs);
    double total_error = b.squaredNorm () + (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ();

    if ((total_error_prev - total_error) / total_error_prev < 1e-5)
    {
      weight_damping *= 10.;
//      PCL_ERROR ("-> damping increased to ");
      if (weight_damping > 1e4)
        break;

      /// Take back the angles and the transformation
      for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
        nodes[n_i]->angles -= vars.block<3, 1> (3 * n_i, 0);
      getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars.block<3, 1> (3 * num_joints, 0);
    }
    else
    {
      weight_damping /= 2.;
//      PCL_ERROR ("-> damping decreased to ");

      total_error_prev = total_error;
    }
//    PCL_ERROR ("%f\n", weight_damping);

//    PCL_ERROR ("Reg zero error: %f\n", (weight_reg_zero * J_reg_zero * vars_abs).squaredNorm ());
//    PCL_ERROR ("Total error: %f\n", total_error);


    lhs += weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                       3 * num_joints + 3);
    vars = lhs.ldlt ().solve (rhs);

    /// Put back the angles and the transformation
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
      nodes[n_i]->angles += vars.block<3, 1> (3 * n_i, 0);
    getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars.block<3, 1> (3 * num_joints, 0);
  }
}*/



int total_ikjoint_iterations = 0;
int count_ikjoint_calls = 0;
void
af::Skeleton::IKWithJointPositions (const std::vector<Eigen::Vector3f> &joint_pos,
                                    const double weight_reg_zero)
{
  int num_joints = nodes.size ();
  PCL_ERROR ("IKWITHJOINTSPOS --- num_joints %d\n", num_joints);

  double weight_damping = 1e3;
  double total_error_prev = std::numeric_limits<double>::max ();

  Eigen::MatrixXf J_reg_zero (Eigen::MatrixXf::Identity (3 * num_joints + 3, 3 * num_joints + 3));
  J_reg_zero.block<3, 3> (0, 0).setZero ();
  J_reg_zero.block<3, 3> (3 * num_joints, 3 * num_joints).setZero ();

  Eigen::MatrixXf lhs (Eigen::MatrixXf::Zero (3 * num_joints + 3, 3 * num_joints + 3));
  Eigen::VectorXf rhs (Eigen::VectorXf::Zero (3 * num_joints + 3));
  Eigen::VectorXf vars_inc (Eigen::VectorXf::Zero (3 * num_joints + 3));
  Eigen::VectorXf vars_abs (Eigen::VectorXf::Zero (3 * num_joints + 3));
  Eigen::MatrixXf J (Eigen::MatrixXf::Zero (joint_pos.size () * 3, num_joints * 3 + 3));
  Eigen::VectorXf b (Eigen::VectorXf::Zero (joint_pos.size () * 3));
  Eigen::Matrix4f T_prev, T_next;

  double total_joints_error = 0.;
  /// Initialize the variables
  size_t iter_i = 0;
  for (; iter_i < 1000; ++iter_i)
  {
    /// Compute the possible variable increment
    vars_inc.setZero ();
    if (iter_i != 0)
      vars_inc = (lhs + weight_damping * Eigen::MatrixXf::Identity (3 * num_joints + 3,
                                                                    3 * num_joints + 3)).ldlt ().solve (rhs);

    /// Update all the nodes temporarily
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
      nodes[n_i]->angles += vars_inc.block<3, 1> (3 * n_i, 0);
    getNode ("hip")->local_transformation.block<3, 1> (0, 3) += vars_inc.block<3, 1> (3 * num_joints, 0);
    computeGlobalTransformations ();

    /// Compute the error
    double total_error = 0.;
    for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
    {
      Node* node_end = getNode (map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      total_error += (node_end_position.block<3, 1> (0, 0) - joint_pos[c_i]).squaredNorm ();
    }

    total_joints_error = total_error;
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
      vars_abs.block<3, 1> (3 * n_i, 0) = nodes[n_i]->angles;
    vars_abs.block<3, 1> (3 * num_joints, 0) = getNode ("hip")->local_transformation.block<3, 1> (0, 3);

    if (iter_i > 15)
      total_error += weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ();
    else
      total_error += 1e3 * (J_reg_zero * vars_abs).squaredNorm ();

    if (std::isnan (total_error))
    {
      PCL_ERROR ("Total error is nan!!!\n");
      std::cerr << vars_abs.transpose () << std::endl;
      for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
      {
        Node* node_end = getNode (map_num_to_name_ [c_i]);
        Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
        std::cerr << "node " << c_i << ": " << node_end_position.block<3, 1> (0, 0).transpose () << " " << "joint pos: " << joint_pos[c_i].transpose () << std::endl;
      }

      std::cerr << weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm () << std::endl;
    }

    /// Check if the error is good
    if ((total_error_prev - total_error) / total_error_prev < 1e-8) //5e-3) //1e-4)
    {
      /// Error bad

      /// Take out the last incremental updates
      for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
        nodes[n_i]->angles -= vars_inc.block<3, 1> (3 * n_i, 0);
      getNode ("hip")->local_transformation.block<3, 1> (0, 3) -= vars_inc.block<3, 1> (3 * num_joints, 0);

      /// Increase the damping
      weight_damping *= 10.;

      /// Stop if we are already damping too much
      if (weight_damping > 1e11)
        break;

      /// Go back to solving the problem with a different damping weight
      continue;
    }
    else
    {
      /// Error good
      weight_damping /= 5.;
      total_error_prev = total_error;
      if (weight_damping < 1e-7)
        weight_damping = 1e-7;

      /// Continue with computing new jacobians
    }

//    PCL_ERROR ("Iteration: %zu - error %f,   damping %f\n", iter_i, total_error, weight_damping);

    J.setZero ();
    b.setZero ();

    /// Set up the constraints
//    PCL_ERROR ("joint_pos size: %ld\n", joint_pos.size ());
    for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
    {
      /// Get the end effector for this constraint
      Node* node_end = getNode (map_num_to_name_ [c_i]);

//      PCL_ERROR ("joint %ld node end size: %ld\n", c_i, node_end->path_from_root.size ());

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        Node *node_on_path = node_end->path_from_root[n_i];
        int node_id = node_on_path->id;

        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        Eigen::Matrix3f angles_rot_inv = (Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()) *
                                          Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()) *
                                          Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ())).matrix ().inverse ();
        Eigen::Matrix3f rot_mat_dx = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () * Rxdx *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
        Eigen::Matrix3f rot_mat_dy = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () * Rydy *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix ();
        Eigen::Matrix3f rot_mat_dz = angles_rot_inv *
                                     Eigen::AngleAxisf (node_on_path->angles (0), Eigen::Vector3f::UnitX ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (1), Eigen::Vector3f::UnitY ()).matrix () *
                                     Eigen::AngleAxisf (node_on_path->angles (2), Eigen::Vector3f::UnitZ ()).matrix () * Rzdz;

        J.block<3, 1> (3 * c_i, 3 * node_id + 0) += T_prev.block<3, 3> (0, 0) * rot_mat_dx * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 1) += T_prev.block<3, 3> (0, 0) * rot_mat_dy * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 2) += T_prev.block<3, 3> (0, 0) * rot_mat_dz * (T_next * Eigen::Vector4f (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      b.block<3, 1> (3 * c_i, 0) = node_end_position.block<3, 1> (0, 0) - joint_pos[c_i];

      J.block<3, 3> (3 * c_i, 3 * num_joints) = Eigen::Matrix3f::Identity ();
    }

#if DEBUG_JACOBIANS
    /// Debug the Jacobian with numerical differentiation
    Eigen::MatrixXf J_numdiff (Eigen::MatrixXf::Zero (3 * joint_pos.size (), 3 * num_joints + 3));
    const double epsilon = 1e-5;
    computeGlobalTransformations ();
    for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
    {
      Node* node_end = getNode (map_num_to_name_ [c_i]);
      Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
      for (size_t v_i = 0; v_i < 3 * num_joints + 3; v_i ++)
        J_numdiff.block <3, 1> (3 * c_i, v_i) -= (node_end_position.block<3, 1> (0, 0) - joint_pos[c_i]) / epsilon;
    }

    Eigen::VectorXf pose_current = this->getCurrentPose ();
    for (size_t v_i = 0; v_i < 3 * num_joints + 3; ++v_i)
    {
      Eigen::VectorXf pose_numdiff = pose_current;
      pose_numdiff (v_i) += epsilon;
      this->applyPose (pose_numdiff);
      computeGlobalTransformations ();

      for (size_t c_i = 0; c_i < joint_pos.size (); ++c_i)
      {
        Node* node_end = getNode (map_num_to_name_ [c_i]);
        Eigen::Vector4f node_end_position = node_end->global_transformation * Eigen::Vector4f (0., 0., 0., 1.);
        J_numdiff.block <3, 1> (3 * c_i, v_i) += (node_end_position.block<3, 1> (0, 0) - joint_pos[c_i]) / epsilon;
      }
    }
    PCL_ERROR ("DEBUG Jacobian IK skeleton diff: %f ratio %f --- J norm %f, J_numdiff norm %f\n",
               (J - J_numdiff).squaredNorm (), (J - J_numdiff).squaredNorm () / J_numdiff.squaredNorm (),
               J.squaredNorm (), J_numdiff.squaredNorm ());
    this->applyPose (pose_current);
#endif

    lhs = J.transpose () * J;
    rhs = J.transpose () * (-b);
    if (iter_i > 15)
    {
      lhs += weight_reg_zero * J_reg_zero;
      rhs += weight_reg_zero * J_reg_zero * (-vars_abs);
    }
    else
    {
      lhs += 1e3 * J_reg_zero;
      rhs += 1e3 * J_reg_zero * (-vars_abs);
    }

//    PCL_ERROR ("Sanity check %f vs %f\n",
//               total_error_prev, b.squaredNorm () + weight_reg_zero * (J_reg_zero * vars_abs).squaredNorm ());
  }

  total_ikjoint_iterations += iter_i;
  count_ikjoint_calls ++;
  PCL_ERROR ("IKJOINT average number of iterations: %f, error: %f, only joints: %f\n",
             static_cast<double> (total_ikjoint_iterations) / static_cast<double> (count_ikjoint_calls),
             total_error_prev, total_joints_error);
}
#endif

#if 0
void
af::Skeleton::IKWithJointPositions (std::vector<Eigen::Vector3f> &joint_pos)
{
  int num_joints = nodes.size ();

  /// TODO add regularization to keep the angles as close to the rest pose as possible
  Eigen::Matrix4f Rxdx, Rydy, Rzdz;
  Rxdx << 0., 0., 0., 0.,
          0., 0., -1., 0.,
          0., 1., 0., 0.,
          0., 0., 0., 0.;
  Rydy << 0., 0., 1., 0.,
          0., 0., 0., 0.,
          -1., 0., 0., 0.,
          0., 0., 0., 0.;
  Rzdz << 0.,-1., 0., 0.,
          1., 0., 0., 0.,
          0., 0., 0., 0.,
          0., 0., 0., 0.;
  std::cerr << "rxdx...:\n" << Rxdx << "\n" << Rydy << "\n" << Rzdz << "\n";
  std::vector<int> indices = {0, 9, 10};

  double error = 0.;
  for (size_t iter_i = 0; iter_i < 1000; ++iter_i)
  {
    computeGlobalTransformations ();

    Eigen::MatrixXd J (Eigen::MatrixXd::Zero (joint_pos.size () * 3, num_joints * 3 + 3));
    Eigen::VectorXd b (Eigen::VectorXd::Zero (joint_pos.size () * 3));

    /// Set up the constraints
    for (size_t i_i = 0; i_i < indices.size (); ++i_i)
    {
      size_t c_i = indices[i_i];
      /// Get the end effector for this constraint
      Node* node_end = getNode (map_num_to_name_ [c_i]);
//      std::cerr << "\nSetting up the jacobian for node " << node_end->name << " : ";

      /// Go through the kinematic chain up to this node and set the jacobians
      for (size_t n_i = 0; n_i < node_end->path_from_root.size (); ++n_i)
      {
        Node *node_on_path = node_end->path_from_root[n_i];
//        std::cerr << node_on_path->name << " -> ";

        int node_id = node_on_path->id;//map_name_to_num_[node_on_path->name];
        Eigen::Matrix4d T_prev, T_next;
        T_prev = node_on_path->global_transformation;
        T_next = node_on_path->global_transformation.inverse () * node_end->global_transformation;

        J.block<3, 1> (3 * c_i, 3 * node_id + 0) = (T_prev * Rxdx * T_next * Eigen::Vector4d (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 1) = (T_prev * Rydy * T_next * Eigen::Vector4d (0., 0., 0., 1.)).block<3, 1> (0, 0);
        J.block<3, 1> (3 * c_i, 3 * node_id + 2) = (T_prev * Rzdz * T_next * Eigen::Vector4d (0., 0., 0., 1.)).block<3, 1> (0, 0);
      }

      Eigen::Vector4d node_end_position = node_end->global_transformation * Eigen::Vector4d (0., 0., 0., 1.);
      b.block<3, 1> (3 * c_i, 0) = node_end_position.block<3, 1> (0, 0) - joint_pos[c_i];

      J.block<3, 3> (3 * c_i, 3 * num_joints) = Eigen::Matrix3d::Identity ();
    }
    error = b.squaredNorm ();
    if (iter_i == 0)
      PCL_ERROR ("IK with Joint positions optimization, initial error: %f ", error);

    PCL_ERROR ("IK error - iteration %zu -> %f\n", iter_i, b.squaredNorm ());

    /// Solve the system
    const double weight_damping = 1.;
    Eigen::VectorXd angles_inc = (J.transpose () * J +
                                  weight_damping * Eigen::MatrixXd::Identity (3 * num_joints + 3,
                                                                              3 * num_joints + 3)).ldlt ().solve (-J.transpose () * b);
//    std::cerr << "angles_inc: " << angles_inc.transpose () << std::endl;

    /// Put back the incremental angles
    for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
    {
      int node_id = n_i; //map_name_to_num_[nodes[n_i]->name];
      Eigen::Matrix3d rot_inc = (Eigen::AngleAxisf (angles_inc (3 * node_id + 0), Eigen::Vector3f::UnitX ()) *
                                 Eigen::AngleAxisf (angles_inc (3 * node_id + 1), Eigen::Vector3f::UnitY ()) *
                                 Eigen::AngleAxisf (angles_inc (3 * node_id + 2), Eigen::Vector3f::UnitZ ())).matrix ();

      nodes[n_i]->local_transformation.block<3, 3> (0, 0) = nodes[n_i]->local_transformation.block<3, 3> (0, 0) * rot_inc;
    }
    getNode ("hip")->local_transformation.block<3, 1> (0, 3) += angles_inc.block<3, 1> (3 * num_joints, 0);
  }

  PCL_ERROR (", final: %f\n", error);
}
#endif


bool
af::Skeleton::closeToNeutralPose (float threshold)
{
  float total_error = 0.;
  for (size_t n_i = 1; n_i < nodes.size (); ++n_i)
  {
    /// Get the minimal rotation angle from the stored angles for this node
    Eigen::Matrix3f rot_mat = (Eigen::AngleAxisf (nodes[n_i]->angles (0), Eigen::Vector3f::UnitX ()) *
                               Eigen::AngleAxisf (nodes[n_i]->angles (1), Eigen::Vector3f::UnitY ()) *
                               Eigen::AngleAxisf (nodes[n_i]->angles (2), Eigen::Vector3f::UnitZ ())).matrix ();
    Eigen::AngleAxisf angle_axis;
    angle_axis.fromRotationMatrix (rot_mat);
    float angle = angle_axis.angle ();

    total_error += angle * angle;
  }

  total_error /= static_cast<float> (nodes.size ());
  total_error = sqrt (total_error);

  PCL_ERROR ("Pose averaged deviation from zero: %f\n", total_error);

  return (total_error < threshold);
}


/*
/// Old version
void
af::Skeleton::scaleNiteToSLSkeleton (const std::vector<af::SkeletonState*> &skeleton_frames,
                                     af::SkeletonMesh &skeleton_mesh_scaled)
{
  /// Compute the average limb lengths and scale our rest skeleton
  std::vector<double> average_limb_size (skeleton_frames.front ()->joint_pos.size (), 0.);
  double average_shoulder_width = 0.,
         average_hip_width = 0.,
         average_trunk_side_length = 0.;

  int count_valid_frames = 0;
  for (size_t frame_i = 0; frame_i < skeleton_frames.size (); ++frame_i)
  {
    if (skeleton_frames[frame_i]->joint_pos.size () != 0)
      count_valid_frames ++;
    for (size_t i = 0; i < skeleton_frames[frame_i]->link.size (); ++i)
      average_limb_size[i] += (skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (0)] -
                               skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (1)]).norm ();

    average_shoulder_width += (skeleton_frames[frame_i]->joint_pos[2] -
                               skeleton_frames[frame_i]->joint_pos[3]).norm ();
    average_hip_width += (skeleton_frames[frame_i]->joint_pos[9] -
                          skeleton_frames[frame_i]->joint_pos[10]).norm ();
    average_trunk_side_length += ((skeleton_frames[frame_i]->joint_pos[1] -
                                   skeleton_frames[frame_i]->joint_pos[9]).norm () +
                                   (skeleton_frames[frame_i]->joint_pos[1] -
                                    skeleton_frames[frame_i]->joint_pos[10]).norm ()) / 2.;
  }
  for (size_t i = 0; i < skeleton_frames.back ()->link.size (); ++i)
    average_limb_size[i] /= static_cast<double> (count_valid_frames);
  average_shoulder_width /= static_cast<double> (count_valid_frames);
  average_hip_width /= static_cast<double> (count_valid_frames);
  average_trunk_side_length /= static_cast<double> (count_valid_frames);

  /// Scale the height of the body
  double average_trunk_height = sqrt (average_trunk_side_length * average_trunk_side_length - average_hip_width * average_hip_width / 4.);
  double template_trunk_height = fabs (skeleton_mesh_scaled.skeleton_.getNode ("abdomen")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("chest")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("lCollar")->local_transformation.block<3, 1> (0, 3).norm ());
//  double scale_height = 1.3 * (average_trunk_height - fabs (skeleton_mesh_scaled.skeleton_.getNode ("neck")->local_transformation (1, 3)))
//      / template_trunk_height;
//  double scale_height = 1.3 * average_trunk_height / template_trunk_height;
  double scale_height = 1.2 * average_trunk_height / template_trunk_height;
  PCL_ERROR ("---> body trunk scale height: %f\n", scale_height);

  skeleton_mesh_scaled.mesh_rest_->vertices_ *= scale_height;
  af::Mesh<float>::writeMeshOBJ ("mesh_vertices_scaled.obj", *skeleton_mesh_scaled.mesh_rest_);
  for (size_t n_i = 0; n_i < skeleton_mesh_scaled.skeleton_rest_.nodes.size (); ++n_i)
    skeleton_mesh_scaled.skeleton_rest_.nodes[n_i]->local_transformation.block<3, 1> (0, 3) *= scale_height;
  skeleton_mesh_scaled.skeleton_ = skeleton_mesh_scaled.skeleton_rest_;


  /// Scale only specific limbs
  af::Skeleton::Node *node = skeleton_mesh_scaled.skeleton_.getNode ("rForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[4] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rHand");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[6] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[5] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lHand");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[7] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShin");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[11] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rFoot");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[13] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShin");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[12] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lFoot");
  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[14] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("head");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[8] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("neck");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[1] / node->local_transformation.block<3, 1> (0, 3).norm ();

  /// Scale the width of the body - between lshoulder and rshoulder and between lThigh and rThigh
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShldr");
  node->local_transformation.block<3, 1> (0, 3) *= average_shoulder_width / 2.1 / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShldr");
  node->local_transformation.block<3, 1> (0, 3) *= average_shoulder_width / 2.1 / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lThigh");
  node->local_transformation.block<3, 1> (0, 3) *= average_hip_width / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rThigh");
  node->local_transformation.block<3, 1> (0, 3) *= average_hip_width / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();

//  node = skeleton_mesh_scaled.skeleton_.getNode ("abdomen");
//  node->local_transformation (1, 3) *= scale_height;
//  node = skeleton_mesh_scaled.skeleton_.getNode ("chest");
//  node->local_transformation (1, 3) *= scale_height;
//  node = skeleton_mesh_scaled.skeleton_.getNode ("head");
//  node->local_transformation (1, 3) *= scale_height;

//  node = skeleton_mesh_scaled.skeleton_.getNode ("lCollar");
//  node->local_transformation (1, 3) *= scale_height;
//  node = skeleton_mesh_scaled.skeleton_.getNode ("rCollar");
//  node->local_transformation (1, 3) *= scale_height;
//  node = skeleton_mesh_scaled.skeleton_.getNode ("neck");
//  node->local_transformation (1, 3) *= scale_height;
}
*/

/*
/// Working backup scaling
void
af::Skeleton::scaleNiteToSLSkeleton (const std::vector<af::SkeletonState*> &skeleton_frames,
                                     af::SkeletonMesh &skeleton_mesh_scaled)
{
  /// Compute the median limb lengths and scale our rest skeleton
  std::vector<std::vector<float> > limb_sizes (skeleton_frames.front ()->joint_pos.size ());
  std::vector<float> shoulder_widths, hip_widths, trunk_side_lengths;

  int count_valid_frames = 0;
  for (size_t frame_i = 0; frame_i < skeleton_frames.size (); ++frame_i)
  {
    if (skeleton_frames[frame_i]->joint_pos.size () != 0)
      count_valid_frames ++;

    for (size_t i = 0; i < skeleton_frames[frame_i]->link.size (); ++i)
      limb_sizes[i].push_back ((skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (0)] -
                               skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (1)]).norm ());

    shoulder_widths.push_back ((skeleton_frames[frame_i]->joint_pos[2] -
                               skeleton_frames[frame_i]->joint_pos[3]).norm ());
    hip_widths.push_back ((skeleton_frames[frame_i]->joint_pos[9] -
                          skeleton_frames[frame_i]->joint_pos[10]).norm ());
    trunk_side_lengths.push_back (((skeleton_frames[frame_i]->joint_pos[1] -
                                   skeleton_frames[frame_i]->joint_pos[9]).norm () +
                                  (skeleton_frames[frame_i]->joint_pos[1] -
                                   skeleton_frames[frame_i]->joint_pos[10]).norm ()) / 2.);
  }

  /// Sort all the data and take out the medians
  for (size_t i = 0; i < limb_sizes.size (); ++i)
    std::sort (limb_sizes[i].begin (), limb_sizes[i].end ());
  std::sort (shoulder_widths.begin (), shoulder_widths.end ());
  std::sort (hip_widths.begin (), hip_widths.end ());
  std::sort (trunk_side_lengths.begin (), trunk_side_lengths.end ());

  const size_t median_index = count_valid_frames / 2;

  /// Scale the height of the body
  double average_trunk_height = sqrt (trunk_side_lengths[median_index] * trunk_side_lengths[median_index] -
                                      hip_widths[median_index] * hip_widths[median_index] / 4.);
  double template_trunk_height = fabs (skeleton_mesh_scaled.skeleton_.getNode ("abdomen")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("chest")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("lCollar")->local_transformation.block<3, 1> (0, 3).norm ());
//  double scale_height = 1.3 * (average_trunk_height - fabs (skeleton_mesh_scaled.skeleton_.getNode ("neck")->local_transformation (1, 3)))
//      / template_trunk_height;
//  double scale_height = 1.3 * average_trunk_height / template_trunk_height;
  double scale_height = 1.2 * average_trunk_height / template_trunk_height;
  PCL_ERROR ("---> body trunk scale height: %f\n", scale_height);

  skeleton_mesh_scaled.mesh_rest_->vertices_ *= scale_height;
  af::Mesh<float>::writeMeshOBJ ("mesh_vertices_scaled.obj", *skeleton_mesh_scaled.mesh_rest_);
  for (size_t n_i = 0; n_i < skeleton_mesh_scaled.skeleton_rest_.nodes.size (); ++n_i)
    skeleton_mesh_scaled.skeleton_rest_.nodes[n_i]->local_transformation.block<3, 1> (0, 3) *= scale_height;
  skeleton_mesh_scaled.skeleton_ = skeleton_mesh_scaled.skeleton_rest_;


  /// Scale only specific limbs
  af::Skeleton::Node *node = skeleton_mesh_scaled.skeleton_.getNode ("rForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[4][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rHand");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[6][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[5][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lHand");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[7][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShin");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[11][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rFoot");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[13][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShin");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[12][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lFoot");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[14][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("head");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[8] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("neck");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[1] / node->local_transformation.block<3, 1> (0, 3).norm ();

  /// Scale the width of the body - between lshoulder and rshoulder and between lThigh and rThigh
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShldr");
  node->local_transformation.block<3, 1> (0, 3) *= shoulder_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShldr");
  node->local_transformation.block<3, 1> (0, 3) *= shoulder_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lThigh");
  node->local_transformation.block<3, 1> (0, 3) *= hip_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rThigh");
  node->local_transformation.block<3, 1> (0, 3) *= hip_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
}
*/



void
af::Skeleton::scaleNiteToSLSkeleton (const std::vector<af::SkeletonState*> &skeleton_frames,
                                     af::SkeletonMesh &skeleton_mesh_scaled)
{
  /// Compute the median limb lengths and scale our rest skeleton
  std::vector<std::vector<float> > limb_sizes (skeleton_frames.front ()->joint_pos.size ());
  std::vector<float> shoulder_widths, hip_widths, trunk_side_lengths;

  std::map<std::string, int> &map_name_to_num = skeleton_mesh_scaled.skeleton_.map_name_to_num_;

  int count_valid_frames = 0;
  for (size_t frame_i = 0; frame_i < skeleton_frames.size (); ++frame_i)
  {
    if (skeleton_frames[frame_i]->joint_pos.size () != 0)
      count_valid_frames ++;

    for (size_t i = 0; i < skeleton_frames[frame_i]->link.size (); ++i)
      limb_sizes[i].push_back ((skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (0)] -
                               skeleton_frames[frame_i]->joint_pos[skeleton_frames[frame_i]->link[i] (1)]).norm ());

    shoulder_widths.push_back ((skeleton_frames[frame_i]->joint_pos[map_name_to_num["lShldr"]] -
                               skeleton_frames[frame_i]->joint_pos[map_name_to_num["rShldr"]]).norm ());
    hip_widths.push_back ((skeleton_frames[frame_i]->joint_pos[map_name_to_num["lThigh"]] -
                          skeleton_frames[frame_i]->joint_pos[map_name_to_num["rThigh"]]).norm ());
    trunk_side_lengths.push_back (((skeleton_frames[frame_i]->joint_pos[map_name_to_num["neck"]] -
                                   skeleton_frames[frame_i]->joint_pos[map_name_to_num["lThigh"]]).norm () +
                                  (skeleton_frames[frame_i]->joint_pos[map_name_to_num["neck"]] -
                                   skeleton_frames[frame_i]->joint_pos[map_name_to_num["rThigh"]]).norm ()) / 2.);
  }

  /// Sort all the data and take out the medians
  for (size_t i = 0; i < limb_sizes.size (); ++i)
    std::sort (limb_sizes[i].begin (), limb_sizes[i].end ());
  std::sort (shoulder_widths.begin (), shoulder_widths.end ());
  std::sort (hip_widths.begin (), hip_widths.end ());
  std::sort (trunk_side_lengths.begin (), trunk_side_lengths.end ());

  const size_t median_index = count_valid_frames / 2;

  /// Scale the height of the body
  double average_trunk_height = sqrt (trunk_side_lengths[median_index] * trunk_side_lengths[median_index] -
                                      hip_widths[median_index] * hip_widths[median_index] / 4.);
  double template_trunk_height = fabs (skeleton_mesh_scaled.skeleton_.getNode ("abdomen")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("chest")->local_transformation.block<3, 1> (0, 3).norm ()) +
                                 fabs (skeleton_mesh_scaled.skeleton_.getNode ("lCollar")->local_transformation.block<3, 1> (0, 3).norm ());
//  double scale_height = 1.3 * (average_trunk_height - fabs (skeleton_mesh_scaled.skeleton_.getNode ("neck")->local_transformation (1, 3)))
//      / template_trunk_height;
//  double scale_height = 1.3 * average_trunk_height / template_trunk_height;
  double scale_height = 1.2 * average_trunk_height / template_trunk_height;
  PCL_ERROR ("---> body trunk scale height: %f\n", scale_height);

  skeleton_mesh_scaled.mesh_rest_->vertices_ *= scale_height;
  af::Mesh<float>::writeMeshOBJ ("mesh_vertices_scaled.obj", *skeleton_mesh_scaled.mesh_rest_);
  for (size_t n_i = 0; n_i < skeleton_mesh_scaled.skeleton_rest_.nodes.size (); ++n_i)
    skeleton_mesh_scaled.skeleton_rest_.nodes[n_i]->local_transformation.block<3, 1> (0, 3) *= scale_height;
  skeleton_mesh_scaled.skeleton_ = skeleton_mesh_scaled.skeleton_rest_;


  /// Scale only specific limbs
  af::Skeleton::Node *node = skeleton_mesh_scaled.skeleton_.getNode ("rForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["rForeArm"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rHand");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["rHand"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lForeArm");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["lForeArm"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lHand");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["lHand"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShin");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["rShin"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rFoot");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["rFoot"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShin");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["lShin"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lFoot");
  node->local_transformation.block<3, 1> (0, 3) *= limb_sizes[map_name_to_num["lFoot"]][median_index] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("head");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[8] / node->local_transformation.block<3, 1> (0, 3).norm ();
//  node = skeleton_mesh_scaled.skeleton_.getNode ("neck");
//  node->local_transformation.block<3, 1> (0, 3) *= average_limb_size[1] / node->local_transformation.block<3, 1> (0, 3).norm ();

  /// Scale the width of the body - between lshoulder and rshoulder and between lThigh and rThigh
  node = skeleton_mesh_scaled.skeleton_.getNode ("lShldr");
  node->local_transformation.block<3, 1> (0, 3) *= shoulder_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rShldr");
  node->local_transformation.block<3, 1> (0, 3) *= shoulder_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("lThigh");
  node->local_transformation.block<3, 1> (0, 3) *= hip_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
  node = skeleton_mesh_scaled.skeleton_.getNode ("rThigh");
  node->local_transformation.block<3, 1> (0, 3) *= hip_widths[median_index] / 2. / node->local_transformation.block<3, 1> (0, 3).norm ();
}







void
af::Skeleton::applyPose (Eigen::VectorXf &pose)
{
  if (pose.size () != nodes.size () * 3 + 3)
  {
    PCL_ERROR ("ERROR: Pose values size %zu != 3 * nodes size %zu\n", pose.rows (), 3 * nodes.size () + 3);
    return;
  }

  for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
    nodes[n_i]->angles = pose.block<3, 1> (3 * n_i, 0);

  getNode ("hip")->local_transformation.block<3, 1> (0, 3) = pose.block<3, 1> (3 * nodes.size (), 0);
}


Eigen::VectorXf
af::Skeleton::getCurrentPose ()
{
  Eigen::VectorXf pose (3 * nodes.size () + 3);
  for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
    pose.block<3, 1> (3 * n_i, 0) = nodes[n_i]->angles;

  pose.block<3, 1> (3 * nodes.size (), 0) = getNode ("hip")->local_transformation.block<3, 1> (0, 3);
  return (pose);
}


void
af::Skeleton::applyAngles (Eigen::VectorXf &angles)
{
  if (angles.size () != (nodes.size () - 1) * 3)
  {
    PCL_ERROR ("ERROR: Pose values size %zu != 3 * nodes size %zu\n", angles.rows (), 3 * (nodes.size () - 1));
    return;
  }

  for (size_t n_i = 1; n_i < nodes.size (); ++n_i)
    nodes[n_i]->angles = angles.block<3, 1> (3 * (n_i - 1), 0);
}


af::Mesh<float>::Ptr
af::Skeleton::generateMesh ()
{
  computeGlobalTransformations ();

  af::Mesh<float>::Ptr mesh_result (new af::Mesh<float> ());
  af::Mesh<float> mesh_aux;

  /// Generate spheres
  for (size_t n_i = 0; n_i < nodes.size (); ++n_i)
  {
    af::Mesh<float>::Ptr mesh_sphere = af::generateSphere (nodes[n_i]->global_transformation.block<3, 1> (0, 3));
    af::Mesh<float>::concatenateMeshes (*mesh_sphere, *mesh_result, mesh_aux);
    *mesh_result = mesh_aux;
  }

  /// Generate cylinders in-between
  for (size_t n_i = 1; n_i < nodes.size (); ++n_i)
  {
    Eigen::Vector3f orig = nodes[n_i]->global_transformation.block<3, 1> (0, 3);
    Eigen::Vector3f tip = nodes[n_i]->parent->global_transformation.block<3, 1> (0, 3);
    af::Mesh<float>::Ptr mesh_cyl = af::generateCylinder (orig, tip);
    af::Mesh<float>::concatenateMeshes (*mesh_cyl, *mesh_result, mesh_aux);
    *mesh_result = mesh_aux;
  }
  mesh_result->computeNormalsUsingConnectivity ();

  return (mesh_result);
}
