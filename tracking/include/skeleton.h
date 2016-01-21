#pragma once

#include <Eigen/Eigen>
#include <deque>

#include "bodies_utils.h"

namespace af
{

class SkeletonMesh;

class Skeleton{
    friend class SkeletonMesh;
    enum NodeError
    {
      ROOT_MUST_BE_INSERTED_FIRST = 0,
      NODE_MUST_HAVE_A_NAME = 0,
      ATTEMPTED_INSERT_DUPLICATE_NODE = 0
    };

public:
    class Node{
      public:
        std::string name;
        int id;
        Node* parent;
        std::vector<Node*> children;

        Eigen::Matrix4f local_transformation, global_transformation;
        Eigen::Vector3f angles;

        Eigen::Vector3f min_scaling, max_scaling, min_rotation, max_rotation, min_translation, max_translation;
        std::vector<int> vertex_indices;
        std::vector<float> vertex_weights;


        std::vector<Node*> path_from_root;

        Node () {
          angles = Eigen::Vector3f::Zero ();
          local_transformation = global_transformation = Eigen::Matrix4f::Identity ();
        }

        /**
         * @brief Node - copies all the members except for the parent node (sets it to NULL)
         * @param copy
         */
        Node (const Node &copy);


        Node (std::string _name,
              Eigen::Matrix4f &_local_transformation,
              Eigen::Vector3f &_min_scaling,
              Eigen::Vector3f &_max_scaling,
              Eigen::Vector3f &_min_rotation,
              Eigen::Vector3f &_max_rotation,
              Eigen::Vector3f &_min_translation,
              Eigen::Vector3f &_max_translation);

        std::string getNodeName(){return name;}
        std::vector<Node*> getAllChildren(){ return children; }
        Node* getChildAt(int i){return children[i];}
        void insertChild(Node *node){children.push_back(node);}

        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    };

//private:
    std::vector<Node*> nodes;
    std::map<std::string,Node*> name_to_node;

public:
    Skeleton ();

    Skeleton& operator =(const Skeleton &copy);

    Skeleton (const Skeleton &copy);

    Node* getRoot(){ return (nodes.size()==0)?NULL:nodes[0]; }

    Node* getNode(std::string name){
        std::map<std::string,Node*>::iterator it = name_to_node.find(name);
        assert(it!=name_to_node.end());
        return it->second;
    }
    
    bool hasNode(std::string name){
        std::map<std::string,Node*>::iterator it = name_to_node.find(name);
        return (it!=name_to_node.end());
    }
    
    void insertNode(Node* new_node, std::string parent_name="");

    void
    computeGlobalTransformations ();

    void
    computePathsDFS (Node *node);

    void
    clearPaths ();

    void
    bindJointNames();

    void
    IKWithJointPositions (const std::vector<Eigen::Vector3f> &joint_pos,
                          const double weight_reg_zero = 1e-4);

    static Eigen::Matrix4f
    jacrot (const Eigen::Vector4d &p);


    bool
    closeToNeutralPose (float threshold = 10. * M_PI / 180.);


    static void
    scaleNiteToSLSkeleton (const std::vector<af::SkeletonState*> &skeleton_frames,
                           af::SkeletonMesh &skeleton_mesh_scaled);

    /**
     * @brief applyPose - apply a complete pose, including the root rotation and translation
     * @param pose
     */
    void
    applyPose (Eigen::VectorXf &pose);

    void
    applyPose (Eigen::VectorXd &pose)
    {
      Eigen::VectorXf pose_f = pose.cast<float> ();
      applyPose (pose_f);
    }

    Eigen::VectorXf
    getCurrentPose ();

    /**
     * @brief applyAngles - apply only the internal angles (i.e., do not touch the root node)
     * @param angles
     */
    void
    applyAngles (Eigen::VectorXf &angles);

    void
    applyAngles (Eigen::VectorXd &angles)
    {
      Eigen::VectorXf angles_f = angles.cast<float> ();
      applyAngles (angles_f);
    }

    af::Mesh<float>::Ptr
    generateMesh ();

    std::map<int, std::string> map_num_to_name_;
    std::map<std::string, int> map_name_to_num_;

    static Eigen::Matrix3f Rxdx, Rydy, Rzdz;
};


}

