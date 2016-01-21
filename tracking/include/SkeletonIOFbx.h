#pragma once
#include "FbxIOBasic.h"
#include "skeleton.h"

namespace af
{

class SkeletonIOFbx : public FbxIOBasic{
private:
    Skeleton *skeletonTree;
public:
    SkeletonIOFbx(Skeleton& t){
        skeletonTree = &t;
    }

    void read (const std::string& filename)
    {
        bool initialize = initializeFBX(filename);
        bool internal_read = readInternal ();
        bool destroy = destroyFBX();

        if (!internal_read || !initialize || !destroy)
            std::cerr<<"FBX file " << filename.c_str () << " open failed";
    }

    void write(const String& /*path*/){
        bool OBJ_WRITE_NOT_IMPLEMENTED = false;
        assert(OBJ_WRITE_NOT_IMPLEMENTED);
    }


private:

    bool
    dfsRead (FbxNode *node, Skeleton* skeleton_tree)
    {
//      std::cerr << "node " << node->GetName () << " is of type " << node->GetTypeName ();

      /// Check if the node is already in the skeleton tree (not the case for roots)
      if (!skeleton_tree->hasNode (node->GetName ()))
      {
          FbxDouble3 fbx_vec;
#if 1
        fbx_vec = node->LclScaling.Get ();
        Eigen::Vector3f scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->LclRotation.Get ();
        Eigen::Vector3f rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->LclTranslation.Get ();
        Eigen::Vector3f translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);
#endif
        fbx_vec = node->GetScalingLimits ().GetMin ();
        Eigen::Vector3f min_scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->GetScalingLimits ().GetMax ();
        Eigen::Vector3f max_scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->GetRotationLimits ().GetMin ();
        Eigen::Vector3f min_rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->GetRotationLimits ().GetMax ();
        Eigen::Vector3f max_rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->GetTranslationLimits ().GetMin ();
        Eigen::Vector3f min_translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        fbx_vec = node->GetTranslationLimits ().GetMax ();
        Eigen::Vector3f max_translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        FbxAMatrix transf_fbx = node->EvaluateLocalTransform ();
        Eigen::Matrix4f transf;
        for (size_t r = 0; r < 4; ++r)
          for (size_t c = 0; c < 4; ++c)
            transf (r, c) = transf_fbx.Get (c, r);
//        transf.block<3, 1> (0, 3) /= 10.;

        fbx_vec = node->PreRotation.Get ();
        Eigen::Vector3f rot_pre (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

        Eigen::Matrix3f rotation_mat (Eigen::Matrix3f::Identity ());
        rotation_mat = Eigen::AngleAxisf (rot_pre (0) * M_PI / 180., Eigen::Vector3f::UnitX ()) *
                       Eigen::AngleAxisf (rot_pre (1) * M_PI / 180., Eigen::Vector3f::UnitY ()) *
                       Eigen::AngleAxisf (rot_pre (2) * M_PI / 180., Eigen::Vector3f::UnitZ ()) *
                       Eigen::AngleAxisf (rotation (0) * M_PI / 180., Eigen::Vector3f::UnitX ()) *
                       Eigen::AngleAxisf (rotation (1) * M_PI / 180., Eigen::Vector3f::UnitY ()) *
                       Eigen::AngleAxisf (rotation (2) * M_PI / 180., Eigen::Vector3f::UnitZ ());


        /// HACK
//        translation /= 10.;
        Skeleton::Node *skeleton_node = new Skeleton::Node (node->GetName (),
                                                            transf,
                                                            min_scaling, max_scaling,
                                                            min_rotation, max_rotation,
                                                            min_translation, max_translation);
        skeleton_tree->insertNode (skeleton_node);



//        std::cerr << "----> Node name: " << node->GetName () << std::endl;
//        printf ("scaling %f %f %f, rotation %f %f %f, translation %f %f %f\nmin_scaling %f %f %f, max_scaling %f %f %f\nmin_rotation %f %f %f, max_rotation %f %f %f\nmin_translation %f %f %f, max_translation %f %f %f\n\n",
//                scaling[0], scaling[1], scaling[2],
//                rotation[0], rotation[1], rotation[2],
//                translation[0], translation[1], translation[2],
//                min_scaling[0], min_scaling[1], min_scaling[2], max_scaling[0], max_scaling[1], max_scaling[2],
//                min_rotation[0], min_rotation[1], min_rotation[2], max_rotation[0], max_rotation[1], max_rotation[2],
//                min_translation[0], min_translation[1], min_translation[2], max_translation[0], max_translation[1], max_translation[2]);

        fbx_vec = node->PostRotation.Get ();
        Eigen::Vector3f rot_post (fbx_vec[0], fbx_vec[1], fbx_vec[2]);
//        std::cerr << "rotation pre: " << rot_pre.transpose () << std::endl;
//        std::cerr << "rotation post: " << rot_post.transpose () << std::endl;
//        std::cerr << "transform:\n" << transf << std::endl;
      }

      if (node->GetNodeAttribute () &&
          node->GetNodeAttribute ()->GetAttributeType () == FbxNodeAttribute::eSkeleton)
      {
        int count = node->GetChildCount ();
        for (size_t i = 0; i < count; ++i)
        {
          FbxNode *child_node = node->GetChild (i);
//          std::cerr << "node " << node->GetName () << " has child " << child_node->GetName ();

          FbxDouble3 fbx_vec = child_node->LclScaling.Get ();
          Eigen::Vector3f scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->LclRotation.Get ();
          Eigen::Vector3f rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->LclTranslation.Get ();
          Eigen::Vector3f translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetScalingLimits ().GetMin ();
          Eigen::Vector3f min_scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetScalingLimits ().GetMax ();
          Eigen::Vector3f max_scaling (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetRotationLimits ().GetMin ();
          Eigen::Vector3f min_rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetRotationLimits ().GetMax ();
          Eigen::Vector3f max_rotation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetTranslationLimits ().GetMin ();
          Eigen::Vector3f min_translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          fbx_vec = child_node->GetTranslationLimits ().GetMax ();
          Eigen::Vector3f max_translation (fbx_vec[0], fbx_vec[1], fbx_vec[2]);


          FbxAMatrix transf_fbx = child_node->EvaluateLocalTransform ();
          Eigen::Matrix4f transf;
          for (size_t r = 0; r < 4; ++r)
            for (size_t c = 0; c < 4; ++c)
              transf (r, c) = transf_fbx.Get (c, r);
//          transf.block<3, 1> (0, 3) /= 10.;

          fbx_vec = node->PreRotation.Get ();
          Eigen::Vector3f rot_pre (fbx_vec[0], fbx_vec[1], fbx_vec[2]);

          Eigen::Matrix3f rotation_mat (Eigen::Matrix3f::Identity ());
          rotation_mat = Eigen::AngleAxisf (rot_pre (0) * M_PI / 180., Eigen::Vector3f::UnitX ()) *
                         Eigen::AngleAxisf (rot_pre (1) * M_PI / 180., Eigen::Vector3f::UnitY ()) *
                         Eigen::AngleAxisf (rot_pre (2) * M_PI / 180., Eigen::Vector3f::UnitZ ()) *
                         Eigen::AngleAxisf (rotation (0) * M_PI / 180., Eigen::Vector3f::UnitX ()) *
                         Eigen::AngleAxisf (rotation (1) * M_PI / 180., Eigen::Vector3f::UnitY ()) *
                         Eigen::AngleAxisf (rotation (2) * M_PI / 180., Eigen::Vector3f::UnitZ ());


//          std::cerr << "----> Node name: " << child_node->GetName () << std::endl;
//          printf ("scaling %f %f %f, rotation %f %f %f, translation %f %f %f\nmin_scaling %f %f %f, max_scaling %f %f %f\nmin_rotation %f %f %f, max_rotation %f %f %f\nmin_translation %f %f %f, max_translation %f %f %f\n\n",
//                  scaling[0], scaling[1], scaling[2],
//                  rotation[0], rotation[1], rotation[2],
//                  translation[0], translation[1], translation[2],
//                  min_scaling[0], min_scaling[1], min_scaling[2], max_scaling[0], max_scaling[1], max_scaling[2],
//                  min_rotation[0], min_rotation[1], min_rotation[2], max_rotation[0], max_rotation[1], max_rotation[2],
//                  min_translation[0], min_translation[1], min_translation[2], max_translation[0], max_translation[1], max_translation[2]);
//          std::cerr << "transform\n" << transf << std::endl;
//          std::cerr << "rotation pre: " << rot_pre.transpose () << std::endl;


          /// HACK
//          translation /= 10.;
          Skeleton::Node *skeleton_node = new Skeleton::Node (child_node->GetName (),
                                                              transf,
                                                              min_scaling, max_scaling,
                                                              min_rotation, max_rotation,
                                                              min_translation, max_translation);
          skeleton_tree->insertNode (skeleton_node, node->GetName ());
          dfsRead (child_node, skeleton_tree);
        }
      }
	  return true;
    }


    bool
    readInternal ()
    {
        /// Getting the root note and the mesh
        FbxNode* rootNode = scene->GetRootNode();

        /// For Guillaume's kinectSDK skeleton
//        FbxMesh* meshFbx = rootNode->FindChild ("C_rig_GRP")->FindChild ("c_model_GRP")->FindChild ("C_body_GEO")->GetMesh();
//        FbxGeometry* geo = rootNode->FindChild ("C_rig_GRP")->FindChild ("c_model_GRP")->FindChild ("C_body_GEO")->GetGeometry();

        FbxMesh* meshFbx = rootNode->FindChild ("pythonMesh")->GetMesh();
        FbxGeometry* geo = rootNode->FindChild ("pythonMesh")->GetGeometry();


        /// Read in the joint tree
        dfsRead (rootNode->FindChild ("python")->FindChild ("hip"), skeletonTree);
//        dfsRead (rootNode->FindChild ("C_rig_GRP")->FindChild ("C_skeleton_GRP")->FindChild ("HipCenter"), skeletonTree);


        /// Now look into the deformations
//        qDebug()<<"Deformers in mesh:"<<meshFbx->GetDeformerCount();

        for (size_t d_i = 0; d_i < meshFbx->GetDeformerCount (); ++d_i)
        {
          FbxDeformer *deformer = geo->GetDeformer (d_i);

          if (deformer->GetDeformerType () == FbxDeformer::eSkin)
          {
            /// We checked above, so we can safely cast
            FbxSkin *skin = (FbxSkin*) deformer;
//            std::cerr << "skin deformer has #clusters = " << skin->GetClusterCount ();
//            std::cerr << "and control points: " << skin->GetControlPointIndicesCount ();

            for (size_t c_i = 0; c_i < skin->GetClusterCount (); ++c_i)
            {
              FbxCluster *cluster = skin->GetCluster (c_i);
//              std::cerr << "cluster " << c_i << " has #control points = " << cluster->GetControlPointIndicesCount ();
              int* control_point_indices = cluster->GetControlPointIndices ();
              double *control_point_weights = cluster->GetControlPointWeights ();

              FbxNode *node = cluster->GetLink ();
//              std::cerr << "the cluster is linked to the node with name: " << node->GetName () << " and type: " << node->GetTypeName () << std::endl;

              Skeleton::Node *skeleton_node = skeletonTree->getNode (node->GetName ());
              skeleton_node->vertex_indices = std::vector<int> (control_point_indices, control_point_indices + cluster->GetControlPointIndicesCount ());
              std::vector<double> aux = std::vector<double> (control_point_weights, control_point_weights + cluster->GetControlPointIndicesCount ());
              skeleton_node->vertex_weights.resize (aux.size ());
              for (size_t i = 0; i < aux.size (); ++i)
                skeleton_node->vertex_weights[i] = aux[i];
            }
          }
        }

        return (true);
    }
};


}

