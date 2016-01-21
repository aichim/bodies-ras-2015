#pragma once

#include <map>
#include <pcl/common/common.h>
#include <pcl/console/print.h>

namespace af
{

template <class Scalar>
struct Material
{
  Material ()
    : color_ambient (Eigen::Matrix<Scalar, 3, 1> (0.5, 0.5, 0.5))
    , color_diffuse (Eigen::Matrix<Scalar, 3, 1> (1.0, 1.0, 1.0))
    , color_specular (Eigen::Matrix<Scalar, 3, 1> (1.0, 1.0, 1.0))
    , specular_coeff (0.)
    , transparency (0.)
    , name ()
  {
  }

  Eigen::Matrix<Scalar, 3, 1> color_ambient, color_diffuse, color_specular;
  Scalar specular_coeff, transparency;

  std::string name;
  std::string texture_diffuse, texture_specular, texture_normal;
};

struct MeshSegment
{
  std::vector<size_t> tri_indices;
  std::vector<size_t> quad_indices;

  inline void clear ()
  { tri_indices.clear (); quad_indices.clear (); }

  inline size_t size () const
  { return (tri_indices.size ()); }
};


struct Texel
{
  Texel (int u_arg, int v_arg, int tri_id_arg, Eigen::Vector3f &bary_arg)
    : u (u_arg), v (v_arg)
    , tri_id (tri_id_arg)
    , bary (bary_arg)
  {}

  int u, v;
  int tri_id;
  Eigen::Vector3f bary;
};


template <class Scalar = double>
class Mesh
{
public:
  typedef typename boost::shared_ptr<Mesh<Scalar> > Ptr;
  typedef typename boost::shared_ptr<const Mesh<Scalar> > ConstPtr;

  Mesh ();

  Mesh (const Mesh<Scalar> &mesh);


  /**
   * @brief computeNormalsUsingConnectivity
   * @note Warning! Will overwrite the original normals if any
   */
  void
  computeNormalsUsingConnectivity ();

  void
  printInfo () const;

  void
  computeNormals ();

  Ptr
  makeShared ()
  {
    return (Ptr (new Mesh<Scalar> (*this)));
  }

  void
  centerIt ();


  static bool
  readMeshOBJ (const std::string &filename,
               Mesh<Scalar> &mesh);

  static bool
  writeMeshOBJ (const std::string &filename,
                const Mesh<Scalar> &mesh);

  static bool
  readMaterialFile (const std::string &filename,
                    std::map<std::string, Material<Scalar> > &materials);

  static bool
  writeMaterialFile (const std::string &filename,
                     const std::map<std::string, Material<Scalar> > &materials);

  static void
  convertToTriangleMesh (const Mesh<Scalar> &mesh,
                         Mesh<Scalar> &result);

  static void
  computeTexelGeometry (const Mesh<Scalar> &mesh,
                        std::vector<af::Texel> &texels,
                        int texture_size);

  void
  stripOutNormals ();

  static void
  concatenateMeshes (const Mesh<Scalar> &mesh_a, const Mesh<Scalar> &mesh_b, Mesh<Scalar> &result);


//  protected:
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> vertices_;
    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> normals_;
    Eigen::Matrix<Scalar, 2, Eigen::Dynamic> tex_coords_;


    Eigen::Matrix<size_t, 3, Eigen::Dynamic> tri_vertices_;
    Eigen::Matrix<size_t, 3, Eigen::Dynamic> tri_normals_;
    Eigen::Matrix<size_t, 3, Eigen::Dynamic> tri_tex_coords_;

    Eigen::Matrix<size_t, 4, Eigen::Dynamic> quad_vertices_;
    Eigen::Matrix<size_t, 4, Eigen::Dynamic> quad_normals_;
    Eigen::Matrix<size_t, 4, Eigen::Dynamic> quad_tex_coords_;

    std::map<std::string, MeshSegment> segments_;
    std::map<std::string, Material<Scalar> > materials_;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}

#include "impl/mesh.hpp"
