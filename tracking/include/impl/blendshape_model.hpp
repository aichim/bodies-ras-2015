#pragma once

#include "common.h"
#include "blendshape_model.h"
#include <boost/filesystem.hpp>


/// Some re-declaration hacks
namespace af
{
void
getFilesFromDirectory (const std::string path_dir,
                       const std::string extension,
                       std::vector<std::string> &files);

std::string
getBaseName (const std::string &filename);
}



template <typename Scalar>
af::BlendshapeModel<Scalar>::BlendshapeModel ()
{
}


template <typename Scalar> bool
af::BlendshapeModel<Scalar>::load (const std::string &folder)
{
  std::vector<std::string> obj_files;
  af::getFilesFromDirectory (folder, ".OBJ", obj_files);

  bool found_neutral = false;
  for (size_t i = 0; i < obj_files.size (); ++i)
    if (af::getBaseName (obj_files[i]) == "neutral")
    {
      typename af::Mesh<Scalar>::Ptr mesh (new af::Mesh<Scalar> ());
      af::Mesh<Scalar>::readMeshOBJ (obj_files[i], *mesh);

      PCL_INFO ("Read neutral with %d vertices\n",  mesh->vertices_.cols ());

      neutral_mesh_ = mesh;
      neutral_ = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (neutral_mesh_->vertices_.data (),
                                                                                     neutral_mesh_->vertices_.cols () * 3, 1);
      exps_.resize (neutral_.rows (), obj_files.size () - 1);
      exp_names_.resize (obj_files.size () - 1);
      found_neutral = true;
      break;
    }

  if (found_neutral == false)
  {
    PCL_ERROR ("Did not find a neutral mesh in the given blendshape folder: %s\n", folder.c_str ());
    return (false);
  }

  int exp_index = 0;
  for (size_t i = 0; i < obj_files.size (); ++i)
  {
    if (af::getBaseName (obj_files[i]) != "neutral")
    {
      typename af::Mesh<Scalar>::Ptr mesh (new af::Mesh<Scalar> ());
      af::Mesh<Scalar>::readMeshOBJ (obj_files[i], *mesh);

      Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>  diff = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (mesh->vertices_.data (),
                                                                                                                                        exps_.rows (), 1) - neutral_;
      exps_.block (0, exp_index, exps_.rows (), 1) = diff;

      exp_names_[exp_index] = af::getBaseName (obj_files[i]);
      exp_index ++;

//      PCL_INFO ("Read expression %s with %d vertices\n", exp_names_[exp_index-1].c_str (), mesh->vertices_.cols ());
    }
  }

  PCL_INFO ("--- Blendshape exps rows: %d, cols %d\n", exps_.rows (), exps_.cols ());

  num_pts_ = neutral_mesh_->vertices_.cols ();
  num_exps_ = exp_names_.size ();

  return (true);
}


template <typename Scalar> bool
af::BlendshapeModel<Scalar>::loadBinary (const std::string &filename)
{
  std::ifstream file (filename.c_str (), std::ios::binary);
  if (!file.is_open ())
  {
    PCL_ERROR ("%s: could not open file to load blendshape model: %s\n", __func__, filename.c_str ());
    return (false);
  }

  neutral_mesh_.reset (new af::Mesh<Scalar> ());

  /// Read the neutral mesh
  int num_vertices;
  file.read (reinterpret_cast<char*> (&num_vertices), sizeof (int));
  neutral_mesh_->vertices_.resize (Eigen::NoChange_t (), num_vertices);
  file.read (reinterpret_cast<char*> (neutral_mesh_->vertices_.data ()), 3 * num_vertices * sizeof (double));

  int num_tris;
  file.read (reinterpret_cast<char*> (&num_tris), sizeof (int));
  neutral_mesh_->tri_vertices_.resize (Eigen::NoChange_t (), num_tris);
  file.read (reinterpret_cast<char*> (neutral_mesh_->tri_vertices_.data ()), 3 * num_tris * sizeof (size_t));

  int num_tex_tris;
  file.read (reinterpret_cast<char*> (&num_tex_tris), sizeof (int));
  neutral_mesh_->tri_tex_coords_.resize (Eigen::NoChange_t (), num_tex_tris);
  file.read (reinterpret_cast<char*> (neutral_mesh_->tri_tex_coords_.data ()), 3 * num_tex_tris * sizeof (size_t));

  int num_tex_coords;
  file.read (reinterpret_cast<char*> (&num_tex_coords), sizeof (int));
  neutral_mesh_->tex_coords_.resize (Eigen::NoChange_t (), num_tex_coords);
  file.read (reinterpret_cast<char*> (neutral_mesh_->tex_coords_.data ()), 2 * num_tex_coords * sizeof (double));

  neutral_ = Eigen::Map<Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> > (neutral_mesh_->vertices_.data (),
                                          neutral_mesh_->vertices_.cols () * 3, 1);

  /// Now read all the blendshape diffs
  file.read (reinterpret_cast<char*> (&num_exps_), sizeof (int));
  exps_.resize (3 * num_vertices, num_exps_);
  file.read (reinterpret_cast<char*> (exps_.data ()), 3 * num_vertices * num_exps_ * sizeof (double));

  /// And the blendshape names
  for (size_t bs_i = 0; bs_i < num_exps_; ++bs_i)
  {
    int string_len;
    file.read (reinterpret_cast<char*> (&string_len), sizeof (int));
    std::string name;
    name.resize (string_len);
    file.read (reinterpret_cast<char*> (&name[0]), string_len * sizeof (char));
    exp_names_.push_back (name);
  }

  return (true);
}


template <typename Scalar> bool
af::BlendshapeModel<Scalar>::save (const std::string &folder) const
{
  boost::filesystem::create_directory (folder);

  af::Mesh<Scalar>::writeMeshOBJ (folder + "/neutral.obj", *neutral_mesh_);

  for (size_t i = 0; i < exp_names_.size (); ++i)
  {
    af::Mesh<Scalar> mesh = *neutral_mesh_;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> weights = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>::Zero (num_exps_);
    weights (i) = 1.;
    evaluate (weights, mesh.vertices_);
    mesh.computeNormalsUsingConnectivity ();
    af::Mesh<Scalar>::writeMeshOBJ (folder + "/" + exp_names_[i] + ".obj", mesh);
  }


  /// TODO finish this tomorrow, the lover wants to sleep ...

  return (true);
}


template <typename Scalar> bool
af::BlendshapeModel<Scalar>::saveBinary (const std::string &filename) const
{
  std::ofstream file (filename.c_str (), std::ios::trunc | std::ios::binary);
  if (!file.is_open ())
  {
    PCL_ERROR ("%s: could not open file to write blendshape model: %s\n", __func__, filename.c_str ());
    return (false);
  }

  /// Write the neutral mesh first.
  int num_vertices = neutral_mesh_->vertices_.cols ();
  file.write (reinterpret_cast<const char*> (&num_vertices), sizeof (int));
  file.write (reinterpret_cast<const char*> (neutral_mesh_->vertices_.data ()), 3 * num_vertices * sizeof (double));

  int num_tris = neutral_mesh_->tri_vertices_.cols ();
  file.write (reinterpret_cast<const char*> (&num_tris), sizeof (int));
  file.write (reinterpret_cast<const char*> (neutral_mesh_->tri_vertices_.data ()), 3 * num_tris * sizeof (size_t));

  int num_tex_tris = neutral_mesh_->tri_tex_coords_.cols ();
  file.write (reinterpret_cast<const char*> (&num_tex_tris), sizeof (int));
  file.write (reinterpret_cast<const char*> (neutral_mesh_->tri_tex_coords_.data ()), 3 * num_tex_tris * sizeof (size_t));

  int num_tex_coords = neutral_mesh_->tex_coords_.cols ();
  file.write (reinterpret_cast<const char*> (&num_tex_coords), sizeof (int));
  file.write (reinterpret_cast<const char*> (neutral_mesh_->tex_coords_.data ()), 2 * num_tex_coords * sizeof (double));


  /// Now write all the blendshape diffs
  file.write (reinterpret_cast<const char*> (&num_exps_), sizeof (int));
  file.write (reinterpret_cast<const char*> (exps_.data ()), 3 * num_vertices * num_exps_ * sizeof (double));

  /// And the blendshape names
  for (size_t bs_i = 0; bs_i < num_exps_; ++bs_i)
  {
    int string_len = exp_names_[bs_i].length ();
    file.write (reinterpret_cast<const char*> (&string_len), sizeof (int));
    file.write (reinterpret_cast<const char*> (exp_names_[bs_i].c_str ()), string_len * sizeof (char));
  }

  return (true);
}



template <typename Scalar> void
af::BlendshapeModel<Scalar>::deepCopy (af::BlendshapeModel<Scalar> &copy) const
{
  copy.neutral_mesh_.reset (new af::Mesh<Scalar> ());
  *copy.neutral_mesh_ = *neutral_mesh_;
  copy.neutral_ = neutral_;

  copy.exp_names_ = exp_names_;
  copy.exps_ = exps_;

  copy.num_exps_ = num_exps_;
  copy.num_pts_ = num_pts_;
}


template <typename Scalar> bool
af::BlendshapeModel<Scalar>::evaluate (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1> &coeffs,
                                       Eigen::Matrix<Scalar, 3, Eigen::Dynamic> &points) const
{
  if (coeffs.rows () != exps_.cols ())
  {
    PCL_ERROR ("Error: the number of coefficients %ld != number of expressions %ld\n",
               coeffs.rows (), exps_.cols ());
    return (false);
  }

  size_t num_pts = neutral_.rows () / 3;
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> aux = neutral_ + exps_ * coeffs;
  points = Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic> > (aux.data (), 3, num_pts);
  return (true);
}


template <typename Scalar> void
af::BlendshapeModel<Scalar>::convertToTriMesh ()
{
  af::Mesh<Scalar> neutral_tri;
  af::Mesh<Scalar>::convertToTriangleMesh (*neutral_mesh_, neutral_tri);
  *neutral_mesh_ = neutral_tri;
}


