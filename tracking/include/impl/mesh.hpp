#pragma once

#include "mesh.h"

#include <fstream>
#include <boost/foreach.hpp>
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <opencv2/opencv.hpp>


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar>
af::Mesh<Scalar>::Mesh ()
{
}


////////////////////////////////////////////////////////////////////////////////
template<typename Scalar>
af::Mesh<Scalar>::Mesh (const Mesh<Scalar> &mesh)
{
  vertices_ = mesh.vertices_;
  normals_ = mesh.normals_;
  tex_coords_ = mesh.tex_coords_;

  tri_vertices_ = mesh.tri_vertices_;
  tri_normals_ = mesh.tri_normals_;
  tri_tex_coords_ = mesh.tri_tex_coords_;

  quad_vertices_ = mesh.quad_vertices_;
  quad_normals_ = mesh.quad_normals_;
  quad_tex_coords_ = mesh.quad_tex_coords_;

  segments_ = mesh.segments_;
  materials_ = mesh.materials_;
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> void
af::Mesh<Scalar>::computeNormalsUsingConnectivity ()
{
  // Update normals
  normals_.resize (3, vertices_.cols ());
  normals_.setZero (normals_.rows (), normals_.cols ());

  tri_normals_ = tri_vertices_;
  quad_normals_ = quad_vertices_;

  // For triangles
  for (int tri_i = 0; tri_i < tri_vertices_.cols (); ++tri_i)
  {
    Eigen::Matrix<Scalar, 3, 1> v1 = vertices_.col (tri_vertices_ (1, tri_i))  - vertices_.col (tri_vertices_ (0, tri_i));
    Eigen::Matrix<Scalar, 3, 1> v2 = vertices_.col (tri_vertices_ (2, tri_i)) - vertices_.col (tri_vertices_ (0, tri_i));
    Eigen::Matrix<Scalar, 3, 1> normal = (v1).cross (v2);

    normals_.col (tri_vertices_ (0, tri_i)) += normal;
    normals_.col (tri_vertices_ (1, tri_i)) += normal;
    normals_.col (tri_vertices_ (2, tri_i)) += normal;
  }


  // For quads
  for (int quad_i = 0; quad_i < quad_vertices_.cols (); ++quad_i)
  {
    Eigen::Matrix<Scalar, 3, 1> v1 = vertices_.col (quad_vertices_ (1, quad_i)) - vertices_.col (quad_vertices_ (0, quad_i));
    Eigen::Matrix<Scalar, 3, 1> v2 = vertices_.col (quad_vertices_ (3, quad_i)) - vertices_.col (quad_vertices_ (0, quad_i));
    Eigen::Matrix<Scalar, 3, 1> normal = (v1).cross (v2);

    normals_.col (quad_vertices_ (0, quad_i)) += normal * 2.;
    normals_.col (quad_vertices_ (1, quad_i)) += normal * 2.;
    normals_.col (quad_vertices_ (2, quad_i)) += normal * 2.;
    normals_.col (quad_vertices_ (3, quad_i)) += normal * 2.;
  }

  // Normalize the normals
  for (int n_i = 0; n_i < normals_.cols (); ++n_i)
    normals_.col (n_i).normalize ();
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> void
af::Mesh<Scalar>::centerIt ()
{
  Eigen::Matrix<Scalar, 3, 1> centroid;
  for (size_t i = 0; i < 3; ++i)
    centroid (i) = vertices_.row (i).sum () / static_cast<Scalar> (vertices_.cols ());

  for (size_t v_i = 0; v_i < vertices_.cols (); ++v_i)
    vertices_.col (v_i) -= centroid;
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> bool
af::Mesh<Scalar>::readMeshOBJ (const std::string &filename,
                               af::Mesh<Scalar> &mesh)
{
  typedef Eigen::Matrix<Scalar, 3, 1> Vector3;
  typedef Eigen::Matrix<Scalar, 2, 1> Vector2;
  typedef Eigen::Matrix<size_t, 3, 1> Vector3_size_t;
  typedef Eigen::Matrix<size_t, 4, 1> Vector4_size_t;
  typedef typename std::vector<Vector3, Eigen::aligned_allocator<Vector3> > VectorOfVector3;
  typedef typename std::vector<Vector2, Eigen::aligned_allocator<Vector2> > VectorOfVector2;
  typedef typename std::vector<Vector3_size_t, Eigen::aligned_allocator<Vector3_size_t> > VectorOfVector3i;
  typedef typename std::vector<Vector4_size_t, Eigen::aligned_allocator<Vector4_size_t> > VectorOfVector4i;

  std::string folder_name = filename.substr (0, filename.find_last_of ("/")) + "/";

  std::ifstream file (filename.c_str ());
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening mesh file %s\n", filename.c_str ());
    return (false);
  }


  VectorOfVector3 vertices_aux;
  VectorOfVector3 normals_aux;
  VectorOfVector2 tex_coords_aux;

  VectorOfVector3i tri_vertices_aux;
  VectorOfVector3i tri_normals_aux;
  VectorOfVector3i tri_tex_coords_aux;

  VectorOfVector4i quad_vertices_aux;
  VectorOfVector4i quad_normals_aux;
  VectorOfVector4i quad_tex_coords_aux;


  std::string line, line_type;
  std::string current_material_name = "";
  bool have_material = false;
  MeshSegment current_segment;
  while (1)
  {
    std::getline (file, line);
    if (!file.good ())
      break;

    if (line == "")
      continue;

    std::stringstream line_ss (line);
    line_ss >> line_type;
    boost::algorithm::to_lower (line_type);

    if (line_type == "#")
      continue;
    else if (line_type == "mtllib")
    {
      std::string material_filename;
      line_ss >> material_filename;
      material_filename = folder_name + material_filename;

      readMaterialFile (material_filename, mesh.materials_);
    }
    else if (line_type == "usemtl")
    {
      if (have_material)
      {
        mesh.segments_.insert (std::make_pair (current_material_name, current_segment));
        current_segment.clear ();
      }
      line_ss >> current_material_name;
      have_material = true;
    }
    else if (line_type == "v")
    {
      Vector3 vertex;
      line_ss >> vertex[0] >> vertex[1] >> vertex[2];
      vertices_aux.push_back (vertex);
    }
    else if (line_type == "vn")
    {
      Vector3 normal;
      line_ss >> normal[0] >> normal[1] >> normal[2];
      normals_aux.push_back (normal);
    }
    else if (line_type == "vt")
    {
      Vector2 tex_coord;
      line_ss >> tex_coord[0] >> tex_coord[1];
      tex_coords_aux.push_back (tex_coord);
    }
    else if (line_type == "f")
    {
      boost::char_separator<char> sep (" ");
      std::string remainder;
      std::getline (line_ss, remainder);
      boost::tokenizer<boost::char_separator<char> > tokens (remainder, sep);
      int tokens_size = std::distance (tokens.begin (), tokens.end ());
      if (tokens_size != 3 &&
          tokens_size != 4)
      {
        PCL_ERROR ("Error: non-triangular or non-quad face found, #vertices = %d.\n", tokens_size);
        return (false);
      }

      if (tokens_size == 3)
      {
        Eigen::Matrix<size_t, 3, 1> vertex_indices (0, 0, 0),
            normal_indices (0, 0, 0),
            tex_coords_indices (0, 0, 0);

        size_t index = 0;
        BOOST_FOREACH (const std::string &word, tokens)
        {
          boost::char_separator<char> property_sep ("/");
          boost::tokenizer<boost::char_separator<char> > property_tokens (word, property_sep);

          boost::tokenizer<boost::char_separator<char> >::const_iterator p_it = property_tokens.begin ();

          size_t property_tokens_size = std::distance (property_tokens.begin (), property_tokens.end ());

          if (property_tokens_size == 1)
            /// Only vertex indices
          {
            std::stringstream ss (*p_it);
            ss >> vertex_indices[index];
          }
          else if (property_tokens_size == 2)
            /// Vertex and texture indices
          {
            std::stringstream ss (*p_it);
            ss >> vertex_indices[index];

            ss.clear ();
            ss.str (*(++p_it));
            ss >> tex_coords_indices[index];
          }
          else if (property_tokens_size == 3)
            /// Vertex, texture, and normal indices
          {
            std::stringstream ss;
            if (p_it->size () != 0)
            {
              ss.str (*p_it);
              ss >> vertex_indices[index];
            }

            p_it ++;
            if (p_it->size () != 0)
            {
              ss.clear ();
              ss.str (*p_it);
              ss >> tex_coords_indices[index];
            }

            p_it ++;
            if (p_it->size () != 0)
            {
              ss.clear ();
              ss.str (*p_it);
              ss >> normal_indices[index];
            }
          }

//          if (have_material)
            index ++;
        }

        tri_vertices_aux.push_back (vertex_indices);

        /// If we read valid normal indices
        if (normal_indices[0] != 0 &&
            normal_indices[1] != 0 &&
            normal_indices[2] != 0)
          tri_normals_aux.push_back (normal_indices);

        /// If we read valid texture coordinate indices
        if (tex_coords_indices[0] != 0 &&
            tex_coords_indices[1] != 0 &&
            tex_coords_indices[2] != 0)
          tri_tex_coords_aux.push_back (tex_coords_indices);

        /// Add the current face to its corresponding segment
        current_segment.tri_indices.push_back (tri_vertices_aux.size () - 1);
      }
      else if (tokens_size == 4)
      {
        Eigen::Matrix<size_t, 4, 1> vertex_indices (0, 0, 0, 0),
            normal_indices (0, 0, 0, 0),
            tex_coords_indices (0, 0, 0, 0);

        size_t index = 0;
        BOOST_FOREACH (const std::string &word, tokens)
        {
          boost::char_separator<char> property_sep ("/");
          boost::tokenizer<boost::char_separator<char> > property_tokens (word, property_sep);

          boost::tokenizer<boost::char_separator<char> >::const_iterator p_it = property_tokens.begin ();

          size_t property_tokens_size = std::distance (property_tokens.begin (), property_tokens.end ());

          if (property_tokens_size == 1)
            /// Only vertex indices
          {
            std::stringstream ss (*p_it);
            ss >> vertex_indices[index];
          }
          else if (property_tokens_size == 2)
            /// Vertex and texture indices
          {
            std::stringstream ss (*p_it);
            ss >> vertex_indices[index];

            ss.clear ();
            p_it ++;
            if (p_it->size () != 0)
            {
              ss.str (*p_it);
              ss >> tex_coords_indices[index];
            }
          }
          else if (property_tokens_size == 3)
            /// Vertex, texture, and normal indices
          {
            std::stringstream ss;
            if (p_it->size () != 0)
            {
              ss.str (*p_it);
              ss >> vertex_indices[index];
            }

            p_it ++;
            if (p_it->size () != 0)
            {
              ss.clear ();
              ss.str (*p_it);
              ss >> tex_coords_indices[index];
            }

            p_it ++;
            if (p_it->size () != 0)
            {
              ss.clear ();
              ss.str (*p_it);
              ss >> normal_indices[index];
            }
          }

//          if (have_material)
            index ++;
        }

        quad_vertices_aux.push_back (vertex_indices);

        /// If we read valid normal indices
        if (normal_indices[0] != 0 &&
            normal_indices[1] != 0 &&
            normal_indices[2] != 0 &&
            normal_indices[3] != 0)
          quad_normals_aux.push_back (normal_indices);

        /// If we read valid texture coordinate indices
        if (tex_coords_indices[0] != 0 &&
            tex_coords_indices[1] != 0 &&
            tex_coords_indices[2] != 0 &&
            tex_coords_indices[3] != 0)
          quad_tex_coords_aux.push_back (tex_coords_indices);

        /// Add the current face to its corresponding segment
        current_segment.quad_indices.push_back (quad_vertices_aux.size () - 1);
      }
    }
  }

  if (have_material)
    mesh.segments_.insert (std::make_pair (current_material_name, current_segment));

  file.close ();


  /// Now convert everything to the Eigen format
  mesh.vertices_ = Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic> > (&vertices_aux[0][0], 3, vertices_aux.size ());
  mesh.normals_ = Eigen::Map<Eigen::Matrix<Scalar, 3, Eigen::Dynamic> > (&normals_aux[0][0], 3, normals_aux.size ());
  mesh.tex_coords_ = Eigen::Map<Eigen::Matrix<Scalar, 2, Eigen::Dynamic> > (&tex_coords_aux[0][0], 2, tex_coords_aux.size ());


  mesh.tri_vertices_ = Eigen::Map<Eigen::Matrix<size_t, 3, Eigen::Dynamic> > (&tri_vertices_aux[0][0], 3, tri_vertices_aux.size ());
  mesh.tri_normals_ = Eigen::Map<Eigen::Matrix<size_t, 3, Eigen::Dynamic> > (&tri_normals_aux[0][0], 3, tri_normals_aux.size ());
  mesh.tri_tex_coords_ = Eigen::Map<Eigen::Matrix<size_t, 3, Eigen::Dynamic> > (&tri_tex_coords_aux[0][0], 3, tri_tex_coords_aux.size ());

  mesh.quad_vertices_ = Eigen::Map<Eigen::Matrix<size_t, 4, Eigen::Dynamic> > (&quad_vertices_aux[0][0], 4, quad_vertices_aux.size ());
  mesh.quad_normals_ = Eigen::Map<Eigen::Matrix<size_t, 4, Eigen::Dynamic> > (&quad_normals_aux[0][0], 4, quad_normals_aux.size ());
  mesh.quad_tex_coords_ = Eigen::Map<Eigen::Matrix<size_t, 4, Eigen::Dynamic> > (&quad_tex_coords_aux[0][0], 4, quad_tex_coords_aux.size ());

  /// Subtract 1 from all the indices, as the obj format starts counting from 1
  mesh.tri_vertices_ -= Eigen::Matrix<size_t, 3, Eigen::Dynamic>::Ones (3, mesh.tri_vertices_.cols ());
  mesh.tri_normals_ -= Eigen::Matrix<size_t, 3, Eigen::Dynamic>::Ones (3, mesh.tri_normals_.cols ());
  mesh.tri_tex_coords_ -= Eigen::Matrix<size_t, 3, Eigen::Dynamic>::Ones (3, mesh.tri_tex_coords_.cols ());
  mesh.quad_vertices_ -= Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Ones (4, mesh.quad_vertices_.cols ());
  mesh.quad_normals_ -= Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Ones (4, mesh.quad_normals_.cols ());
  mesh.quad_tex_coords_ -= Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Ones (4, mesh.quad_tex_coords_.cols ());

  return (true);
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> bool
af::Mesh<Scalar>::writeMeshOBJ (const std::string &filename,
                                const Mesh<Scalar> &mesh)
{
  std::ofstream file (filename.c_str ());

  /// Write material file first
  if (mesh.materials_.size () != 0)
  {
    std::string filename_material = filename.substr (0, filename.find_last_of (".")) + ".mtl";
    writeMaterialFile (filename_material.c_str (), mesh.materials_);

    file << "mtllib " << filename_material << "\n";
  }


  /// Write the vertices
  for (int v_i = 0; v_i < mesh.vertices_.cols (); ++v_i)
    file << "v " << mesh.vertices_(0, v_i) << " " << mesh.vertices_(1, v_i) << " " << mesh.vertices_(2, v_i) << "\n";

  /// Write the normals
  for (int n_i = 0; n_i < mesh.normals_.cols (); ++n_i)
    file << "vn " << mesh.normals_(0, n_i) << " " << mesh.normals_(1, n_i) << " " << mesh.normals_(2, n_i) << "\n";

  /// Write the texture coordinates
  for (int t_i = 0; t_i < mesh.tex_coords_.cols (); ++t_i)
    file << "vt " << mesh.tex_coords_(0, t_i) << " " << mesh.tex_coords_(1, t_i) << "\n";

  if (mesh.segments_.size () == 0)
  {
    /// Just write out the faces

    /// Triangular faces
    if (mesh.tri_normals_.cols () != 0 &&
        mesh.tri_tex_coords_.cols () != 0)
    {
      for (int f_i = 0; f_i < mesh.tri_vertices_.cols (); ++f_i)
        file << "f " << mesh.tri_vertices_(0, f_i) + 1 << "/" << mesh.tri_tex_coords_(0, f_i) + 1 << "/" << mesh.tri_normals_(0, f_i) + 1 << " " <<
                mesh.tri_vertices_(1, f_i) + 1 << "/" << mesh.tri_tex_coords_(1, f_i) + 1 << "/" << mesh.tri_normals_(1, f_i) + 1 << " " <<
                mesh.tri_vertices_(2, f_i) + 1 << "/" << mesh.tri_tex_coords_(2, f_i) + 1 << "/" << mesh.tri_normals_(2, f_i) + 1 << "\n";
    }
    else if (mesh.tri_normals_.cols () == 0 &&
             mesh.tri_tex_coords_.cols () != 0)
    {
      for (int f_i = 0; f_i < mesh.tri_vertices_.cols (); ++f_i)
        file << "f " << mesh.tri_vertices_(0, f_i) + 1 << "/" << mesh.tri_tex_coords_(0, f_i) + 1 << " " <<
                mesh.tri_vertices_(1, f_i) + 1 << "/" << mesh.tri_tex_coords_(1, f_i) + 1 << " " <<
                mesh.tri_vertices_(2, f_i) + 1 << "/" << mesh.tri_tex_coords_(2, f_i) + 1 << "\n";
    }
    else
    {
      for (int f_i = 0; f_i < mesh.tri_vertices_.cols (); ++f_i)
        file << "f " << mesh.tri_vertices_(0, f_i) + 1 << " " << mesh.tri_vertices_(1, f_i) + 1 << " " << mesh.tri_vertices_(2, f_i) + 1 << "\n";
    }

    /// Quad faces
    if (mesh.quad_normals_.cols () != 0 &&
        mesh.quad_tex_coords_.cols () != 0)
    {
      for (int f_i = 0; f_i < mesh.quad_vertices_.cols (); ++f_i)
        file << "f " << mesh.quad_vertices_(0, f_i) + 1 << "/" << mesh.quad_tex_coords_(0, f_i) + 1 << "/" <<  mesh.quad_normals_(0, f_i) + 1 << " " <<
                mesh.quad_vertices_(1, f_i) + 1 << "/" << mesh.quad_tex_coords_(1, f_i) + 1 << "/" << mesh.quad_normals_(1, f_i) + 1 << " " <<
                mesh.quad_vertices_(2, f_i) + 1 << "/" << mesh.quad_tex_coords_(2, f_i) + 1 << "/" << mesh.quad_normals_(2, f_i) + 1 << " " <<
                mesh.quad_vertices_(3, f_i) + 1 << "/" << mesh.quad_tex_coords_(3, f_i) + 1 << "/" << mesh.quad_normals_(3, f_i) + 1 << "\n";
    }
    else if (mesh.quad_normals_.cols () == 0 &&
             mesh.quad_tex_coords_.cols () != 0)
    {
      for (int f_i = 0; f_i < mesh.quad_vertices_.cols (); ++f_i)
        file << "f " << mesh.quad_vertices_(0, f_i) + 1 << "/" << mesh.quad_tex_coords_(0, f_i) + 1 << " " <<
                mesh.quad_vertices_(1, f_i) + 1 << "/" << mesh.quad_tex_coords_(1, f_i) + 1 << " " <<
                mesh.quad_vertices_(2, f_i) + 1 << "/" << mesh.quad_tex_coords_(2, f_i) + 1 << " " <<
                mesh.quad_vertices_(3, f_i) + 1 << "/" << mesh.quad_tex_coords_(3, f_i) + 1 << "\n";
    }
    else
    {
      for (int f_i = 0; f_i < mesh.quad_vertices_.cols (); ++f_i)
        file << "f " << mesh.quad_vertices_(0, f_i) + 1 << " " << mesh.quad_vertices_(1, f_i) + 1 << " " << mesh.quad_vertices_(2, f_i) + 1 << " " << mesh.quad_vertices_(3, f_i) + 1 << "\n";
    }
  }
  else
  {
    for (std::map<std::string, MeshSegment>::const_iterator s_it = mesh.segments_.begin ();
         s_it != mesh.segments_.end (); ++s_it)
    {
      file << "usemtl " << s_it->first << "\n";

      /// Triangles
      if (mesh.tri_normals_.cols () != 0 &&
          mesh.tri_tex_coords_.cols () != 0)
      {
        for (size_t tri_i = 0; tri_i < s_it->second.tri_indices.size (); ++tri_i)
        {
          size_t f_i = s_it->second.tri_indices[tri_i];
          file << "f " << mesh.tri_vertices_(0, f_i) + 1 << "/" << mesh.tri_tex_coords_(0, f_i) + 1 << "/" << mesh.tri_normals_ (0, f_i) + 1 << " " <<
                  mesh.tri_vertices_(1, f_i) + 1 << "/" << mesh.tri_tex_coords_(1, f_i) + 1 << "/" << mesh.tri_normals_ (1, f_i) + 1 << " " <<
                  mesh.tri_vertices_(2, f_i) + 1 << "/" << mesh.tri_tex_coords_(2, f_i) + 1 << "/" << mesh.tri_normals_ (2, f_i) + 1 << "\n";
        }
      }
      else if (mesh.tri_normals_.cols () == 0 &&
               mesh.tri_tex_coords_.cols () != 0)
      {
        for (size_t tri_i = 0; tri_i < s_it->second.tri_indices.size (); ++tri_i)
        {
          size_t f_i = s_it->second.tri_indices[tri_i];
          file << "f " << mesh.tri_vertices_(0, f_i) + 1 << "/" << mesh.tri_tex_coords_(0, f_i) + 1 << " " <<
                  mesh.tri_vertices_(1, f_i) + 1 << "/" << mesh.tri_tex_coords_(1, f_i) + 1 << " " <<
                  mesh.tri_vertices_(2, f_i) + 1 << "/" << mesh.tri_tex_coords_(2, f_i) + 1 << "\n";
        }
      }
      else
      {
        for (size_t tri_i = 0; tri_i < s_it->second.tri_indices.size (); ++tri_i)
        {
          size_t f_i = s_it->second.tri_indices[tri_i];
          file << "f " << mesh.tri_vertices_(0, f_i) + 1 << " " << mesh.tri_vertices_(1, f_i) + 1 << " " << mesh.tri_vertices_(2, f_i) + 1 << "\n";
        }
      }

      /// Quads
      if (mesh.quad_normals_.cols () != 0 &&
          mesh.quad_tex_coords_.cols () != 0)
      {
        for (size_t quad_i = 0; quad_i < s_it->second.quad_indices.size (); ++quad_i)
        {
          size_t f_i = s_it->second.quad_indices[quad_i];
          file << "f " << mesh.quad_vertices_(0, f_i) + 1 << "/" << mesh.quad_tex_coords_(0, f_i) + 1 << "/" <<  mesh.quad_normals_(0, f_i) + 1 << " " <<
                  mesh.quad_vertices_(1, f_i) + 1 << "/" << mesh.quad_tex_coords_(1, f_i) + 1 << "/" << mesh.quad_normals_(1, f_i) + 1 << " " <<
                  mesh.quad_vertices_(2, f_i) + 1 << "/" << mesh.quad_tex_coords_(2, f_i) + 1 << "/" << mesh.quad_normals_(2, f_i) + 1 << " " <<
                  mesh.quad_vertices_(3, f_i) + 1 << "/" << mesh.quad_tex_coords_(3, f_i) + 1 << "/" << mesh.quad_normals_(3, f_i) + 1 << "\n";
        }
      }
      else if (mesh.quad_normals_.cols () == 0 &&
               mesh.quad_tex_coords_.cols () != 0)
      {
        for (size_t quad_i = 0; quad_i < s_it->second.quad_indices.size (); ++quad_i)
        {
          size_t f_i = s_it->second.quad_indices[quad_i];
          file << "f " << mesh.quad_vertices_(0, f_i) + 1 << "/" << mesh.quad_tex_coords_(0, f_i) + 1 << " " <<
                  mesh.quad_vertices_(1, f_i) + 1 << "/" << mesh.quad_tex_coords_(1, f_i) + 1 << " " <<
                  mesh.quad_vertices_(2, f_i) + 1 << "/" << mesh.quad_tex_coords_(2, f_i) + 1 << " " <<
                  mesh.quad_vertices_(3, f_i) + 1 << "/" << mesh.quad_tex_coords_(3, f_i) + 1 << "\n";
        }
      }
      else
      {
        for (size_t quad_i = 0; quad_i < s_it->second.quad_indices.size (); ++quad_i)
        {
          size_t f_i = s_it->second.quad_indices[quad_i];
          file << "f " << mesh.quad_vertices_(0, f_i) + 1 << " " << mesh.quad_vertices_(1, f_i) + 1 << " " << mesh.quad_vertices_(2, f_i) + 1 << " " << mesh.quad_vertices_ (3, f_i) + 1 << "\n";
        }
      }
    }
  }

  return (true);
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> bool
af::Mesh<Scalar>::readMaterialFile (const std::string &filename,
                                    std::map<std::string, Material<Scalar> > &materials)
{
  std::ifstream file (filename.c_str ());
  if (!file.is_open ())
  {
    PCL_ERROR ("Error opening material file %s\n", filename.c_str ());
    return (false);
  }

  std::string line, line_type;
  Material<Scalar> *current_material = 0;
  while (1)
  {
    std::getline (file, line);
    if (!file.good ())
      break;

    std::stringstream line_ss (line);
    line_ss >> line_type;
    boost::algorithm::to_lower (line_type);

    if (line_type == "newmtl")
    {
      if (current_material != 0)
        materials.insert (std::make_pair (current_material->name, *current_material));

      current_material = new Material<Scalar> ();
      line_ss >> current_material->name;
    }
    else if (current_material != 0)
    {
      if (line_type == "ka")
        line_ss >> current_material->color_ambient[0] >> current_material->color_ambient[1] >>  current_material->color_ambient[2];
      else if (line_type == "kd")
        line_ss >> current_material->color_diffuse[0] >> current_material->color_diffuse[1] >> current_material->color_diffuse[2];
      else if (line_type == "ks")
        line_ss >> current_material->color_specular[0] >> current_material->color_specular[1] >> current_material->color_specular[2];
      else if (line_type == "ns")
        line_ss >> current_material->specular_coeff;
      else if (line_type == "d" || line_type == "Tr")
        line_ss >> current_material->transparency;
      else if (line_type == "map_kd")
        line_ss >> current_material->texture_diffuse;
      else if (line_type == "map_ks")
        line_ss >> current_material->texture_specular;
      else if (line_type == "map_bump" ||
               line_type == "bump")
        line_ss >> current_material->texture_normal;
    }
  }

  if (current_material != 0)
    materials.insert (std::make_pair (current_material->name, *current_material));

  file.close ();

  return (true);
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> bool
af::Mesh<Scalar>::writeMaterialFile (const std::string &filename,
                                     const std::map<std::string, Material<Scalar> > &materials)
{
  std::ofstream file (filename.c_str ());
  for (typename std::map<std::string, Material<Scalar> >::const_iterator m_it = materials.begin ();
       m_it != materials.end (); ++m_it)
  {
    file << "newmtl " << m_it->second.name << "\n";
    file << "Ka " << m_it->second.color_ambient[0] << " " << m_it->second.color_ambient[1] << " " << m_it->second.color_ambient[2] << "\n";
    file << "Kd " << m_it->second.color_diffuse[0] << " " << m_it->second.color_diffuse[1] << " " << m_it->second.color_diffuse[2] << "\n";
    file << "Ks " << m_it->second.color_specular[0] << " " << m_it->second.color_specular[1] << " " << m_it->second.color_specular[2] << "\n";
    file << "Ns " << m_it->second.specular_coeff << "\n";
    file << "Tr " << m_it->second.transparency << "\n";

    if (m_it->second.texture_diffuse != "")
      file << "map_Kd " << m_it->second.texture_diffuse << "\n";
    if (m_it->second.texture_specular != "")
      file << "map_Ks " << m_it->second.texture_specular << "\n";
    if (m_it->second.texture_normal != "")
      file << "map_bump " << m_it->second.texture_normal << "\n";
  }

  file.close ();

  return (true);
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> void
af::Mesh<Scalar>::printInfo () const
{
  PCL_INFO ("Mesh stats:\n   #vertices: %ld\n   #normals: %ld\n   #tex_coords: %ld\n   #tri_vertices: %ld\n   #tri_normals: %ld\n   #tri_tex_coords: %ld\n   #quad_vertices: %ld\n   #quad_normals: %ld\n   #quad_tex_coords: %ld\n   #segments: %ld\n   #materials: %ld\n",
            vertices_.cols (), normals_.cols (), tex_coords_.cols (),
            tri_vertices_.cols (), tri_normals_.cols (), tri_tex_coords_.cols (),
            quad_vertices_.cols (), quad_normals_.cols (), quad_tex_coords_.cols (),
            segments_.size (), materials_.size ());

  PCL_INFO ("Segments: %ld of sizes: ", segments_.size ());
  for (std::map<std::string, MeshSegment>::const_iterator s_it = segments_.begin ();
       s_it != segments_.end (); ++s_it)
    PCL_INFO ("%ld ", s_it->second.size ());
  PCL_INFO ("\n");

  PCL_INFO ("Materials: %ld\n", materials_.size ());
  for (typename std::map<std::string, Material<Scalar> >::const_iterator m_it = materials_.begin ();
       m_it != materials_.end (); ++m_it)
    PCL_INFO (" Material name: %s\n   ambient: %f %f %f\n   diffuse: %f %f %f\n   specular: %f %f %f\n   specular_coeff: %f\n   transparency: %f\n   diffuse texture: %s\n   specular texture: %s\n   bump map: %s\n\n",
              m_it->second.name.c_str (),
              m_it->second.color_ambient[0], m_it->second.color_ambient[1], m_it->second.color_ambient[2],
        m_it->second.color_diffuse[0], m_it->second.color_diffuse[1], m_it->second.color_diffuse[2],
        m_it->second.color_specular[0], m_it->second.color_specular[1], m_it->second.color_specular[2],
        m_it->second.specular_coeff,
        m_it->second.transparency,
        m_it->second.texture_diffuse.c_str (), m_it->second.texture_specular.c_str (), m_it->second.texture_normal.c_str ());
}


////////////////////////////////////////////////////////////////////////////////
template <typename Scalar> void
af::Mesh<Scalar>::convertToTriangleMesh (const Mesh<Scalar> &mesh, Mesh<Scalar> &result)
{
  /// Copy all the triangles to the result
  result = mesh;
  result.tri_tex_coords_.resize (3, mesh.tri_vertices_.cols () + 2 * mesh.quad_vertices_.cols ());
  result.tri_vertices_.resize (3, mesh.tri_vertices_.cols () + 2 * mesh.quad_vertices_.cols ());

  result.tri_vertices_.block (0, 0, mesh.tri_vertices_.rows (), mesh.tri_vertices_.cols ()) = mesh.tri_vertices_;
  result.tri_tex_coords_.block (0, 0, mesh.tri_tex_coords_.rows (), mesh.tri_tex_coords_.cols ()) = mesh.tri_tex_coords_;

  result.quad_normals_.resize (4, 0);
  result.quad_tex_coords_.resize (4, 0);
  result.quad_vertices_.resize (4, 0);

  /// HACK put all the triangles in a single segment
  result.segments_.clear ();
  MeshSegment segment;
  for (int i = 0; i < result.tri_vertices_.cols (); ++i)
    segment.tri_indices.push_back (i);
  result.segments_.insert (std::make_pair ("segment", segment));


  int tri_index = mesh.tri_vertices_.cols ();


  /// For each quad
  for (int quad_i = 0; quad_i < mesh.quad_vertices_.cols (); ++quad_i)
  {
    /// Triangle 1
    Eigen::Matrix<Scalar, 3, 1> a = mesh.vertices_.col (mesh.quad_vertices_ (0, quad_i));
    Eigen::Matrix<Scalar, 3, 1> b = mesh.vertices_.col (mesh.quad_vertices_ (1, quad_i));
    Eigen::Matrix<Scalar, 3, 1> c = mesh.vertices_.col (mesh.quad_vertices_ (2, quad_i));
    Eigen::Matrix<Scalar, 3, 1> d = mesh.vertices_.col (mesh.quad_vertices_ (3, quad_i));


    /// Edge b-d
    double energy_1 = fabs (60. - (d-a).normalized ().dot ((b-a).normalized ())) +
        fabs (60. - (a-b).normalized ().dot ((d-b).normalized ())) +
        fabs (60. - (a-d).normalized ().dot ((b-d).normalized ())) +
        fabs (60. - (d-c).normalized ().dot ((b-c).normalized ())) +
        fabs (60. - (c-b).normalized ().dot ((d-b).normalized ())) +
        fabs (60. - (c-d).normalized ().dot ((b-d).normalized ()));

    /// Edge a-c
    double energy_2 = fabs (60. - (b-a).normalized ().dot ((c-a).normalized ())) +
        fabs (60. - (a-b).normalized ().dot ((c-b).normalized ())) +
        fabs (60. - (a-c).normalized ().dot ((b-c).normalized ())) +
        fabs (60. - (d-a).normalized ().dot ((c-a).normalized ())) +
        fabs (60. - (a-c).normalized ().dot ((d-c).normalized ())) +
        fabs (60. - (a-d).normalized ().dot ((c-d).normalized ()));

    if (energy_1 < energy_2)
    {
      result.tri_vertices_ (0, tri_index) = mesh.quad_vertices_ (0, quad_i);
      result.tri_vertices_ (1, tri_index) = mesh.quad_vertices_ (1, quad_i);
      result.tri_vertices_ (2, tri_index) = mesh.quad_vertices_ (3, quad_i);

      if (mesh.quad_tex_coords_.cols () > 0)
      {
        result.tri_tex_coords_ (0, tri_index) = mesh.quad_tex_coords_ (0, quad_i);
        result.tri_tex_coords_ (1, tri_index) = mesh.quad_tex_coords_ (1, quad_i);
        result.tri_tex_coords_ (2, tri_index) = mesh.quad_tex_coords_ (3, quad_i);
      }

//      result.tri_normals_ (0, tri_index) = mesh.quad_normals_ (0, quad_i);
//      result.tri_normals_ (1, tri_index) = mesh.quad_normals_ (1, quad_i);
//      result.tri_normals_ (2, tri_index) = mesh.quad_normals_ (3, quad_i);


      tri_index ++;
      result.tri_vertices_ (0, tri_index) = mesh.quad_vertices_ (1, quad_i);
      result.tri_vertices_ (1, tri_index) = mesh.quad_vertices_ (2, quad_i);
      result.tri_vertices_ (2, tri_index) = mesh.quad_vertices_ (3, quad_i);
  
      if (mesh.quad_tex_coords_.cols () > 0)
      {
        result.tri_tex_coords_ (0, tri_index) = mesh.quad_tex_coords_ (1, quad_i);
        result.tri_tex_coords_ (1, tri_index) = mesh.quad_tex_coords_ (2, quad_i);
        result.tri_tex_coords_ (2, tri_index) = mesh.quad_tex_coords_ (3, quad_i);
      }

//      result.tri_normals_ (0, tri_index) = mesh.quad_normals_ (1, quad_i);
//      result.tri_normals_ (1, tri_index) = mesh.quad_normals_ (2, quad_i);
//      result.tri_normals_ (2, tri_index) = mesh.quad_normals_ (3, quad_i);
      tri_index ++;
    }
    else
    {
      result.tri_vertices_ (0, tri_index) = mesh.quad_vertices_ (0, quad_i);
      result.tri_vertices_ (1, tri_index) = mesh.quad_vertices_ (1, quad_i);
      result.tri_vertices_ (2, tri_index) = mesh.quad_vertices_ (2, quad_i);

      if (mesh.quad_tex_coords_.cols () > 0)
      {
        result.tri_tex_coords_ (0, tri_index) = mesh.quad_tex_coords_ (0, quad_i);
        result.tri_tex_coords_ (1, tri_index) = mesh.quad_tex_coords_ (1, quad_i);
        result.tri_tex_coords_ (2, tri_index) = mesh.quad_tex_coords_ (2, quad_i);
      }

//      result.tri_normals_ (0, tri_index) = mesh.quad_normals_ (0, quad_i);
//      result.tri_normals_ (1, tri_index) = mesh.quad_normals_ (1, quad_i);
//      result.tri_normals_ (2, tri_index) = mesh.quad_normals_ (2, quad_i);


      tri_index ++;
      result.tri_vertices_ (0, tri_index) = mesh.quad_vertices_ (2, quad_i);
      result.tri_vertices_ (1, tri_index) = mesh.quad_vertices_ (3, quad_i);
      result.tri_vertices_ (2, tri_index) = mesh.quad_vertices_ (0, quad_i);
      
      if (mesh.quad_tex_coords_.cols () > 0)
      {
        result.tri_tex_coords_ (0, tri_index) = mesh.quad_tex_coords_ (2, quad_i);
        result.tri_tex_coords_ (1, tri_index) = mesh.quad_tex_coords_ (3, quad_i);
        result.tri_tex_coords_ (2, tri_index) = mesh.quad_tex_coords_ (0, quad_i);
      }
      
//      result.tri_normals_ (0, tri_index) = mesh.quad_normals_ (2, quad_i);
//      result.tri_normals_ (1, tri_index) = mesh.quad_normals_ (3, quad_i);
//      result.tri_normals_ (2, tri_index) = mesh.quad_normals_ (0, quad_i);
      tri_index ++;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////
template<typename Scalar> void
af::Mesh<Scalar>::computeTexelGeometry (const Mesh<Scalar> &mesh,
                                        std::vector<af::Texel> &texels,
                                        int texture_size)
{
  /// For each triangle
  cv::Mat texture_colors = cv::Mat::zeros (texture_size, texture_size, CV_8UC3);
  for (int tri_i = 0; tri_i < mesh.tri_tex_coords_.cols (); ++tri_i)
  {
    /// Get the triangle texture coordinates
    Eigen::Matrix<Scalar, 2, 1> p0 (mesh.tex_coords_.col (mesh.tri_tex_coords_ (0, tri_i))),
                                p1 (mesh.tex_coords_.col (mesh.tri_tex_coords_ (1, tri_i))),
                                p2 (mesh.tex_coords_.col (mesh.tri_tex_coords_ (2, tri_i)));

    cv::Point tri_points[1][3];
    tri_points[0][0] = cv::Point (p0[0] * texture_size, p0[1] * texture_size);
    tri_points[0][1] = cv::Point (p1[0] * texture_size, p1[1] * texture_size);
    tri_points[0][2] = cv::Point (p2[0] * texture_size, p2[1] * texture_size);
    const cv::Point* ppt[1] = { tri_points[0] };

    int vertex_counter[1];
    vertex_counter[0] = 3;

    unsigned char r, g, b;
    size_t temp = tri_i + 1; /// to avoid confusion with the black background
    r = temp / (256 * 256);
    g = (temp % (256 * 256)) / 256;
    b = temp % 256;

    cv::fillPoly (texture_colors, ppt, vertex_counter, 1, cv::Scalar (r, g, b));
  }


  for (int x = 0; x < texture_size; ++x)
    for (int y = 0; y < texture_size; ++y)
    {
      cv::Vec3b pix = texture_colors.at<cv::Vec3b> (cv::Point (x, y));
      int tri_id = pix[0] * 256 * 256 + pix[1] * 256 + pix[2];
      if (tri_id == 0)
        continue;

      tri_id --;

      Eigen::Matrix<Scalar, 3, 3> U;
      U (0, 0) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (0, tri_id))[0];
      U (1, 0) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (0, tri_id))[1];
      U (2, 0) = 1.;
      U (0, 1) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (1, tri_id))[0];
      U (1, 1) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (1, tri_id))[1];
      U (2, 1) = 1.;
      U (0, 2) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (2, tri_id))[0];
      U (1, 2) = mesh.tex_coords_.col (mesh.tri_tex_coords_ (2, tri_id))[1];
      U (2, 2) = 1.;

      Eigen::Matrix<Scalar, 3, 1> A = U.inverse () * Eigen::Matrix<Scalar, 3, 1> (static_cast<double> (x) / static_cast<double> (texture_size),
                                                                                  static_cast<double> (y) / static_cast<double> (texture_size),
                                                                                  1.);

      texels.push_back (af::Texel (x, texture_size - 1 - y, tri_id, A));
    }

//  cv::imshow ("texture colors", texture_colors);
//  cv::waitKey ();
}


////////////////////////////////////////////////////////////////////////////////
template<typename Scalar> void
af::Mesh<Scalar>::concatenateMeshes (const Mesh<Scalar> &mesh_a, const Mesh<Scalar> &mesh_b, Mesh<Scalar> &result)
{
  result.vertices_.resize (3, mesh_a.vertices_.cols () + mesh_b.vertices_.cols ());
  result.vertices_.block (0, 0, 3, mesh_a.vertices_.cols ()) = mesh_a.vertices_;
  result.vertices_.block (0, mesh_a.vertices_.cols (), 3, mesh_b.vertices_.cols ()) = mesh_b.vertices_;

  result.tex_coords_.resize (2, mesh_a.tex_coords_.cols () + mesh_b.tex_coords_.cols ());
  result.tex_coords_.block (0, 0, 2, mesh_a.tex_coords_.cols ()) = mesh_a.tex_coords_;
  result.tex_coords_.block (0, mesh_a.tex_coords_.cols (), 2, mesh_b.tex_coords_.cols ()) = mesh_b.tex_coords_;

  /// Fix the faces
  result.tri_vertices_.resize (3, mesh_a.tri_vertices_.cols () + mesh_b.tri_vertices_.cols ());
  result.tri_vertices_.block (0, 0, 3, mesh_a.tri_vertices_.cols ()) = mesh_a.tri_vertices_;
  result.tri_vertices_.block (0, mesh_a.tri_vertices_.cols (), 3, mesh_b.tri_vertices_.cols ()) =
      mesh_b.tri_vertices_ + Eigen::Matrix<size_t, 3, Eigen::Dynamic>::Ones (3, mesh_b.tri_vertices_.cols ()) * mesh_a.vertices_.cols ();

  result.quad_vertices_.resize (4, mesh_a.quad_vertices_.cols () + mesh_b.quad_vertices_.cols ());
  result.quad_vertices_.block (0, 0, 4, mesh_a.quad_vertices_.cols ()) = mesh_a.quad_vertices_;
  result.quad_vertices_.block (0, mesh_a.quad_vertices_.cols (), 4, mesh_b.quad_vertices_.cols ()) =
      mesh_b.quad_vertices_ + Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Ones (4, mesh_b.quad_vertices_.cols ()) * mesh_a.vertices_.cols ();

  result.tri_tex_coords_.resize (3, mesh_a.tri_tex_coords_.cols () + mesh_b.tri_tex_coords_.cols ());
  result.tri_tex_coords_.block (0, 0, 3, mesh_a.tri_tex_coords_.cols ()) = mesh_a.tri_tex_coords_;
  result.tri_tex_coords_.block (0, mesh_a.tri_tex_coords_.cols (), 3, mesh_b.tri_tex_coords_.cols ()) =
      mesh_b.tri_tex_coords_ + Eigen::Matrix<size_t, 3, Eigen::Dynamic>::Ones (3, mesh_b.tri_tex_coords_.cols ()) * mesh_a.tex_coords_.cols ();

  result.quad_tex_coords_.resize (4, mesh_a.quad_tex_coords_.cols () + mesh_b.quad_tex_coords_.cols ());
  result.quad_tex_coords_.block (0, 0, 4, mesh_a.quad_tex_coords_.cols ()) = mesh_a.quad_tex_coords_;
  result.quad_tex_coords_.block (0, mesh_a.quad_tex_coords_.cols (), 4, mesh_b.quad_tex_coords_.cols ()) =
      mesh_b.quad_tex_coords_ + Eigen::Matrix<size_t, 4, Eigen::Dynamic>::Ones (4, mesh_b.quad_tex_coords_.cols ()) * mesh_a.tex_coords_.cols ();
}
