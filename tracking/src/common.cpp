#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/optional.hpp>

#include "common.h"

#include <fstream>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>



void
af::CameraParams::print () const
{
  PCL_INFO ("Camera parameters:\n   image_path: %s\n   height: %d\n   width: %d\n   fx: %f\n   fy: %f\n   cx: %f\n   cy: %f\n",
            image_path.c_str (), height, width, fx, fy, cx, cy);
  PCL_INFO ("   skew: %f\n   k1: %f\n   k2: %f\n   k3: %f\n   p1: %f\n   p2: %f\n",
            skew, k1, k2, k3, p1, p2);
  std::cout << "   pose:\n" << pose.matrix () << std::endl;
}

void
af::clamp (double &val, const double min, const double max)
{
  if (val < min)
    val = min;
  else if (val > max)
    val = max;
}


void
af::getFilesFromDirectory (const std::string path_dir,
                           const std::string extension,
                           std::vector<std::string> &files)
{
  if (path_dir != "" && boost::filesystem::exists (path_dir))
  {
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr (path_dir); itr != end_itr; ++itr)
      if (!is_directory (itr->status ()) &&
          boost::algorithm::to_upper_copy (boost::filesystem::extension (itr->path ())) == boost::algorithm::to_upper_copy (extension))
      {
        files.push_back (itr->path ().string ());
      }
  }
  sort (files.begin (), files.end ());
}


void
af::getFilesFromDirectory (const std::string path_dir,
                           std::vector<std::string> &files)
{
  if (path_dir != "" && boost::filesystem::exists (path_dir))
  {
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr (path_dir); itr != end_itr; ++itr)
      if (!is_directory (itr->status ()))
        files.push_back (itr->path ().string ());
  }
  sort (files.begin (), files.end ());
}


void
af::getDirectoriesFromDirectory (const std::string path_dir,
                                 std::vector<std::string> &dirs)
{
  if (path_dir != "" && boost::filesystem::exists (path_dir))
  {
    boost::filesystem::directory_iterator end_itr;
    for (boost::filesystem::directory_iterator itr (path_dir); itr != end_itr; ++itr)
      if (is_directory (itr->status ()))
        dirs.push_back (itr->path ().string ());
  }
  sort (dirs.begin (), dirs.end ());
}


std::string
af::getBaseName (const std::string &filename)
{
  size_t slash_pos = filename.find_last_of ('/');
  size_t dot_pos = filename.find_last_of (".");
  std::string base_name = filename.substr (slash_pos + 1, dot_pos - slash_pos - 1);

  return (base_name);
}

std::string
af::getExtension (const std::string &filename)
{
  size_t dot_pos = filename.find_last_of (".");
  std::string extension = filename.substr (dot_pos + 1, filename.length () - dot_pos);

  return (boost::algorithm::to_upper_copy (extension));
}


bool
af::isMatrixHealthy (const Eigen::MatrixXd &matrix)
{
  for (int i = 0; i < matrix.rows (); ++i)
    for (int j = 0; j < matrix.cols (); ++j)
      if (!pcl_isfinite (matrix (i, j)))
        return (false);

  return (true);
}



bool
af::readIndicesFile (const std::string &filename, std::vector<size_t> &indices)
{
  std::ifstream file (filename.c_str ());
  if (!file.is_open ())
  {
    PCL_ERROR ("Error reading indices file %s.\n", filename.c_str ());
    return (false);
  }
  size_t points_size;
  file >> points_size;
  for (size_t i = 0; i < points_size; ++i)
  {
    size_t v;
    file >> v;
    indices.push_back (v);
  }
  file.close ();
  return (true);
}


bool
af::writeIndicesFile (const std::string &filename, const std::vector<size_t> &indices)
{
  std::ofstream file (filename.c_str ());
  if (!file.is_open ())
  {
    PCL_ERROR ("Error writing indices file %s.\n", filename.c_str ());
    return (false);
  }

  file << indices.size () << "\n";
  for (size_t i = 0; i < indices.size (); ++i)
    file << indices[i] << " ";

  file.close ();
  return (true);
}


af::Mesh<float>::Ptr
af::generateCylinder (const double radius, const double length,
                      const int num_circle_samples, const int num_length_samples)
{
  /// First generate a unit length, unit radius cylinder pointing upwards
  af::Mesh<float>::Ptr mesh (new af::Mesh<float> ());
  mesh->vertices_.resize (3, (num_length_samples + 1) * num_circle_samples);
  mesh->tri_vertices_.resize (3, 2 * (num_length_samples + 1 - 1) * num_circle_samples + 2 * (num_circle_samples - 2));
  int tri_i = 0;
  for (int l_i = 0; l_i <= num_length_samples; ++l_i)
  {
    for (int c_i = 0; c_i < num_circle_samples; ++c_i)
    {
      /// Generate the point
      Eigen::Vector3f p;
      p (0) = radius * cos (2. * M_PI / static_cast<double> (num_circle_samples) * static_cast<double> (c_i));
      p (1) = radius * sin (2. * M_PI / static_cast<double> (num_circle_samples) * static_cast<double> (c_i));
      p (2) = length * static_cast<double> (l_i) / static_cast<double> (num_length_samples);

      mesh->vertices_.col (l_i * num_circle_samples + c_i) = p;

      /// Form two triangles
      if (l_i != 0)
      {
        mesh->tri_vertices_ (0, tri_i) = l_i * num_circle_samples + c_i;
        mesh->tri_vertices_ (1, tri_i) = (l_i - 1) * num_circle_samples + c_i;
        mesh->tri_vertices_ (2, tri_i) = l_i * num_circle_samples + (c_i + 1) % num_circle_samples;
        tri_i ++;

        mesh->tri_vertices_ (0, tri_i) = l_i * num_circle_samples + c_i;
        mesh->tri_vertices_ (2, tri_i) = (l_i - 1) * num_circle_samples + c_i;
        mesh->tri_vertices_ (1, tri_i) = (l_i - 1) * num_circle_samples + (c_i - 1 + num_circle_samples) % num_circle_samples;
        tri_i ++;
      }
      if (l_i == 0 && c_i != 0 && c_i != num_circle_samples - 1)
      {
        mesh->tri_vertices_ (0, tri_i) = 0;
        mesh->tri_vertices_ (2, tri_i) = c_i;
        mesh->tri_vertices_ (1, tri_i) = (c_i + 1) % num_circle_samples;
        tri_i ++;
      }
      if (l_i == num_length_samples && c_i != 0 && c_i != num_circle_samples - 1)
      {
        mesh->tri_vertices_ (0, tri_i) = l_i * num_circle_samples + 0;
        mesh->tri_vertices_ (1, tri_i) = l_i * num_circle_samples + c_i;
        mesh->tri_vertices_ (2, tri_i) = l_i * num_circle_samples + (c_i + 1) % num_circle_samples;
        tri_i ++;
      }
    }
  }

  return (mesh);
}


af::Mesh<float>::Ptr
af::generateCylinder (Eigen::Vector3f &orig, Eigen::Vector3f &tip,
                      const int num_circle_samples, const int num_length_samples,
                      const double radius)
{
  double length = (orig - tip).norm ();
  af::Mesh<float>::Ptr mesh = generateCylinder (radius, length, num_circle_samples, num_length_samples);

  /// Rotate the cylinder
  float angle = acos (Eigen::Vector3f::UnitZ ().dot ((tip - orig).normalized ()));
  Eigen::Vector3f axis = (Eigen::Vector3f::UnitZ ().cross ((tip - orig).normalized ())).normalized ();
  Eigen::Matrix3f rot = Eigen::AngleAxisf (angle, axis).matrix ();
  std::cerr << "rotation matrix:\n" << rot << std::endl;

  for (size_t v_i = 0; v_i < mesh->vertices_.cols (); ++v_i)
    mesh->vertices_.col (v_i) = rot * mesh->vertices_.col (v_i) + orig;

  return (mesh);
}


af::Mesh<float>::Ptr
af::generateSphere (const Eigen::Vector3f &center,
                    const float radius,
                    const int num_samples)
{
  af::Mesh<float>::Ptr mesh (new af::Mesh<float> ());
  mesh->vertices_.resize (3, (num_samples + 1) * (num_samples + 1));
  mesh->tri_vertices_.resize (3, 2 * (num_samples + 1) * (num_samples - 1 + 1));

//  mesh->vertices_ = Eigen::Matrix3Xf::Zero (3, num_samples * (num_samples + 1));
  int tri_i = 0;
  for (int phi_i = 0; phi_i <= num_samples; ++phi_i)
    for (int theta_i = 0; theta_i <= num_samples; ++theta_i)
    {
      Eigen::Vector3f p;
      p (0) = radius * cos (2. * M_PI * static_cast<double> (theta_i) / static_cast<double> (num_samples)) *
                    sin (M_PI * static_cast<double> (phi_i) / static_cast<double> (num_samples));
      p (1) = radius * sin (2. * M_PI * static_cast<double> (theta_i) / static_cast<double> (num_samples)) *
                    sin (M_PI * static_cast<double> (phi_i) / static_cast<double> (num_samples));
      p (2) = radius * cos (M_PI * static_cast<double> (phi_i) / static_cast<double> (num_samples));
      mesh->vertices_.col (phi_i * (num_samples + 1) + theta_i) = p + center;

      /// Form two triangles
      if (phi_i != 0)
      {
        mesh->tri_vertices_ (0, tri_i) = phi_i * (num_samples + 1) + theta_i;
        mesh->tri_vertices_ (2, tri_i) = (phi_i - 1) * (num_samples + 1) + theta_i;
        mesh->tri_vertices_ (1, tri_i) = phi_i * (num_samples + 1) + (theta_i + 1) % num_samples;
        tri_i ++;

        mesh->tri_vertices_ (0, tri_i) = phi_i * (num_samples + 1) + theta_i;
        mesh->tri_vertices_ (1, tri_i) = (phi_i - 1) * (num_samples + 1) + theta_i;
        mesh->tri_vertices_ (2, tri_i) = (phi_i - 1) * (num_samples + 1) + (theta_i - 1 + num_samples) % num_samples;
        tri_i ++;
      }
    }

  return (mesh);
}



void
af::displayBlendshapeWeights (const Eigen::VectorXf &values,
                              cv::Mat &img,
                              const int bar_height, const int bar_width)
{
  int num_bs = values.rows ();
  PCL_ERROR ("Num bs: %d\n", num_bs);

  img = cv::Mat::zeros (bar_height * 2, bar_width * num_bs, CV_8UC3);
  img.setTo (cv::Scalar (255, 255, 255));

  cv::Point origin (0, bar_height);

  PCL_ERROR ("Origin: %d %d\n", origin.x, origin.y);

//  cv::rectangle (img, cv::Rect (origin, origin + cv::Point (num_bs * bar_width, -10 - 50)), cv::Scalar (50, 50, 50), CV_FILLED);
//  cv::rectangle (img, cv::Rect (origin, origin + scale * cv::Point (14 * 4, -10 - 50)), cv::Scalar (15, 230, 230), CV_FILLED);
//  cv::rectangle (img, cv::Rect (origin, origin + scale * cv::Point (num_bs * 4, -10 - 50)), cv::Scalar (180, 180, 180));
  for (size_t i = 0; i < num_bs; ++i)
  {
    cv::Point location = origin + cv::Point (bar_width * i, 0);
//    char str[4];
//    sprintf (str, "%zu", i);
//    cv::putText (img, str, location, cv::FONT_HERSHEY_PLAIN, 0.35, cv::Scalar (255, 255, 255));

    if (fabs (values (i)) < 0.1)
    {
#ifdef OPENCV_3_0
		cv::rectangle(img, cv::Rect(location + cv::Point(0, -10),
			location + cv::Point(bar_width, -10 - 5)), cv::Scalar(0, 0, 255), -1);
#else
		cv::rectangle(img, cv::Rect(location + cv::Point(0, -10),
			location + cv::Point(bar_width, -10 - 5)), cv::Scalar(0, 0, 255), CV_FILLED);
#endif
      cv::rectangle (img, cv::Rect (location + cv::Point (0, -10),
                                    location + cv::Point (bar_width, -10 - 5)), cv::Scalar (150, 150, 150));
    }
    else
    {
#ifdef OPENCV_3_0
		cv::rectangle(img, cv::Rect(location + cv::Point(0, -10),
			location + cv::Point(bar_width, -10 - values(i) * bar_height)), cv::Scalar(255, 0, 0), 2);
#else
		cv::rectangle(img, cv::Rect(location + cv::Point(0, -10),
			location + cv::Point(bar_width, -10 - values(i) * bar_height)), cv::Scalar(255, 0, 0), CV_FILLED);
#endif
      cv::rectangle (img, cv::Rect (location + cv::Point (0, -10),
                                    location + cv::Point (bar_width, -10 - values (i) * bar_height)), cv::Scalar (150, 150, 150));
    }
  }
}


void
af::saveTextForWordCloud (const std::string &filename,
                          const af::BlendshapeModel<float>::ConstPtr bs_model,
                          const Eigen::VectorXf &bs_weights)
{
  std::ofstream file (filename);
  for (size_t bs_i = 0; bs_i < bs_model->num_exps_; ++bs_i)
  {
    int count = fabs (bs_weights (bs_i)) * 40.;
    for (size_t i = 0; i < count; ++i)
      file << bs_model->exp_names_[bs_i] << std::endl;
  }
  file.close ();
}



bool
af::saveMatrix (const std::string filename, const Eigen::MatrixXf &matrix)
{
  std::ofstream file (filename.c_str (), std::ios::binary | std::ios::out);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error writing matrix to file %s.", filename.c_str ());
    return (false);
  }

  size_t num_rows = matrix.rows ();
  size_t num_cols = matrix.cols ();
  file.write (reinterpret_cast<char*> (&num_rows), sizeof (num_rows));
  file.write (reinterpret_cast<char*> (&num_cols), sizeof (num_cols));

  for (size_t i = 0; i < num_rows; ++i)
    for (size_t j = 0; j < num_cols; ++j)
    {
      float val = matrix (i, j);
      file.write (reinterpret_cast<char*> (&val), sizeof (val));
    }
  file.close ();

  return (true);
}


bool
af::loadMatrix (const std::string filename, Eigen::MatrixXf &matrix)
{
  std::ifstream file (filename.c_str (), std::ios::binary | std::ios::in);
  if (!file.is_open ())
  {
    PCL_ERROR ("Error reading matrix from file %s.\n", filename.c_str ());
    return (false);
  }

  size_t num_rows, num_cols;
  file.read (reinterpret_cast<char*> (&num_rows), sizeof (num_rows));
  file.read (reinterpret_cast<char*> (&num_cols), sizeof (num_cols));

  printf ("Reading matrix of size %zu %zu from %s\n", num_rows, num_cols, filename.c_str ());

  matrix.resize (num_rows, num_cols);
  for (size_t i = 0; i < num_rows; ++i)
    for (size_t j = 0; j < num_cols; ++j)
    {
      float val;
      file.read (reinterpret_cast<char*> (&val), sizeof (val));
      matrix (i, j) = val;
    }
  file.close ();

  return (true);
}
