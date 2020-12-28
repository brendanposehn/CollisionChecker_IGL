#include <igl/read_triangle_mesh.h> 
#include <igl/get_seconds.h>
#include <igl/material_colors.h>
#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/swept_volume.h>
#include <igl/copyleft/cgal/mesh_boolean.h> //bb added
#include <igl/copyleft/cgal/intersect_other.h> //bb added
#include <igl/opengl/glfw/Viewer.h>
#include <igl/readPLY.h> //bb added 
#include <igl/PI.h>
#include <Eigen/Core>
#include <Eigen/Dense>

#include <iostream>

#include "collision_helpers.h"

#include <set>

int main(int argc, char * argv[])
{
  using namespace std;
  using namespace igl;
  
  Eigen::MatrixXi F,SF,Fp,FC;
  Eigen::MatrixXd V,SV,Vp,VC;
  bool show_swept_volume = false;
  // Define a rigid motion
  const auto & transform_lame = [](const double t)->Eigen::Affine3d
  {
    Eigen::Affine3d T = Eigen::Affine3d::Identity();
    T.rotate(Eigen::AngleAxisd(t*2.*igl::PI,Eigen::Vector3d(0,1,0)));
    T.translate(Eigen::Vector3d(0,0.125*cos(2.*igl::PI*t),0));
    return T;
  };

  const auto & transform_rot_x = [](const double t)->Eigen::Affine3d
  {
    Eigen::Affine3d T = Eigen::Affine3d::Identity();
    T.rotate(Eigen::AngleAxisd(t*2.*igl::PI,Eigen::Vector3d(1,0,0)));
    return T;
  };

  // Read in inputs as double precision floating point meshes
  // read_triangle_mesh("/home/brend/BCCancer/CollisionChecker_Tester/bunny.off",V,F);
  //read_triangle_mesh("/home/brend/BCCancer/meshes/patient_and_carbon.ply", V, F); //should work with .ply files?
  // read_triangle_mesh("/home/brend/BCCancer/meshes/1e3pts_just_linac.ply",V,F);
  read_triangle_mesh("/home/brend/BCCancer/meshes/1e4pts_just_linac.ply",V,F);


  const double center_x = -1945.849976;
  const double center_y = 788.838989;
  const double center_z = 1120.229980; //as defined by the marker object

  Eigen::MatrixXd TV;
  std::cout << "translating mesh" << std::endl;
  TV = translate_mesh(V, -center_x, -center_y, -center_z);
  std::cout << "done translating" << std::endl;

  //these are the flat pieces to either side of the slightly tilted part^
  Eigen::Vector3d pt_on_flat_1;
  pt_on_flat_1 << -1041.829956, 708.521973, 1798.189941;
  Eigen::Vector3d pt_on_flat_2;
  pt_on_flat_2 << -1362.209961, 1425.619995, 1881.880005;
  Eigen::Vector3d pt_on_flat_3;
  pt_on_flat_3 << -1332.599976, 1352.079956, 2158.459961;
  double theta_z = std::atan((pt_on_flat_1(0) - pt_on_flat_2(0))/(pt_on_flat_1(1) - pt_on_flat_2(1)));


  // all on the main circular part that should all be at the same z
  Eigen::Vector3d pt_on_disc_1;
  pt_on_disc_1 << -2096.889893, 1140.979980, 715.638977;
  Eigen::Vector3d pt_on_disc_2;
  pt_on_disc_2 << -1924.589966, 532.711975, 629.797974;
  Eigen::Vector3d pt_on_disc_3;
  pt_on_disc_3 << -1678.760010, 1000.619995, 676.361023;

  Eigen::MatrixXd RV_1, RV_2;

  Eigen::Vector3d n_approx_z;
  n_approx_z = (pt_on_disc_2 - pt_on_disc_1).cross(pt_on_disc_3 - pt_on_disc_1);
  Eigen::Vector3d n_z;
  n_z << 0, 0, 1;
  Eigen::Vector3d rot_to_z_axis;
  rot_to_z_axis = n_approx_z.cross(n_z);
  rot_to_z_axis.normalize();

  double theta_rot_z = std::acos(n_z.dot(n_approx_z)/(n_approx_z.norm()));

  Eigen::Matrix3d rot_mat_z;
  rot_mat_z = Eigen::AngleAxisd(theta_rot_z, rot_to_z_axis); //does the rot_to_z_axis need to be normalized?

  RV_1 = rotate_mesh(TV, rot_mat_z);

  Eigen::Vector3d n_approx_x;
  n_approx_x = (pt_on_flat_2 - pt_on_flat_1).cross(pt_on_flat_3 - pt_on_flat_1);
  Eigen::Vector3d n_x;
  n_x << 1, 0, 0;
  Eigen::Vector3d rot_to_x_axis;
  rot_to_x_axis = n_approx_x.cross(n_x);
  rot_to_x_axis.normalize();

  double theta_rot_x = std::acos(n_x.dot(n_approx_x)/(n_approx_x.norm()));

  Eigen::Matrix3d rot_mat_x;
  rot_mat_x = Eigen::AngleAxisd(theta_rot_x, rot_to_x_axis); //does the rot_to_z_axis need to be normalized?

  RV_2 = rotate_mesh(RV_1, rot_mat_x);

  const int grid_size = 20;
  const int time_steps = 100;
  const double isolevel = 0; //defalt was 0.1
  std::cerr<<"Computing swept volume...";
  igl::copyleft::swept_volume(RV_2,F,transform_rot_x,time_steps,grid_size,isolevel,SV,SF);
  std::cerr<<" finished."<< std::endl;

  std::cout << "reading patient" << std::endl;
  read_triangle_mesh("/home/brend/BCCancer/meshes/1e5pts_patient_and_carbon.ply",Vp,Fp);
  //has trouble reading the  w 1e5 pts?

  std::cout << "done reading patient" << std::endl;

  Eigen::MatrixXd TVp, RV_1p, RV_2p, swept_trans_linac;
  TVp = translate_mesh(Vp, -center_x, -center_y, -center_z);

  std::cout << "rotating patient 2x" << std::endl;
  RV_1p = rotate_mesh(TVp, rot_mat_z);
  RV_2p = rotate_mesh(RV_1p, rot_mat_x);

  double move_linac = -300; // [mm]
  swept_trans_linac = translate_mesh(SV, move_linac, 0, 0);
  
  std::cout << "checking intersection...";
  Eigen::MatrixXi IF;
  const bool b = false;
  bool a = igl::copyleft::cgal::intersect_other(swept_trans_linac,SF,RV_2p,Fp, b, IF); //this will give a 'true' for intersection
  // bool a = igl::copyleft::cgal::intersect_other(SV,SF,RV_2p,Fp, b, IF); //this will give a 'false' for intersection
  std::cout << "done checking intersection" << std::endl;
  std::cout << a << std::endl;

  // are the IF faces the intersecting faces? is that corresponding to the first or second mesh? would assume first 

  std::cout << "here is IF:" << std::endl;
  std::cout << IF << std::endl;

  Eigen::MatrixXd VI;
  Eigen::MatrixXi FI;
  std::cout << "Getting intersecting faces" << std::endl;
  IF_to_VandF(IF, VI, FI, swept_trans_linac, SF, true);
  std::cout << "Finished intersection work" << std::endl;

  igl::opengl::glfw::Viewer viewer2;

  viewer2.data().set_mesh(RV_2p,Fp);
  viewer2.data().show_lines = false;
  viewer2.data().set_face_based(true);
  add_viewer_axes(viewer2);

  viewer2.append_mesh();

  viewer2.data().set_mesh(swept_trans_linac,SF);
  viewer2.data().show_lines = false;
  viewer2.data().set_face_based(true);
  add_viewer_axes(viewer2);

  //set colors randombly for the viewer
  for (auto &data : viewer2.data_list){
      data.set_colors(0.5*Eigen::RowVector3d::Random().array() + 0.5);
  } 

  viewer2.launch();

}