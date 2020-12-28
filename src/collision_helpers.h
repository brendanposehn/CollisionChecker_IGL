#ifndef COLLISION_HELPERS_H
#define COLLISION_HELPERS_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <set>

#include <igl/opengl/glfw/Viewer.h>


Eigen::MatrixXd translate_mesh(Eigen::MatrixXd V, double x, double y, double z);

void add_viewer_axes(igl::opengl::glfw::Viewer &viewer);

Eigen::MatrixXd rotate_mesh(Eigen::MatrixXd V, double theta, char dir);
Eigen::MatrixXd rotate_mesh(Eigen::MatrixXd V, Eigen::Matrix3d rot_mat);

void add_meshes(Eigen::MatrixXd VA, Eigen::MatrixXi FA, Eigen::MatrixXd VB, Eigen::MatrixXi FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F);

void IF_to_VandF(Eigen::MatrixXi IF, Eigen::MatrixXd &VI, Eigen::MatrixXi &FI, Eigen::MatrixXd V, Eigen::MatrixXi F, bool first);

#endif