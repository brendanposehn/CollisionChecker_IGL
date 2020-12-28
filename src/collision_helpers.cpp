#include "collision_helpers.h"

Eigen::MatrixXd translate_mesh(Eigen::MatrixXd V, double x, double y, double z){
  //V is Zx3 where Z is numpber of vertices

  Eigen::Vector3d translation;
  translation << x, y, z;

  Eigen::MatrixXd VT;
  VT = V;

// iterate over every vertex and apply translation
  for (int i = 0; i < V.rows(); i++){  
    Eigen::Vector3d currRow;
    currRow = VT.row(i);
    // VT.row(i) = VT.row(i) + translation; //sizes are the same?
    VT.row(i) = currRow + translation; //this works!
  }

  return VT;

} 

void add_viewer_axes(igl::opengl::glfw::Viewer &viewer){
// void add_viewer_axes(int a){ //works fine when the input is just an integer so stuff gets messed up when we are using the viewer

  Eigen::Vector3d origin;
  origin << 0, 0, 0;
  Eigen::Vector3d unit_x;
  unit_x << 1000, 0, 0;
  Eigen::Vector3d unit_y;
  unit_y << 0, 1000, 0;
  Eigen::Vector3d unit_z;
  unit_z << 0, 0, 1000;

  Eigen::Vector3d red;
  red << 1, 0, 0;
  Eigen::Vector3d green;
  green << 0, 1, 0;
  Eigen::Vector3d blue;
  blue << 0, 0, 1;

  Eigen::Matrix3d starts;
  starts << origin, origin, origin;

  Eigen::Matrix3d ends;
  ends << unit_x, unit_y, unit_z;

  Eigen::Matrix3d colors;
  colors << red, green, blue;

  viewer.data().add_edges(starts, ends, colors);

}

Eigen::MatrixXd rotate_mesh(Eigen::MatrixXd V, double theta, char dir){

  Eigen::Matrix3d rot_mat;

  if (dir == 'x'){
    rot_mat <<  1,         0,                0,
                0, std::cos(theta), -std::sin(theta),
                0, std::sin(theta),  std::cos(theta);
    std::cout << "rot x" << std:: endl;
  }else if (dir == 'y'){
    rot_mat << std::cos(theta), 0, std::sin(theta),
               0,               1,               0,
              -std::sin(theta), 0, std::cos(theta);
    std::cout << "rot y" << std::endl;
  }else{
    rot_mat << std::cos(theta), -std::sin(theta), 0,
               std::sin(theta),  std::cos(theta), 0,
               0,                0,               1;
    std::cout << "rot_z" << std::endl;
  }

  Eigen::MatrixXd VR;
  VR = V;

// iterate over every vertex and apply rotation
  for (int i = 0; i < V.rows(); i++){  
    Eigen::Vector3d currRow;
    currRow = VR.row(i).transpose();
    VR.row(i) = (rot_mat * currRow).transpose(); 
  }

  return VR;

} 

Eigen::MatrixXd rotate_mesh(Eigen::MatrixXd V, Eigen::Matrix3d rot_mat){

  Eigen::MatrixXd VR;
  VR = V;

// iterate over every vertex and apply rotation
  for (int i = 0; i < V.rows(); i++){  
    Eigen::Vector3d currRow;
    currRow = VR.row(i).transpose();
    VR.row(i) = (rot_mat * currRow).transpose(); 
  }

  return VR;

}

void add_meshes(Eigen::MatrixXd VA, Eigen::MatrixXi FA, Eigen::MatrixXd VB, Eigen::MatrixXi FB, Eigen::MatrixXd &V, Eigen::MatrixXi &F){
  V.resize(VA.rows()+VB.rows(),VA.cols());
  V<<VA,VB;
  F.resize(FA.rows()+FB.rows(),FA.cols());
  F<<FA,(FB.array()+VA.rows());
}

void IF_to_VandF(Eigen::MatrixXi IF, Eigen::MatrixXd &VI, Eigen::MatrixXi &FI, Eigen::MatrixXd V, Eigen::MatrixXi F, bool first){
  //IF contains the faces that intersect from both
  //V and F are the mesh we want to plot intersections of 
  int c = 1;
  //if we want the intersecting faces from the first entry
  if (first){
    c = 0;
  }

  //if want to change to list:
  //https://stackoverflow.com/questions/24139428/check-if-element-is-in-the-list-contains

  std::cout << "IF rows: " << IF.rows() << ". IF cols: " << IF.cols() << std::endl;


  //faces are the row indices of F that intersected
  Eigen::VectorXi int_faces = IF.col(c);
  
  std::set<int> int_faces_set;
  for(int i=0; i<int_faces.size(); i++){
    int_faces_set.insert(int_faces(i));
  }

  std::set<int> vert_inds;

  //METHOD 1: Using the list (lots of extras)

  //need storage for which vertices are together especially consdiering could already have been added in
  // Eigen::MatrixXi F_full(int_faces.size(), 3); //vertex indices wrt full matrix
  // for(int i=0; i<int_faces.size(); i++){
  //   Eigen::Vector3i r = F.row(int_faces(i));
  //   F_full.row(i) = r; //the old vertex indices
  //   vert_inds.insert(r(0));
  //   vert_inds.insert(r(1));
  //   vert_inds.insert(r(2));
  // }//the set is stored in increasing order, they dont keep their order


  //METHOD 2: Using the set 

  Eigen::MatrixXi F_full(int_faces_set.size(), 3); //vertex indices wrt full matrix
  int i = 0;
  for (auto it = int_faces_set.begin(); it != int_faces_set.end(); ++it){
    Eigen::Vector3i r = F.row(*it);
    F_full.row(i) = r;
    vert_inds.insert(r(0));
    vert_inds.insert(r(1));
    vert_inds.insert(r(2));
    i++;
  }
  //end METHOD 2

  for (auto it = vert_inds.begin(); it != vert_inds.end(); ++it){
    std::cout << *it << std::endl;
  }

  if (vert_inds.empty()){
    std::cout << "No intersections provided" << std::endl; 
    return;
  }

  std::cout << "Number vertices in intersecting faces: " << vert_inds.size() << std::endl;
  std::cout << "Total number of vertices in full mesh: " << V.rows() << std::endl;

  VI.resize(vert_inds.size(), 3);

  //make the new V matrix with all the vertices
  // std::set<int, std::greater<int> >::iterator itr; //what does the greater do???
  std::set<int>::iterator itr;
  int j = 0;
  for (itr = vert_inds.begin(); itr != vert_inds.end(); ++itr){
    // // *itr gives current element
    // std::cout << "VI access" << std::endl;
    // std::cout << VI.row(i) << std::endl;
    // std::cout << "V access" << std::endl;
    // std::cout << *itr << std::endl;
    // std::cout << V.row(*itr) << std::endl;
    VI.row(j) = V.row(*itr);
    i++;
  }

  //make new F matrix with proper references to VI
  //this is going to be very computationally complex LIKELY EXISTS BETTER WAY TO DO THIS
  //we can use the same ordering as in the set
  FI.resize(int_faces_set.size(), 3);
  //just want to go through vertex indices

  //these values are all waaaay too big
  for(int i=0; i<FI.rows(); i++){
    for(int j=0; j<=3; j++){
      //need to find indices of each vertex element
      int curr_vert = F_full(i, j);//this is a single OLD vertex index
      //the CURR_VERT values are way off, some are negative?
      // std::cout << "curr vert: " << curr_vert << std::endl;
      std::set<int>::iterator itr2;
      for (itr2 = vert_inds.begin(); itr2 != vert_inds.end(); itr2++){
        // std::cout << "itrs2: " << *itr2 << std::endl;
        if( *itr2 == curr_vert){
          int d = std::distance(vert_inds.begin(), itr2);
          //we never get inside of this loop
          // std::cout << << "d: " << d << std::endl;
          FI(i, j) = d;
          break;
        }
      } 
    }
  }  

  std::cout << "Here is FI:\n" << FI << std::endl;

}

