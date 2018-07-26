#include <vector>
#include <array>
#include <cmath>
#include <sstream>
#include <iostream>
#include <memory>

// Represents density of one type of lipids or Proteins 
class RotBox{
private:
  std::array< std::array<int,2>, 3> orgBox; // original box size 
  std::array <std::array<float,2>, 3> rotBox;  // the bigger box which always encompass the rotating original box 
  std::array <float, 3> cellSize; 
  std::array <unsigned int, 3> cellNu;
  // Given the original box dimension, compute the rotating box dimension
  void compRotBox(){
    // Assuming 1) box center is center of rotation, 
    // and 2) the rotation is mainly in xy plane
    // if the ccenter of rotation is to be determined dynamically,
    // then CoRx, CoRy, px, and py should be initialized accordingly. 
    float CoRx = orgBox[0][1]/2.0;
    float CoRy = orgBox[1][1]/2.0;
    float px = std::sqrt(CoRx*CoRx + CoRy*CoRy); // longest distanse between the center of rotation and a boundary point (x-dimension)    
    float py = std::sqrt(CoRx*CoRx + CoRy*CoRy); // longest distanse between the center of rotation and a boundary point (y-dimension)       
    // x-dimension 
    rotBox[0][0] = std::floor(CoRx - px); 
    rotBox[0][1] = std::ceil(CoRx + px);  
    // y-dimension 
    rotBox[1][0] = std::floor(CoRy - py); 
    rotBox[1][1] = std::ceil(CoRy + py);  
    // z-dimension 
    rotBox[2][0] = orgBox[2][0];  
    rotBox[2][1] = orgBox[2][1];  

    std::cout << "Original Box: [" << orgBox[0][1] << " " << orgBox[1][1] << " " << orgBox[2][1] << "]\n"; 
    std::cout << "Rotating Box: x:[" << rotBox[0][0] << " " << rotBox[0][1] << "]"
                           << " y:[" << rotBox[1][0] << " " << rotBox[1][1] << "]"
                           << " z:[" << rotBox[2][0] << " " << rotBox[2][1] << "]\n"; 
  } 
public:
  RotBox(const float (&box) [3][3],const float cells[3]){
    for (char i = 0; i < 3; i++){
      orgBox[i][0] = 0; 
      orgBox[i][1] = std::ceil(box[i][i]);
    }
    compRotBox();
    for (auto i:{0,1,2})
     cellSize[i] = cells[i];
    
    for (auto i:{0,1,2}){
      cellNu[i] = std::ceil((rotBox[i][1] - rotBox[i][0])/cellSize[i]);
    }
  }
  // Given the cell size, compute number of cells in each dimension
  void getCellNumbers(int cnu[3]) const{
    for (auto i:{0,1,2}){
      cnu[i] = cellNu[i];
    }
  }
  //return the lower point of the Rot Box
  void getLowestPoint(float org[3]){
    org [0] = rotBox[0][0];
    org [1] = rotBox[1][0];
    org [2] = rotBox[2][0]; 
  }
  // Given a coordinate, compute its mapping in the grid
  void cellIndex(const float pos [3], unsigned int cellIdx[3]) const{
    for (auto i:{0,1,2}){
      cellIdx[i] = std::floor((pos[i] - rotBox[i][0])/cellSize[i]); 
      if (pos[i] < rotBox[i][0]) {
        cellIdx[i] = 0;  
        std::cout << "out of bound " << i << " : " << pos[i] << "\n"; 
      }
      if (pos[i] > rotBox[i][1]) {
        cellIdx[i] = cellNu[i] - 1; 
        std::cout << "out of bound " << i << " : " << pos[i] << "\n"; 
      }
    }
  } 

};


class DensityGrid {
private:
  std::vector <std::vector<std::vector<int>>> densMap; // 3D array which records number of lipids in each cell over time
  std::unique_ptr <std::vector <unsigned int>> atoms; // List of atom ids   
  int dim [3];
public:
  // Consturctor
  DensityGrid(int const size[3], std::vector<unsigned int>* as)
    : densMap(size[2], std::vector<std::vector<int>>(size[1], std::vector<int>(size[0]))), atoms(as){
    for (auto i:{0,1,2}) dim[i] = size[i]; 
    std::cout << "Lipid density map created with " << atoms->size() << " atoms\n";
  } 
  // Given a list of atom positions and RotBox 
  void updDensity(const float (*ps) [3], const RotBox & box){
    unsigned int idx[3] ={0,0,0};
    
    for (auto a:*atoms){
      //map ps[a] to grid cell and update the value
      box.cellIndex(ps[a], idx);
      densMap[idx[2]][idx[1]][idx[0]] += 1;  
    }
  }
  
  friend std::ostream & operator<<(std::ostream & os, DensityGrid const & grid){
    std::ostringstream buff (std::ios_base::out); 
    for (int i = 0; i < grid.dim [2]; i++){
      for (int j = 0; j < grid.dim[1]; j++){
        for (int k = 0; k < grid.dim[0]; k++)
          buff << grid.densMap[i][j][k] << " ";              
        buff << std::endl; 
      } 
    }    
    os << buff.str();
    return os;
  }  

};


class ProtDensityGrid {
private:
  std::vector <std::vector<std::vector<int>>> densMap; // 3D array which records number of lipids in each cell over time
  unsigned int fstAtomIdx;    
  unsigned int lstAtomIdx;    
  int dim [3];
public:
  // Consturctor
  ProtDensityGrid(int const size[3], unsigned int fIdx, unsigned int lIdx)
    : densMap(size[2], std::vector<std::vector<int>>(size[1], std::vector<int>(size[0]))),fstAtomIdx(fIdx), lstAtomIdx(lIdx)  
  {
     std::cout << "Protein density map created with " << (lstAtomIdx - fstAtomIdx + 1) << " atoms\n";
     for (auto i:{0,1,2}) dim[i] = size[i];
  } 
  // Given a list of atom positions and RotBox 
  void updDensity(const float (*ps) [3], const RotBox & box){
    unsigned int idx[3];
    for (unsigned int i = fstAtomIdx; i <= lstAtomIdx; i++){
      //map ps[a] to grid cell and update the value 
      box.cellIndex(ps[i], idx);
      densMap[idx[2]][idx[1]][idx[0]] += 1;  
    }
  }

  friend std::ostream & operator<<(std::ostream & os, ProtDensityGrid const & grid){
    std::ostringstream buff (std::ios_base::out); 
    for (int i = 0; i < grid.dim [2]; i++){
      for (int j = 0; j < grid.dim[1]; j++){
        for (int k = 0; k < grid.dim[0]; k++)
          buff << grid.densMap[i][j][k] << " ";              
        buff << std::endl; 
      } 
    }    
    os << buff.str();
    return os;
  }  

};
