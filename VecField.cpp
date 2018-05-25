// The necessary inforamtion to build a VECTOR FIELD:
//   * Grid Size (dimensions and steps)
//   * the selected lipids (in each grid cell)
//   * Center of mass for each selected lipid at time t_i
//   * Center of mass for each selected lipid at time t_i+1
//   *
//

#include <utility>
#include <array>
#include <vector>
#include <memory>
#include <map>
#include <numeric>
#include <functional>
#include <iostream>
#include <fstream>

#define DIM 3
typedef std::array <float, DIM> fCoord;
typedef std::array <int, DIM> iCoord;

fCoord operator-(fCoord const & a, fCoord const & b){
  fCoord c;
  for (int i = 0; i < DIM; i++){
    c[i] = a[i] - b[i]; 
  }
  return c; 
}

fCoord operator+(fCoord const & a, fCoord const & b){
  fCoord c;
  for (int i = 0; i < DIM; i++){
    c[i] = a[i] + b[i]; 
  }
  return c; 
}


fCoord operator*(fCoord a, float s){
  fCoord c;
  for (int i = 0; i < DIM; i++){
    c[i] = s*a[i]; 
  }
  return c; 
}

class VecField{
    private:
      struct cellVec{
        fCoord v1;
        fCoord v2;
        unsigned int Nu; 
      };  

      fCoord cell_sz;
      fCoord grid_st; 
      iCoord cell_nu; 
      std::map<unsigned int, cellVec> a2cMapping; 
      std::vector<fCoord> cellVF;  
      unsigned int getCellIdx(fCoord p);  

      void mapAtom2Cell(fCoord const & x_ti, fCoord const & x_tii){
        unsigned int cidx;
        // find the atom cell index
        cidx = getCellIdx(x_ti); 
        // update the cell velocity at ti and tii  
        auto & cvec = a2cMapping[cidx];   
        for (int i = 0; i < DIM; i++){
          cvec.v1[i] += x_ti[i];
          cvec.v2[i] += x_tii[i]; 
        }
        cvec.Nu++;
      }

    public:
      VecField(fCoord cellSize, iCoord cellNu, fCoord st); 
      void getNextFrame(fCoord const * ps_ti, fCoord const * ps_tii,
                        std::vector<std::array<unsigned int, 2>> const & proteins,
                        std::vector<unsigned int> const & atoms); 
      void writeFrame(const char * fn);
//      void print();
};

VecField::VecField(fCoord cellSize, iCoord cellNu, fCoord st) 
  : cell_sz(cellSize), cell_nu(cellNu), grid_st(st){
  a2cMapping.clear();
  std::cout  << "DIMENSIONS " << cell_nu[0] +1<< " " << cell_nu[1] +1<< " " << cell_nu[2] +1 << std::endl;
}


unsigned int VecField::getCellIdx(fCoord p){

  iCoord cidx;
  for (int i = 0; i < DIM; i++){
    cidx[i] = (int) ((p[i] - grid_st[i]) / cell_sz[i]) % cell_nu[i];
//    cidx[i] = (cidx[i] >= cell_nu[i]) ? cell_nu[i] - 1: cidx[i];
  }

  unsigned int cellIdx = cidx[0] + cidx[1]* cell_nu[0] + cidx[2]* cell_nu[0]*cell_nu[1];  //TODO it is a hack 
  return cellIdx; 

}  


void VecField::getNextFrame(fCoord const * ps_ti, fCoord const * ps_tii, 
                            std::vector<std::array<unsigned int, 2>> const & proteins, 
                            std::vector<unsigned int> const & atoms){
 
  a2cMapping.clear();
  cellVF.clear();

  int k = 1 ; //13585;
  std::cout << ps_ti[k][0] << " "  << ps_tii[k][0] << " " <<  ps_ti[k][1] << " " 
            << ps_tii[k][1] << " " << ps_ti[k][2] << " " <<  ps_tii[k][2] << " " << getCellIdx(ps_ti[k]);


  for (auto & p:proteins){ // loop on each protein
    for (int i = p[0]; i <= p[1]; i++){ // loop on each amino-acid atom
      mapAtom2Cell(ps_ti[i], ps_tii[i]);
    }
  }


  for (auto & a:atoms){
    mapAtom2Cell(ps_ti[a], ps_tii[a]);
  }


  cellVec cvec;
  float tmp;
  for (auto & e:a2cMapping){
    cvec = e.second;
    tmp = 1.0 / cvec.Nu;

    cellVF.push_back( (cvec.v2 - cvec.v1) * tmp );  // calculate mv, (v2 - v1)/Nu, per cell
  }
   
} 


// TODO void VecField::print(){
//     std::vector <int> indcies {0, 1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18}; 	
//     for (int i : indcies)
//       std::cout << "ti_CoM[" << i<< "]: " << (*ti_CoM)[i][0] << " " << (*ti_CoM)[i][1] << " " << (*ti_CoM)[i][2] << std::endl; 
//}


void VecField::writeFrame(const char* fn){
   std::ofstream vtkFl; 
   vtkFl.open(fn, std::ofstream::out);
   
   vtkFl << "# vtk DataFile Version 2.0" << std::endl;
   vtkFl << "Center of Mass Trajectories" << std::endl;
   vtkFl << "ASCII" << std::endl; 
   // Geometry/Topology of the dataset
   vtkFl << "DATASET STRUCTURED_POINTS" << std::endl;
   vtkFl << "DIMENSIONS " << cell_nu[0] +1<< " " << cell_nu[1] +1<< " " << cell_nu[2] +1 << std::endl;
   vtkFl << "ORIGIN "<< "0 0 0" << std::endl;
   vtkFl << "SPACING "<< cell_sz[0] << " " << cell_sz[1] << " " << cell_sz[2] << std::endl;
   // Dataset attributes
   unsigned int nrCellData = cell_nu[0]*cell_nu[1]*cell_nu[2]; 
   vtkFl << "CELL_DATA " << nrCellData << std::endl;
   vtkFl << "VECTORS "<< "dp"<< " float" << std::endl;

   unsigned int cellIdx = 0; 
   unsigned int eIdx = 0; 

   for (auto e: a2cMapping){
     for (; cellIdx < e.first; ++cellIdx){
       vtkFl << "0 0 0" << std::endl;  
     }
     vtkFl << cellVF[eIdx][0] << " " << cellVF[eIdx][1] << " " << cellVF[eIdx][2] << std::endl;
//     std::cout << cellVF[eIdx][0] << " " << cellVF[eIdx][1] << " " << cellVF[eIdx][2] << std::endl;
     cellIdx++; 
     eIdx++;
   }

//   std::cout << "After loop " << nrCellData << " " << cellIdx <<"\n";
    
   //Wirte (0,0,0) for the remaining cells
   for (unsigned int i = cellIdx; i < nrCellData; i++) 
       vtkFl << "0 0 0" << std::endl;  
   

   vtkFl.close();
}

