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
#include <algorithm>
#include <numeric>
#include <functional>
#include <iostream>
#include <fstream>

#define DIM 3
typedef std::array <float, DIM> fCoord;
typedef std::array <int, DIM> iCoord;

struct Comp2Atom{
    std::vector<unsigned int> atomIdxs; 
    unsigned int nrAtoms;  
};

struct Complex2Grid{
    int complexId;
    unsigned int cellIdx;
    iCoord cellCoord; 
};

class MakeGrid{
    private:
      fCoord cell_sz;
      fCoord grid_st; 
      iCoord cell_no; 
      std::vector <Comp2Atom> const & c2aMapping; 
      std::unique_ptr<std::vector<fCoord>> ti_CoM; 
      std::unique_ptr<std::vector<fCoord>> tii_CoM; 
      std::vector<fCoord> dp; 
      std::vector<Complex2Grid> c2gMapping; 
      void calCoM (fCoord const * ps);
      void calMvs();      
      void mapClxs2Grid();      
    public:
      MakeGrid(fCoord cellSize, iCoord cellNo, fCoord st, std::vector<Comp2Atom> & c2aM, fCoord const * ps); 
      void getNextFrame(fCoord const * ps); 
      void writeFrame(const char * fn);
      void print();
};

MakeGrid::MakeGrid(fCoord cellSize, iCoord cellNo, fCoord st, std::vector<Comp2Atom> & c2aM, fCoord const * ps) 
  : cell_sz(cellSize), cell_no(cellNo), grid_st(st), dp(c2aM.size()), c2gMapping(c2aM.size()), c2aMapping(c2aM) {
      int complexNo = c2aMapping.size();	   
      std::cout << complexNo << " " << c2aM.size() << std::endl;
      ti_CoM = std::unique_ptr<std::vector<fCoord>>(new std::vector<fCoord>(complexNo)); 
      tii_CoM = std::unique_ptr<std::vector<fCoord>>(new std::vector<fCoord>(complexNo)); 
      calCoM(ps);
      ti_CoM.swap(tii_CoM);     
}

// Complex (start_idx, number of atoms) -> table of positions -> CoM (rvec)  
void MakeGrid::calCoM(fCoord const * ps){
 
/*TODO  std::vector <int> indcies {0, 1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18};
     for (int i : indcies)
       std::cout << "ps[" << i<< "]: " << ps[i][0] << " " << ps[i][1] << " " << ps[i][2] << std::endl; 
*/
    std::transform( std::begin(c2aMapping), std::end(c2aMapping), std::begin(*tii_CoM)
		  , [&ps](Comp2Atom const & e){   
                        int atom_no = e.nrAtoms;
 
			fCoord t{};                                            // Value-initialized to zero
                        t = std::accumulate( std::begin(e.atomIdxs), std::end(e.atomIdxs), t       // Accumulate the positions of all the atoms which form this complex
		                           , [&ps](const fCoord &a, unsigned int const i){ 
                                                 fCoord const & b = *(ps + i); 
		                                 fCoord r; 
		          	        	 std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(r), std::plus<float>());
			  		         return r; 
		 		           });
      		        std::transform( std::begin(t),std::end(t), std::begin(t)
				     , [&atom_no](const float &a){return a / atom_no;}); // Average the accumulated positions to calculate CoM

			return t; 
                  });
}   

void MakeGrid::calMvs(){
    // Compute the displacement (dp) of the complexes over time [t_i,t_i+1]
    std::transform( std::begin(*tii_CoM), std::end(*tii_CoM)  
                  , std::begin(*ti_CoM)
                  , std::begin(dp)
                  , [](const fCoord &p1, const fCoord &p2) -> fCoord {    // Calculate the displacement for one Complex
                        fCoord r;
                        std::transform(std::begin(p1), std::end(p1), std::begin(p2)
                                      , std::begin(r), [](const float &a, const float &b){ return a - b;});
                        return r;                         
                    });
}

void MakeGrid::mapClxs2Grid(){
    // For each complex find the grid cell (where the center of mass is located)
    int i = 0;
    std::transform( std::begin(*ti_CoM), std::end(*ti_CoM) 
                  , std::begin(c2gMapping)
                  , [&i, this](const fCoord c){
                        Complex2Grid r;
                        r.complexId = i; i++;
                        for (int i = 0; i < DIM; i++){
			  r.cellCoord[i] =  ((c[i] - grid_st[i]) / cell_sz[i]) ;
			  r.cellCoord[i] = (r.cellCoord[i] >= cell_no[i]) ? cell_no[i] - 1: r.cellCoord[i];
			}
                        r.cellIdx = r.cellCoord[0] + r.cellCoord[1]* cell_no[0] + r.cellCoord[2]* cell_no[0]*cell_no[1];  //TODO it is a hack 
			return r; 
                  });
    // sort them 
/*    std::sort(std::begin(c2gMapping), std::end(c2gMapping),
		    [](const Complex2Grid & a, const Complex2Grid & b){  // TODO why this comparison does not work?
		       for (int i = DIM -1 ; i > -1 ; i--){
		         if (a.cellCoord[i] > b.cellCoord[i]) 
		            return false; 
		       }
		       return true; 
		    }); */
      
     std::sort( std::begin(c2gMapping), std::end(c2gMapping)
	       , [](const Complex2Grid & a, const Complex2Grid & b){
	            return a.cellIdx < b.cellIdx;
	         });

}

void MakeGrid::getNextFrame(const fCoord * ps ){
    calCoM (ps);
    calMvs();      
    mapClxs2Grid();      
     
    ti_CoM.swap(tii_CoM);      
}


void MakeGrid::print(){
     std::vector <int> indcies {0, 1, 2, 3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18}; 	
     for (int i : indcies)
       std::cout << "ti_CoM[" << i<< "]: " << (*ti_CoM)[i][0] << " " << (*ti_CoM)[i][1] << " " << (*ti_CoM)[i][2] << std::endl; 
     for (int i : indcies)
       std::cout << "tii_CoM[" << i<< "]: " << (*tii_CoM)[i][0] << " " << (*tii_CoM)[i][1] << " " << (*tii_CoM)[i][2] << std::endl; 
     for (int i : indcies)
       std::cout << "dp[" << i<< "]: " << dp[i][0] << " " << dp[i][1] << " " << dp[i][2] << std::endl; 
     for (int i : indcies)
       std::cout << "c2gMapping[" << i<< "]: " << c2gMapping[i].cellIdx << " --- "<< c2gMapping[i].cellCoord[0] << " " << c2gMapping[i].cellCoord[1] << " " << c2gMapping[i].cellCoord[2] << std::endl; 
}


void MakeGrid::writeFrame(const char* fn){
   std::ofstream vtkFl; 
   vtkFl.open(fn, std::ofstream::out);
   
   vtkFl << "# vtk DataFile Version 2.0" << std::endl;
   vtkFl << "Center of Mass Trajectories" << std::endl;
   vtkFl << "ASCII" << std::endl; 
   // Geometry/Topology of the dataset
   vtkFl << "DATASET STRUCTURED_POINTS" << std::endl;
   vtkFl << "DIMENSIONS " << cell_no[0] +1<< " " << cell_no[1] +1<< " " << cell_no[2] +1 << std::endl;
   vtkFl << "ORIGIN "<< "0 0 0" << std::endl;
   vtkFl << "SPACING "<< cell_sz[0] << " " << cell_sz[1] << " " << cell_sz[2] << std::endl;
   // Dataset attributes
   unsigned int nrCellData = cell_no[0]*cell_no[1]*cell_no[2]; 
   vtkFl << "CELL_DATA " << nrCellData << std::endl;
   vtkFl << "VECTORS "<< "dp"<< " float" << std::endl;

   unsigned int cellIdx = 0; 
   unsigned int eIdx = 0; 
   for (Complex2Grid & e: c2gMapping){
       eIdx = e.cellIdx; 
       if (cellIdx > e.cellIdx){
	   continue;
       }
       else{ 
            if (cellIdx < e.cellIdx){
	        for (int i = cellIdx; i < e.cellIdx; i++) 
                     vtkFl << "0 0 0" << std::endl;  
	    }
       }
       vtkFl << dp[e.complexId][0] << " " << dp[e.complexId][1] << " " << dp[e.complexId][2] << std::endl;
       cellIdx=eIdx+1;
   }
//   std::cout << "After loop " << nrCellData << " " << cellIdx <<"\n";
    
   //Wirte (0,0,0) for the remaining cells
   for (unsigned int i = cellIdx; i < nrCellData; i++) 
       vtkFl << "0 0 0" << std::endl;  
   

   vtkFl.close();
}

