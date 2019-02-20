#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>

#include <cstring>
#include "MakeGrid.cpp"
#include "AtomTraj.cpp"
#include "VecField.cpp"
#include "ParseGroFile.cpp"
#include "DensityGrid.cpp"

std::unique_ptr<std::vector<unsigned int>> listAtoms (FilteredTypes const & ftype, std::vector<Complex> const & clxs){

  std::unique_ptr<std::vector<unsigned int>> atoms (new std::vector<unsigned int>());
  
  for (auto & c: clxs){
      if (c.rName.compare(ftype.name) == 0){ // Filter complexes which has to be included in the analysis, and copy them to fClxs              
         for (auto & a: ftype.atoms){ 
           atoms->push_back(a+c.atom_idx);
         } 
      } 
  } // loop on complexes (clxs)   
  return atoms; 
}

void writeVTKFile(const char * fn, const int dim [3], const float origin [3], const float spacing[3], 
                  std::vector<ProtDensityGrid> const & prots , std::vector<DensityGrid> const & lipids, std::vector<FilteredTypes> const & typs, std::vector<std::string> aAcidNames){
   
   std::ofstream vtkFl; 
   vtkFl.open(fn, std::ofstream::out);
   
   vtkFl << "# vtk DataFile Version 2.0" << std::endl;
   vtkFl << "Density Map" << std::endl;
   vtkFl << "ASCII" << std::endl; 
   // Geometry/Topology of the dataset
   vtkFl << "DATASET STRUCTURED_POINTS" << std::endl;
   vtkFl << "DIMENSIONS " << dim[0] +1<< " " << dim[1] +1<< " " << dim[2] +1 << std::endl;
   vtkFl << "ORIGIN "<<  origin[0] << " " << origin[1] << " " << origin[2]  << std::endl; 
   vtkFl << "SPACING "<< spacing[0] << " " << spacing[1] << " " << spacing[2] << std::endl;
   // Dataset attributes
   unsigned int nrCellData = dim[0]*dim[1]*dim[2]; 
  

   vtkFl << "CELL_DATA " << nrCellData << std::endl;
   int countr = 1; 
   for (auto & p:prots){
     for (unsigned char i = 0; i < aAcidNames.size(); i++){
       vtkFl << "SCALARS prot"<< countr << "_" <<  aAcidNames[i] << " int 1" << std::endl;
       vtkFl << "LOOKUP_TABLE default"<< std::endl;
       p.writeAminoAcid(vtkFl, i); 
     }
     countr++; 
   } 

   countr = 0; 
   for (auto & l:lipids){
     vtkFl << "SCALARS "<< typs[countr].name << " float 1" << std::endl;
     vtkFl << "LOOKUP_TABLE default"<< std::endl;
     countr++; 
     vtkFl << l;
   } 

   vtkFl.close();
}

void calDensity(char const * xtcFName, char const * vtkBsFName, int atomNr, unsigned int deltaFrame,  std::vector<FilteredTypes> const & ftypes, GroData const & groData ){
  matrix box; 
  float prec;
  float time;  
  int step; 
  int status; 
//  int atomNr;  
  // Open the trajectory file 
  XDRFILE * xtcFPtr = xdrfile_open(xtcFName, "r");
  if (xtcFPtr == NULL) std::cout << "Failed to Open the input file: " <<  xtcFName << std::endl; 
  
  // Read the first frame
  std::unique_ptr<rvec []>  ps (new rvec[groData.atoms.size()]); 
  status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  if (status != exdrOK) std::cout << "Failed to read the xtc input file: "  << xtcFName << std::endl; 

  // Create RotBox
  float cells[] = {.5, .5,.5}; // TODO ask user for values
  RotBox rbox(box, cells);  


  int cellNu [3]; 
  rbox.getCellNumbers(cellNu);
std::cout << "cellNU " << cellNu[0] << " " << cellNu[1] << " " << cellNu[2]  << " " << atomNr <<" ..\n";
  // For each t in ftypes 1) extract atom indcies, and 2) create a density object
  std::vector<DensityGrid> lipidDenMaps; 
  std::unique_ptr<std::vector<unsigned int>> atomIdxs; 
  for (auto & t:ftypes){
    atomIdxs = listAtoms(t, groData.lipids);
    lipidDenMaps.push_back(DensityGrid(cellNu, atomIdxs.release()));   
  }
  
  // TODO how to output the file -- prepare data needed for output VTK files
  char vtkFName [100]; 
  float origin [3];
  rbox.getLowestPoint(origin);

  std::vector<std::string> aaNames;
  for (auto & a:groData.aatypes){
    aaNames.push_back(a.name);
  }

  // For every protein create a density object
  std::vector<ProtDensityGrid> protDenMaps; 
  for (auto & p:groData.ps){
    protDenMaps.push_back(ProtDensityGrid(cellNu, p.fstAtomIdx, p.lstAtomIdx, p.aminoAcids));   
  } 
  
std::cout << "HERE.  all created "<< groData.aatypes.size()<< " " << groData.ctypes.size() << "\n"; 

  // Calculate density of Proteins using the based frame 
  status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  for (auto & p:protDenMaps)
    p.updDensity(ps.get(), rbox);

  // Iterate over simulation and calculate Lipid denisties in deltaFrames 
  for (int i, dt = 0; status == exdrOK; dt++){
    for (i = 1; i < deltaFrame && status == exdrOK; i++){
//     std::cout << "----- iteration " <<  i << " "  
//                                     << box[0][0] << " " << box[1][1] << " " << box[2][2] << " "
//                                     << time << " " << step << "-----------" << std::endl;
      // process this frame
      for (auto & l:lipidDenMaps)
        l.updDensity(ps.get(), rbox); 

      status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);    
    }

    // Normalize the Lipids density figures
    for (auto & l:lipidDenMaps)
//      l.normalize( dt * deltaFrame + i );
      l.normalize( i );

    sprintf(vtkFName,"%s%d.vtk", vtkBsFName, dt );
    writeVTKFile (vtkFName, cellNu, origin, cells, protDenMaps, lipidDenMaps, ftypes, aaNames); 

    for (auto & l:lipidDenMaps) l.reset();

  }
  xdrfile_close(xtcFPtr); // Close the input xtc file


//std::cout << frameNu << " HERE. before returng callDen \n"; 
}

int main(int argc, char ** argv){
 
  // TODO input and outpu files should be provided as arguments -- add option flags -- test on other platform?
  //char outBsFlName [] = "/usr/not-backed-up/Datasets/MDS/VTK/test.vtk."; 
  char outFlName [100] { "/Users/fouzhanhosseini/Documents/BioVis/Data/VTK/test.vtk"}; 
  //char xtcFlName [] = "/usr/not-backed-up/Datasets/MDS/md.Fouzhan.xtc";
  char xtcFlName [100] = "/Users/fouzhanhosseini/Documents/BioVis/Data/md.Fouzhan.xtc";
  //char groFlName []= "/usr/not-backed-up/Datasets/MDS/starting_frame.Fouzhan.gro";
  char groFlName [100]= "/Users/fouzhanhosseini/Documents/BioVis/Data/starting_frame.Fouzhan.gro";

  if (argc > 1){
    std::strcpy (outFlName, argv[1]);
    std::strcpy (groFlName, argv[2]); 
    std::strcpy (xtcFlName, argv[3]);
  }

  // Inquire about number of proteins included in the file and their first & last atom number
  int np, fa, la;
  std::cout << "Please Enter Number of Proteins" << std::endl; 
  std::cin >> np; 
  std::vector<std::pair<int, int>> pns; 
  for (int i = 0; i < np; i++){
    std::cout << "Please Enter first and last atom number" << std::endl;
    std::cin >> fa >> la ;    
    pns.push_back(std::pair<int, int>(fa, la));
  }
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(),'\n');	
 

  // Read number of atoms from the file and initialize the necessary data structures 
  int status; 
  int atomNr;  
  status = read_xtc_natoms(xtcFlName, &atomNr);
  GroData groData(atomNr); 
  std::unique_ptr<rvec []> ps (new rvec[atomNr]);
  ReadGro rGro;
  rGro.readGroFile(groFlName, groData, ps.get(), pns);
  //testGroReader(groData, ps.get());
  
  auto filteredTypes = filterComplexs(groData);
  
  calDensity(xtcFlName, outFlName, atomNr, std::stoi(argv[4]), *filteredTypes, groData);
}
