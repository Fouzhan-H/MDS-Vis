#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <map>

#include <cstring>

#include "MakeGrid.cpp"
#include "parseGroFile.cpp"

// FilteredTypes => std::vector<Comp2Atom>
std::unique_ptr<std::vector<Comp2Atom>> listComp2Atoms (std::vector<FilteredTypes> const & ftypes, std::vector<Complex> const & clxs){
    
  std::unique_ptr<std::vector<Comp2Atom>> c2aFiltrd (new std::vector<Comp2Atom>());
  
  for (auto & c: clxs){
    for (auto &  t: ftypes){ 
      if (c.rName.compare(t.name) == 0){ // Filter complexes which has to be included in the analysis, and copy them to fClxs              
         Comp2Atom c2a; 
         // TODO c2a.atomIdxs.reserve(t.atoms.size()); 
         for (auto & a: t.atoms) c2a.atomIdxs.push_back(c.atom_idx+a);	 
	 c2a.nrAtoms = t.atoms.size();
         (*c2aFiltrd).push_back(std::move(c2a));
         break; 
      } 
    } // loop on filtered types (ftypes)
  } // loop on complexes (clxs)   

  return c2aFiltrd;
} 

void process_xtc(char * xtcFName, GroData const & groData, const std::array<float, 6> & range, std::unique_ptr<rvec[]> ps, std::vector<FilteredTypes> const & ftypes, char * outBsFlName){
  char outFlName [80];
  int status; 
  int atomNr;  
  
  matrix box; 
  float prec;
  float time;  
  int step; 

  // TODO Allow Filtering based on complex or atom 
  int fltClNr = groData.cs.size(); 
  
  std::unique_ptr<std::vector<Comp2Atom>>  c2aFiltrd;

  c2aFiltrd = listComp2Atoms(ftypes, groData.cs);

/*  std::unique_ptr<std::vector<Comp2Atom>>  c2aFiltrd (new std::vector<Comp2Atom>(fltClNr));
  std::transform( std::begin(groData.cs), std::end(groData.cs), std::begin(*c2aFiltrd)
		, [](const Complex & c){
		     Comp2Atom r;
		     r. = c.atom_idx;
		     r.nrAtoms = c.atom_no; 
		     return r;
		  }); */
  
  // open the trajectory file 
  XDRFILE * xtcFPtr = xdrfile_open(xtcFName, "r");
  if (xtcFPtr == NULL) std::cout << "Failed to Open the input file: " <<  xtcFName << std::endl; 


  fCoord cellSize {.5, .5, .5};  
  iCoord cellNo;
  cellNo[0] = (int) ((range[1] - range[0]) / .5) + 1;  
  cellNo[1] = (int) ((range[3] - range[2]) / .5) + 1;  
  cellNo[2] = (int) ((range[5] - range[4]) / .5) + 1;  
  MakeGrid mkGd(cellSize, cellNo, {range[0], range[2], range[4]}, *c2aFiltrd, reinterpret_cast<fCoord const *>(ps.get())); //TODO: reinterpret_cast can be problematic 
  

  // iteratively read the file and generate the grid 
  status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  for (int i = 0; i < 1000 && status == exdrOK; i++){
     std::cout << "----- iteration " <<  i << " " << time << " " << box[0][0] << " " << box[1][1] << " " << box[2][2] << " " << box[0][1] << " " << box[0][2] << "-----------" << std::endl;
     // TODO process this frame
     mkGd.getNextFrame(reinterpret_cast<fCoord const *>(ps.get()));
// mkGd.print();
     sprintf(outFlName, "%s%d", outBsFlName,  i);
     mkGd.writeFrame(outFlName);
     // TODO write to output file sprintf(ofname, "%s%d", ofname_base,  i); writeCVS(ofname, atomNr, as);
     status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  }

}


int main(int argc, char ** argv){
 
  // TODO input and outpu files should be provided as arguments
  //char outBsFlName [] = "/usr/not-backed-up/Datasets/MDS/VTK/test.vtk."; 
  char outBsFlName [100] { "/Users/fouzhanhosseini/Documents/BioVis/Data/VTK/test.vtk."}; 
  //char xtcFlName [] = "/usr/not-backed-up/Datasets/MDS/md.Fouzhan.xtc";
  char xtcFlName [100] = "/Users/fouzhanhosseini/Documents/BioVis/Data/md.Fouzhan.xtc";
  //char groFlName []= "/usr/not-backed-up/Datasets/MDS/starting_frame.Fouzhan.gro";
  char groFlName [100]= "/Users/fouzhanhosseini/Documents/BioVis/Data/starting_frame.Fouzhan.gro";

  if (argc > 1){
    std::strcpy (groFlName, argv[1]); 
    std::strcpy (xtcFlName, argv[2]); 
    std::strcpy (outBsFlName, argv[3]); 
  }

std::cout << "start " << xtcFlName << " " << groFlName << "\n"; 

  // Read number of atoms from the file and initialize the necessary data structures 
  int status; 
  int atomNr;  
  status = read_xtc_natoms(xtcFlName, &atomNr);
  GroData groData(atomNr); 
  std::unique_ptr<rvec []> ps (new rvec[atomNr]);
 
  readGroFile(groFlName, groData, ps.get());

  auto range = analyseAtomsPos(ps.get(), atomNr);

  auto filteredTypes = filterComplexs(groData);
  
//  process_xtc(xtcFlName, groData, range, std::move(ps), *filteredTypes, outBsFlName);   
}
