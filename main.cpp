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

std::unique_ptr<std::vector<unsigned int>> listAtoms (std::vector<FilteredTypes> const & ftypes, std::vector<Complex> const & clxs, std::vector<std::string> & complexTypes){
  std::unique_ptr<std::vector<unsigned int>> atoms (new std::vector<unsigned int>());
  
  for (auto & c: clxs){
    for (auto &  t: ftypes){ 
      if (c.rName.compare(t.name) == 0){ // Filter complexes which has to be included in the analysis, and copy them to fClxs              
         for (auto & a: t.atoms){ 
           atoms->push_back(a+c.atom_idx);
           complexTypes.push_back(c.rName);   
         } 
         break; 
      } 
    } // loop on filtered types (ftypes)
  } // loop on complexes (clxs)   
  return atoms; 
}

void processXTC(char * xtcFName, GroData const & groData, const std::array<float, 6> & range, std::unique_ptr<rvec[]> ps_prev, std::vector<FilteredTypes> const & ftypes, char * outBsFlName){
  char outFlName [80];
  int status; 
  int atomNr;  
  
  matrix box; 
  float prec;
  float time;  
  int step; 
 
  // Extract atoms composing proteins
  std::vector<std::array<unsigned int, 2>> pAtoms; 
  for (auto & p: groData.ps){
    pAtoms.push_back({p.fstAtomIdx, p.lstAtomIdx});
  }

  // TODO Allow Filtering based on complex or atom 
  int fltClNr = groData.lipids.size(); 
  std::unique_ptr<std::vector<unsigned int>>  fAtoms;
  std::vector <std::string> types;
  fAtoms = listAtoms(ftypes, groData.lipids, types);
  std::cout << "atom size" <<  fAtoms->size() << "\n";


  // open the trajectory file 
  XDRFILE * xtcFPtr = xdrfile_open(xtcFName, "r");
  if (xtcFPtr == NULL) std::cout << "Failed to Open the input file: " <<  xtcFName << std::endl; 


  fCoord cellSize {1.0, 1.0, 1.0};  
  iCoord cellNo;
  cellNo[0] = (int) ((range[1] - range[0]) / 1) + 1;  
  cellNo[1] = (int) ((range[3] - range[2]) / 1) + 1;  
  cellNo[2] = (int) ((range[5] - range[4]) / 1) + 1;  
  VecField vf (cellSize, cellNo, {range[0], range[2], range[4]});
  
  std::unique_ptr<rvec []>  ps (new rvec[groData.atoms.size()]); 
  // iteratively read the file and generate the grid 
  status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  for (int i = 0; status == exdrOK; i++){
     std::cout << "----- iteration " <<  i << " " << time << " " << box[0][0] << " " << box[0][1] << " " << box[0][2] 
                                                          << " " << box[1][0] << " " << box[1][1] << " " << box[1][2]  
                                                          << " " << box[2][0] << " " << box[2][1] << " " << box[2][2] << "-----------" << std::endl;
     // process this frame
     vf.getNextFrame(reinterpret_cast<fCoord const *>(ps_prev.get())
                   , reinterpret_cast<fCoord const *>(ps.get())
                   , pAtoms, *fAtoms);

     sprintf(outFlName, "%s%d.vtk", outBsFlName,  i);
     vf.writeFrame(outFlName);
     ps.swap(ps_prev); 
     status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  }

}

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
  int fltClNr = groData.lipids.size(); 
  
  std::unique_ptr<std::vector<Comp2Atom>>  c2aFiltrd;

  c2aFiltrd = listComp2Atoms(ftypes, groData.lipids);

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
  for (int i = 0; i < 3 && status == exdrOK; i++){
     std::cout << "----- iteration " <<  i << " " << time << " " << box[0][0] << " " << box[1][1] << " " << box[2][2] << " " << box[0][1] << " " << box[0][2] << "-----------" << std::endl;
     // TODO process this frame
     mkGd.getNextFrame(reinterpret_cast<fCoord const *>(ps.get()));
// mkGd.print();
     sprintf(outFlName, "%s%d.vtk", outBsFlName,  i);
     mkGd.writeFrame(outFlName);
     // TODO write to output file sprintf(ofname, "%s%d", ofname_base,  i); writeCVS(ofname, atomNr, as);
     status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  }

}


void xtc2trajs(char * xtcFName, GroData const & groData, std::vector<FilteredTypes> const & ftypes, char * outFlName){
  //Allow Filtering based on complex or atom 
  std::unique_ptr<std::vector<unsigned int>>  fAtoms;
  std::vector <std::string> types;
  fAtoms = listAtoms(ftypes, groData.lipids, types);
  std::cout << "Number of atoms: " <<  fAtoms->size() << "\n";
  // extract atom trajectories 
  matrix box; 
  float prec;
  float time;  
  int step; 
  int status; 
  int atomNr;  
  // open the trajectory file 
  XDRFILE * xtcFPtr = xdrfile_open(xtcFName, "r");
  if (xtcFPtr == NULL) std::cout << "Failed to Open the input file: " <<  xtcFName << std::endl; 
  
  int frameNu = 6000; /*TODO it should either be provided as input or be removed */
  AtomTraj trajReader ( frameNu, fAtoms->size(), groData.box);  
  std::unique_ptr<rvec []>  ps (new rvec[groData.atoms.size()]); 
  // iteratively read the xtc file to extract trajectories  
  status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  for (int i = 0; i < frameNu && status == exdrOK; i++){
     std::cout << "----- iteration " <<  i << " " << time << " "
                                     << box[0][0] << " " << box[1][1] << " " << box[2][2] << " "
                                     << box[0][1] << " " << box[0][2] << " " << time << " " << step << "-----------" << std::endl;
     // process this frame
     trajReader.addFrame(ps.get(), fAtoms.get(), time, static_cast<const float (&)[3][3]> (box));
     status = read_xtc(xtcFPtr, atomNr, &step, &time, box, ps.get(), &prec);
  }

  trajReader.writeTrajs(outFlName, types); // dump out trajectories considering PBC 
}

int main(int argc, char ** argv){
 
  // TODO input and outpu files should be provided as arguments -- add option flags -- test on other platform?
  //char outBsFlName [] = "/usr/not-backed-up/Datasets/MDS/VTK/test.vtk."; 
  char outBsFlName [100] { "/Users/fouzhanhosseini/Documents/BioVis/Data/VTK/test.vtk."}; 
  //char xtcFlName [] = "/usr/not-backed-up/Datasets/MDS/md.Fouzhan.xtc";
  char xtcFlName [100] = "/Users/fouzhanhosseini/Documents/BioVis/Data/md.Fouzhan.xtc";
  //char groFlName []= "/usr/not-backed-up/Datasets/MDS/starting_frame.Fouzhan.gro";
  char groFlName [100]= "/Users/fouzhanhosseini/Documents/BioVis/Data/starting_frame.Fouzhan.gro";

  if (argc > 1){
    std::strcpy (outBsFlName, argv[1]);
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

  // TODO assumption made in this function are not true, 
  //    first the range can change slightly frame by frame  
  //    box value gives the simulation range, i.e. includeing both protein and lipids; the other question is when do we need this funciton?
  // auto range = analyseAtomsPos(ps.get(), atomNr); 

  //testGroReader(groData, ps.get());
  auto filteredTypes = filterComplexs(groData);
  
//TODO  process_xtc(xtcFlName, groData, groData.box, std::move(ps), *filteredTypes, outBsFlName);   
//  processXTC(xtcFlName, groData, groData.box, std::move(ps), *filteredTypes, outBsFlName);   

  xtc2trajs(xtcFlName,groData, *filteredTypes, outBsFlName);
}
