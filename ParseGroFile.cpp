#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <map>
#include <memory>

#include <cstring>
#include <cmath>

#define MAX_GRO_LINE_LENGTH 200 
#define GRO_NAME_SIZE       5  


// Represent one atom as according to Gro file 
struct Atom {
  char name [GRO_NAME_SIZE];  // atom name 
  int num;        // atom number
};

// Represent one Complex according to Gro format, 
// links with atoms composing the complex by keeping
// the index to the first atom in list of atoms and 
// number of atoms composing the complex
struct Complex { 
  Complex(int num, std::string const & name, int a_idx)
     : rNum(num), atom_idx(a_idx), rName(name) {
     }
  int rNum;           // residue number 
  std::string rName;  // residue name 
  int atom_idx;       // Index of the first atom in atoms list 
  int atom_no;        // Number of atoms composing the complex 
}; 


// This struct represetn an complex type, e.g. POPC
struct ComplexType{
  ComplexType(std::string n) : name(n){}; 
  std::string name;                  // Complex type name
  unsigned short nrAtoms;            // Number of atoms
  std::vector<std::string> atoms;    // list of atoms 

  friend std::ostream & operator<<(std::ostream & os, ComplexType const & ctype);
};

struct GroData {
  GroData(int n): atoms(n){}
  std::vector<ComplexType> ctypes; // list of Complex types 	
  std::vector <Atom> atoms;        // array of Atoms 
  std::vector <Complex> cs;        // array of Complexes 
  // std::map<std::string, std::vector<std::string>> a2ctype;  TODO do we need this?
};

struct FilteredTypes{
  int typeIdx; 
  std::string name;                  // Complex type name
  std::vector<int> atoms;
};

std::ostream & operator<<(std::ostream & os, ComplexType const & c){
  int cntr_a = 1; 
  os << c.name << " consists of : " << c.nrAtoms << "\n atoms:";
  for (auto & a:c.atoms){
    os << " " << cntr_a << "."<< a;
    cntr_a++;
  }
  os << "\n"; // TODO do I need a new line or flush?
  return os;
}


void readGroFile(const char * fn, GroData & groData, rvec * ps){
  char line [MAX_GRO_LINE_LENGTH];
  char * token;

  std::ifstream gf; 
  gf.open(fn, std::fstream::in);

  // Read the first lin,title string, ignore it for now 
  gf.ignore(50, '\n');     

  // Read number of atoms 
  gf.getline(line, MAX_GRO_LINE_LENGTH); 
  token = strtok(line, " ");
  int nr = std::atoi(token); 

  // Iterate over all lines in the file (one atom per line)
  // and read atoms information 
  std::string cTypeNamePrev;
  std::string cTypeName;
  int cNuPrev=-1, cNu = 0;
  bool newType = false; 
  char cNum [GRO_NAME_SIZE+1]; 
  cNum [GRO_NAME_SIZE]= 0x00;
  int a_cntr = 0;
  std::vector<std::string> asList; 
  for (int i = 0; i < nr ; i++){
    gf.getline(line,MAX_GRO_LINE_LENGTH);
    // Read the new complex and check whether it is a new complex 
    std::copy(line, line+GRO_NAME_SIZE, cNum); // copy residue number to tmp container 
    cNu = std::atoi (cNum);
    token = std::strtok(line+GRO_NAME_SIZE, " ");
    cTypeName.assign(token);

    if (cNuPrev == cNu && cTypeName.compare(cTypeNamePrev) == 0) 
       a_cntr++;                                              // not a new complex or complex type 
    else{ 
      if (newType){                                           // Set nrAtoms and atoms list for the last recorded complex type 
            groData.ctypes.back().atoms = std::move(asList);
	    groData.ctypes.back().nrAtoms = a_cntr;
	    newType = false;
       }	    
       if (cTypeName.compare(cTypeNamePrev) != 0){            // This is a new complex type 
	  groData.ctypes.push_back(ComplexType(cTypeName));
          newType = true;
       }

       if (groData.cs.size() > 0) groData.cs.back().atom_no = a_cntr;   // set atom_no of the last complex 
       groData.cs.push_back(Complex(cNu, cTypeName, i));                // store this new complex 
       // reset the values for the next round 
       asList.clear(); 
       a_cntr = 1;                             
       cNuPrev = cNu; 
       cTypeNamePrev.assign(cTypeName);  
    }

    // record the read atom information
    token = strtok(NULL, " ");
    asList.push_back(std::string(token));
    std::strcpy(groData.atoms[i].name, token);    // as[i].aName = strtok(NULL, " ");
    groData.atoms[i].num = std::atoi(strtok(NULL, " "));

    // record the atom position
    ps[i][0] =  std::atof(strtok(NULL, " "));
    ps[i][1] =  std::atof(strtok(NULL, " "));
    ps[i][2] =  std::atof(strtok(NULL, " "));
  }

  if (groData.cs.size() > 0) groData.cs.back().atom_no = a_cntr;  // set atom_no of the last complex 
  if (newType){                                                   // Set nrAtoms and atoms list for the last recorded complex type 
     groData.ctypes.back().atoms = std::move(asList);
     groData.ctypes.back().nrAtoms = a_cntr;
     newType = false;
  }	    
  // TODO read the box dimension if needed 

  gf.close();
}

void calHist(float min, float max, float * p, int n){
  float binSize = 0.5; 
  int binNr = (max - min) / binSize;
  binNr++;
  std::vector<int> hist(binNr);
  int tmp; 
  for (int i = 0; i < n; i++, p+=3){
    tmp = (*p - min) / binSize;
    hist[tmp]++; 
  }
  for (int i = 0;  i < binNr; i++)
    std::cout << hist[i]  << " ";
  std::cout << std::endl; 
}
// TODO: 
// possible arguments:
// 1. std::vector<Comp2Atom>
// 2. rvec *
// 3. std::vector <indcies> 
std::array<float, 6> analyseAtomsPos(rvec const * ps, int const nrPoints){
  // Filter Max, Min x,yz pos
  float max_x, min_x, max_y, min_y, max_z, min_z; 
  min_x = max_x = ps[0][0];
  min_y = max_y = ps[0][1];
  min_z = max_z = ps[0][2];
  for (int i = 1; i < nrPoints; i++){
    if (ps[i][0] < min_x) min_x = ps[i][0];    
    if (ps[i][0] > max_x) max_x = ps[i][0];    
    if (ps[i][1] < min_y) min_y = ps[i][1];    
    if (ps[i][1] > max_y) max_y = ps[i][1];    
    if (ps[i][2] < min_z) min_z = ps[i][2];    
    if (ps[i][2] > max_z) max_z = ps[i][2];    
  }
  std::cout << "x range: [" << min_x << ", " << max_x << "] "
	    << "y range: [" << min_y << ", " << max_y << "] "
	    << "z range: [" << min_z << ", " << max_z << "]" << std::endl; 
  std::array <float, 6> range;
  range[0] = floor(min_x); range[1] = ceil(max_x);
  range[2] = floor(min_y); range[3] = ceil(max_y);
  range[4] = floor(min_z); range[5] = ceil(max_z);
  // Histogram information
  calHist(min_x, max_x, (float *)ps, nrPoints );
  calHist(min_y, max_y, (float *)ps + 1, nrPoints );
  calHist(min_z, max_z, (float *)ps + 2, nrPoints );
   
  // divide the atoms based on z coordinate
  return range; 
}



void filterAtoms(ComplexType const & ctype, FilteredTypes & ftype){
  std::string usr_in = "";
  std::string token;
  // Let users to choose  atoms to be included in the analyses
  std::cout << "you have chosen to include " << ctype.name << "in the analyses.\n " 
            << ctype  
            << "Now choose the atoms. Enter the atom indcies and then press Enter"<< std::endl;
  std::getline(std::cin, usr_in);
  std::istringstream iss(usr_in);
  while (std::getline(iss, token, ' ')){
    ftype.atoms.push_back(std::stoi(token)); 
  }
}

std::shared_ptr<std::vector<FilteredTypes>> filterComplexs(GroData const & data){
  // Print out information on various complex types 
  int cntr_c=1;
  std::cout << "There are " << data.ctypes.size() << " types of lipids in the input data.\n\n";   
  for (auto & c: data.ctypes){
      std::cout << cntr_c << "." << c << std::endl;
      cntr_c++;
  }

  // Let users to choose complex types to be included in the analyses
  std::vector<int> ctypeIdxs; 
  std::string usr_in = "";
  cntr_c = 1; 
  for (auto & c: data.ctypes){
    std::cout << "Would you like to include lipid type " << c.name << " in the analyses? ((y)es/(n)o)" << std::endl;
    getline(std::cin, usr_in);
    if (usr_in[0] == 'Y' || usr_in[0] == 'y'){
       ctypeIdxs.push_back(cntr_c-1);
    }
    cntr_c++;
  }
   
  // Filter atoms per complex types
  auto ftyps = std::make_shared<std::vector<FilteredTypes>>() ; 
  for (auto & cIdx: ctypeIdxs){
    FilteredTypes fitem;
    fitem.typeIdx = cIdx;
    fitem.name = data.ctypes[cIdx].name;
    filterAtoms(data.ctypes[cIdx], fitem);
    ftyps->push_back(fitem);
  }
  
  return ftyps;
}



void testGroReader(GroData const & data, rvec const * ps){
  
  std::cout << "Various number of lipid types: " << data.ctypes.size() <<  std::endl;
  for (auto & c: data.ctypes)
    std::cout << c.name << " " << c.nrAtoms << " " << c.atoms[0]  << " " << c.atoms[c.nrAtoms - 1]  << std::endl;
 
  std::cout << "Total number of lipids: " << data.cs.size() <<  std::endl;
  for (int i = 450; i < 460; i++)
    std::cout << data.cs[i].rName << " " << data.cs[i].rNum << " " << data.cs[i].atom_idx << " " << data.cs[i].atom_no << std::endl;
  
  std::cout << "Total number of atoms: " << data.atoms.size() << std::endl;
  for (int i = 4640; i < 4651; i++)
    std::cout <<  data.atoms[i].name  << " " << data.atoms[i].num << " " 
              << ps [i][0] << " " << ps[i][1] << " " << ps[i][2] << std::endl;

}
