#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>

#include <iostream>
#include <istream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <vector>
#include <array>
#include <map>
#include <memory>
#include<utility>

#include <cstring>
#include <string>
#include <cmath>

#define MAX_GRO_LINE_LENGTH 69 // %5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f 
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
  Complex(int num, std::string const & name, int a_idx, int aNr)
     : rNum(num), atom_idx(a_idx), rName(name), atom_no(aNr) {
     }
  int rNum;           // residue number 
  std::string rName;  // residue name 
  int atom_idx;       // Index of the first atom in atoms list 
  int atom_no;        // Number of atoms composing the complex 
}; 


struct Protein{
  Protein(int fai, int lai, int fmi, int lmi)
    : fstAtomIdx(fai), lstAtomIdx(lai), fstMolIdx(fmi), lstMolIdx(lmi){};
  int fstAtomIdx; 
  int lstAtomIdx; 
  int fstMolIdx; 
  int lstMolIdx; 
};

// This struct represetn an complex type, e.g. POPC
struct ComplexType{
  ComplexType(std::string n, int aNr, std::vector<std::string> & as) : name(n), nrAtoms(aNr), atoms(as){}; 
  std::string name;                  // Complex type name
  unsigned short nrAtoms;            // Number of atoms
  std::vector<std::string> atoms;    // list of atoms 

  friend std::ostream & operator<<(std::ostream & os, ComplexType const & ctype);
};

struct GroData {
  GroData(int n): atoms(n){}
  std::vector<ComplexType> ctypes; // list of Complex types 	
  std::vector<ComplexType> aatypes; // list of Complex types 	
  std::vector <Atom> atoms;        // array of Atoms 
  std::vector <Complex> lipids;        // array of Complexes
  std::vector <Protein> ps; 
  std::vector <Complex> aminoAcids; 

  // std::map<std::string, std::vector<std::string>> a2ctype;  TODO do we need this?
  bool hasLipidType(std::string name){
    for (auto it = ctypes.crbegin(); it < ctypes.crend(); it++){
      if ((*it).name == name)
       return true; 
    }
    return false; 
  }

  bool hasAminoAcid(std::string name){
    for (auto & c : aatypes){
      if (c.name == name)
       return true;
    }
    return false;
  }
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

class ReadGro {
private:
  struct GroLine {
    int cNu;
    int aNu;
    std::string cName;
    std::string aName; 
  };

  void readLine(GroLine & gl){
    char line [MAX_GRO_LINE_LENGTH];
    char * lptr; 
    std::string item;  
    int tmp;
    // Read the new line/atom 
    gf.read(line,MAX_GRO_LINE_LENGTH);
    lptr = line; 
    // Extract residue number       
    item.assign(lptr, 5);  
    lptr += 5; 
    gl.cNu = std::stoi(item);
    // Extract residue name
    tmp = 4;  
    while (lptr[tmp] == ' ') tmp--; 
    gl.cName.assign(lptr, tmp+1);
    lptr += 5;
    // Extract atom name 
    for (tmp = 5 ; *lptr == ' '; tmp--) lptr++;
    gl.aName.assign(lptr, tmp);
    lptr += tmp; 
    // Extract atom Number 
    item.assign(lptr, 5);
    lptr += 5;
    gl.aNu = std::stoi(item);
    // record the atom position
    item.assign(lptr, 8);
    lptr+=8;
    pos[lIdx][0] = std::stof(item);
    item.assign(lptr, 8);
    lptr+=8;
    pos[lIdx][1] = std::stof(item);
    item.assign(lptr, 8);
    pos[lIdx][2] = std::stof(item);
    
    lIdx ++;
  }
  
  std::ifstream gf; 
  int lIdx = 0; 
  rvec * pos; 
  void readComplex(GroData & groData, int nr);
  void readProtein(GroData & groData, int nr);
public:
  ReadGro(): lIdx(0){}
  void readGroFile(const char * fn, GroData & groData, rvec * ps, std::vector<std::pair<int, int>> pns);
};



void ReadGro::readGroFile(const char * fn, GroData & groData, rvec * ps, std::vector<std::pair<int, int>> pns){
  pos = ps; 
  char line [MAX_GRO_LINE_LENGTH];
  char * token;

  gf.open(fn, std::fstream::in);

  // Read the first lin,title string, ignore it for now 
  gf.ignore(50, '\n');     

  // Read number of atoms 
  gf.getline(line, MAX_GRO_LINE_LENGTH); 
  token = strtok(line, " ");
  int nr = std::atoi(token); 

  int pidx = 0, tmp; 
  bool checkProtein = pns.size() > 0 ? true :false; 
  while (lIdx < nr){
    // if it belongs to a Protein call readProtein 
    if (checkProtein){
      tmp = std::get<0>(pns[pidx]);
      if (lIdx + 1 == tmp){
        tmp = std::get<1>(pns[pidx]) - tmp + 1; 
        readProtein(groData, tmp);   // Read protein for the next tmp lines 
        pidx++; 
        checkProtein = pns.size() > pidx ? true: false;        
      } else{
         //TODO wrong logic 
         tmp = tmp - 1 -  lIdx;
         readComplex(groData, tmp);  //Read lipids for tmp lines 
      }
    }else {
      tmp  = nr - lIdx; 
      readComplex(groData, tmp);  //Read lipids for tmp lines 
    } 
  }
  gf.close();
}



void ReadGro::readComplex(GroData & groData, int nr){

  // Iterate over all lines in the file (one atom per line)
  // and read atoms information 
  GroLine ndata;

  std::vector<std::string> as; 
  //Save ldata info
  int end = lIdx + nr;
  readLine(ndata); 
  std::strcpy(groData.atoms[lIdx -1 ].name, ndata.aName.c_str()); 
  groData.atoms[lIdx - 1].num = ndata.aNu;
  bool nc;
  int aCntr = 1;
  int cNuPrev = ndata.cNu;
  std::string cNamePrev = ndata.cName; 
  bool nct = groData.hasLipidType(ndata.cName) ? false: true;  // Is it a new complex type? ;
  if (nct) as.push_back(ndata.aName);
  for (int i = lIdx ; i < end; i++){
    // Read the new line/atom 
    readLine(ndata); 
    // save atom info in the atom table
    std::strcpy(groData.atoms[i].name, ndata.aName.c_str()); 
    groData.atoms[i].num = ndata.aNu;
    
    nc = (cNuPrev == ndata.cNu) ? false: true;    // Is it a new complex?  
    if (nc){
      // if it is a new complex,
      //  1. save prev complex info
      groData.lipids.push_back(Complex(cNuPrev, cNamePrev, i - aCntr, aCntr));        
      //  2. if the prev complex was a new type, save the type 
      if (nct) {
        groData.ctypes.push_back(ComplexType(cNamePrev, aCntr, as));
        as.clear();
      }
      //  3. check whether the new complex has a new type
      aCntr = 0;
      cNuPrev = ndata.cNu; 
      cNamePrev = ndata.cName; 
      nct = groData.hasLipidType(ndata.cName) ? false: true;  // Is it a new complex type? 
    } 
    aCntr++;
    // if this atom belongs to a new complex type, save the atom name
    if (nct) as.push_back(ndata.aName);
  }
  
  groData.lipids.push_back(Complex(cNuPrev, cNamePrev, end - aCntr, aCntr));  // save the last complex  
  if (nct)
     groData.ctypes.push_back(ComplexType(cNamePrev, aCntr, as)); // save the last complex type  

}


void ReadGro::readProtein(GroData & groData, int nr){

  int fstAtomIdx = lIdx; 
  int fstMolIdx = groData.aminoAcids.size(); 
  std::vector<std::string> as; 
  //Save ldata info
  int end = lIdx + nr;
  GroLine ndata;
  readLine(ndata); 
  std::strcpy(groData.atoms[lIdx -1 ].name, ndata.aName.c_str()); 
  groData.atoms[lIdx - 1].num = ndata.aNu;
  bool nc;
  int aCntr = 1;
  int cNuPrev = ndata.cNu;
  std::string cNamePrev = ndata.cName; 
  bool nct = groData.hasAminoAcid(ndata.cName) ? false: true;  // Is it a new complex type? ;
  if (nct) as.push_back(ndata.aName);
  // Iterate over all lines in the file (one atom per line)
  // and read atoms information 
  for (int i = lIdx ; i < end; i++){
    // Read the new line/atom 
    readLine(ndata); 
    // save atom info in the atom table
    std::strcpy(groData.atoms[i].name, ndata.aName.c_str()); 
    groData.atoms[i].num = ndata.aNu;
    
    nc = (cNuPrev == ndata.cNu) ? false: true;    // Is it a new complex?  
    if (nc){
      // if it is a new complex,
      //  1. save prev complex info
      groData.aminoAcids.push_back(Complex(cNuPrev, cNamePrev, i - aCntr, aCntr));        
      //  2. if the prev complex was a new type, save the type 
      if (nct) {
        groData.aatypes.push_back(ComplexType(cNamePrev, aCntr, as));
        as.clear();
      }
      //  3. check whether the new complex has a new type
      aCntr = 0;
      cNuPrev = ndata.cNu; 
      cNamePrev = ndata.cName; 
      nct = groData.hasAminoAcid(ndata.cName) ? false: true;  // Is it a new complex type? 
    } 
    aCntr++;
    // if this atom belongs to a new complex type, save the atom name
    if (nct) as.push_back(ndata.aName);
  }
  
  groData.aminoAcids.push_back(Complex(cNuPrev, cNamePrev, end - aCntr, aCntr));  // save the last complex  
  if (nct)
     groData.aatypes.push_back(ComplexType(cNamePrev, aCntr, as)); // save the last complex type 
  
  
  groData.ps.push_back(Protein(fstAtomIdx
                              , lIdx - 1 
                              , fstMolIdx
                              , groData.aminoAcids.size() - 1)); 

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
  
  std::cout << "Number of atoms: " << data.atoms.size() 
            << "\nNumber of non amino acid molecules:" << data.lipids.size() 
            << "\nNumber of amino acids:  " << data.aminoAcids.size()  
            << "\nNumber of Proteins:  " << data.ps.size() 
            << "\nTypes of amino acids: "<< data.aatypes.size() 
            << "\nTypes of non amino acid molecules " << data.ctypes.size() << "\n" ;
  
  std::cout << "Proteins: " <<  std::endl;
  for (auto & p : data.ps)
    std::cout << "Protein first & last atom idx: (" << p.fstAtomIdx << ", " << p.lstAtomIdx << ") Protein first & last Mol idx:(" << p.fstMolIdx << ", " << p.lstMolIdx << ")\n";
  
  std::cout << "Amino acids: " << std::endl;
  for (auto & c: data.aatypes)
    std::cout << c.name << " " << c.nrAtoms << " " << c.atoms[0]  << " " << c.atoms[c.nrAtoms - 1]  << std::endl;

  std::cout << "Various number of lipid types: " << data.ctypes.size() <<  std::endl;
  for (auto & c: data.ctypes)
    std::cout << c.name << " " << c.nrAtoms << " " << c.atoms[0]  << " " << c.atoms[c.nrAtoms - 1]  << std::endl;
 
  std::cout << "amino acids 50 to 59" << std::endl;
  for (int i = 50; i < 60; i++)
    std::cout << data.aminoAcids[i].rName << " " << data.aminoAcids[i].rNum << " " << data.aminoAcids[i].atom_idx << " " << data.aminoAcids[i].atom_no << std::endl;
  
  std::cout << "Lipids 450 to 459" << std::endl;
  for (int i = 450; i < 460; i++)
    std::cout << data.lipids[i].rName << " " << data.lipids[i].rNum << " " << data.lipids[i].atom_idx << " " << data.lipids[i].atom_no << std::endl;
  
  std::cout << "Total number of atoms: " << data.atoms.size() << std::endl;
  for (int i = 4640; i < 4651; i++)
    std::cout <<  data.atoms[i].name  << " " << data.atoms[i].num << " " 
              << ps [i][0] << " " << ps[i][1] << " " << ps[i][2] << std::endl;

}
