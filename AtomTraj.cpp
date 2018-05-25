#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>


class AtomTraj{
private:
  struct STPoint{
    float x;
    float y;
    float z;
    float t;
  
    friend std::ostream & operator << (std::ostream & os, STPoint const & p){
      os << p.x << " " << p.y << " " << p.t << std::endl; 
    }
 
  };

  static const float BOX_MARGIN;  

  unsigned int frameNu;
  unsigned int atomNu;
  unsigned int frameIdx;  
  std::vector <std::vector <STPoint> > trajs; 
  std::vector <std::vector <int> > jumps; 
public:
  AtomTraj(unsigned int fnu, unsigned int aNu)
    : trajs(std::vector<std::vector <STPoint> >( aNu, std::vector<STPoint>(fnu))), frameNu(fnu), atomNu(aNu)
    , jumps (std::vector <std::vector <int>> (aNu) )
    , frameIdx(0) {}
 
  bool isJump(unsigned int fIdx, unsigned int aIdx, const float (&box)[3][3]){
    if (fIdx == 0) return false;

    float bxl = 0, bxr = box[0][0]; 
    float byl = 0, byr = box[1][1]; 

    float px = trajs[aIdx][fIdx - 1].x;
    float x = trajs[aIdx][fIdx].x;

    if (x > bxr - BOX_MARGIN && px < bxl + BOX_MARGIN ) // crossing the left vertical border
      return true;
    if (x < bxl + BOX_MARGIN && px > bxr - BOX_MARGIN ) // crossing the right vertical border
      return true;

    float py = trajs[aIdx][fIdx - 1].y;
    float y = trajs[aIdx][fIdx].y;
    if (y > byr - BOX_MARGIN && py < byl + BOX_MARGIN ) // crossing the bottom horizontal  border
      return true;
    if (y < byl + BOX_MARGIN && py > byr - BOX_MARGIN ) // crossing the top horizontal border
      return true;

    return false; 
  }

  void addFrame(const float * ps[3], std::vector <unsigned int> const * as , const float (&box) [3][3]){
    for (auto  a:*as){
      // extract position
      trajs[a][frameIdx].x = ps[a][0]; 
      trajs[a][frameIdx].y = ps[a][1]; 
      trajs[a][frameIdx].z = ps[a][2]; 
      if (isJump(frameIdx, a, box))
        jumps[a].push_back(frameIdx);
    }
    frameIdx++; 
  }

  void printSubTraj( unsigned int aIdx, unsigned int tIdx_start, unsigned int tIdx_end
                , std::ostringstream & buf){
     for (unsigned i = tIdx_start; i < tIdx_end; i++){
       buf << trajs[aIdx][i];
     }
     buf << "0 0 0" << std::endl;  // end of a trajectory

  }

  void printAtomTraj(int aIdx, std::ostringstream & buf){   
    // Iterate over atom positions and add them to the 'trajectory' 
    // When the atom pass box boundries, create a new trajectory
    auto & ajs = jumps[aIdx];
    auto & ats = trajs[aIdx];
    unsigned int tIdx = 0;  
    for (unsigned short j = 0; j < ajs.size(); j++){ 
      // Each atom trajectory is dividend to a number of trajectories
      // each time the atom passes the box boundaries, a new trajectory starts 
      printSubTraj(aIdx, tIdx, ajs[j], buf); 
      tIdx = ajs[j];  
    }
    printSubTraj(aIdx, tIdx, frameNu, buf);

  }

  void writeTrajs(const char * fn){
    std::ofstream trajFl; 
    trajFl.open (fn, std::ofstream::out);   
    std::ostringstream buffer (std::ios_base::out); 
    // iterate over atoms and dump trajectories
    for (unsigned int i = 0; i < atomNu; i++){
      printAtomTraj(i, buffer);
      trajFl << buffer; 
      buffer.str("");
    }
 
    trajFl.close();
  }

};

const float AtomTraj::BOX_MARGIN = 0.2; 
