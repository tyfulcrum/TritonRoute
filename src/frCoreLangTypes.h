#ifndef _FR_CORE_LANG_TYPES_H_
#define _FR_CORE_LANG_TYPES_H_

#include <climits>

namespace coret {
#define cuDIRBITSIZE 3
#define cuWAVEFRONTBUFFERSIZE 2
#define cuWAVEFRONTBITSIZE (cuWAVEFRONTBUFFERSIZE * cuDIRBITSIZE)
#define cuWAVEFRONTBUFFERHIGHMASK (111 << ((cuWAVEFRONTBUFFERSIZE - 1) * cuDIRBITSIZE))

  using frLayerNum = int;
  using frCoord = int;
  using frUInt4 = unsigned int;
  using frDist  = double;
  using frCost = unsigned int;
  using frMIdx = int; // negative value expected 
  enum class frDirEnum { UNKNOWN = 0, D = 1, S = 2, W = 3, E = 4, N = 5, U = 6 };
  enum frOrientEnum {
      frcR0       = 0, // N
      frcR90      = 1, // W
      frcR180     = 2, // S
      frcR270     = 3, // E
      frcMY       = 4, // FN
      frcMXR90    = 5, // FW
      frcMX       = 6, // FS
      frcMYR90    = 7  // FE
  };
  enum frPrefRoutingDirEnum {
      frcNotApplicablePrefRoutingDir = 0,
      frcNonePrefRoutingDir          = 1,
      frcHorzPrefRoutingDir          = 2,
      frcVertPrefRoutingDir          = 3 
  };
  class cuWavefrontGrid {
    public:
      cuWavefrontGrid(): xIdx(-1), yIdx(-1), zIdx(-1), pathCost(0), cost(0), layerPathArea(0), 
      vLengthX(INT_MAX), 
      vLengthY(INT_MAX), 
      dist(0), prevViaUp(false), 
      tLength(INT_MAX), backTraceBuffer() {}
      cuWavefrontGrid(int xIn, int yIn, int zIn, frCoord layerPathAreaIn, 
          frCoord vLengthXIn, frCoord vLengthYIn,
          bool prevViaUpIn, frCoord tLengthIn,
          frCoord distIn, frCost pathCostIn, frCost costIn/*, frDirEnum preTurnDirIn*/): 
        xIdx(xIn), yIdx(yIn), zIdx(zIn), pathCost(pathCostIn), cost(costIn), 
        layerPathArea(layerPathAreaIn), vLengthX(vLengthXIn), vLengthY(vLengthYIn),
        dist(distIn), prevViaUp(prevViaUpIn), tLength(tLengthIn), backTraceBuffer() {}
      cuWavefrontGrid(int xIn, int yIn, int zIn, frCoord layerPathAreaIn, 
          frCoord vLengthXIn, frCoord vLengthYIn,
          bool prevViaUpIn, frCoord tLengthIn,
          frCoord distIn, frCost pathCostIn, frCost costIn, 
          unsigned int backTraceBufferIn): 
        xIdx(xIn), yIdx(yIn), zIdx(zIn), pathCost(pathCostIn), cost(costIn), 
        layerPathArea(layerPathAreaIn), vLengthX(vLengthXIn), vLengthY(vLengthYIn),
        dist(distIn), prevViaUp(prevViaUpIn), tLength(tLengthIn), backTraceBuffer(backTraceBufferIn) {}
      bool operator<(const cuWavefrontGrid &b) const {
        if (this->cost != b.cost) {
          return this->cost > b.cost; // prefer smaller cost
        } else {
          if (this->dist != b.dist) {
            return this->dist > b.dist; // prefer routing close to pin gravity center (centerPt)
          } else {
            if (this->zIdx != b.zIdx) {
              return this->zIdx < b.zIdx; // prefer upper layer
            } else {
              return this->pathCost < b.pathCost; //prefer larger pathcost, DFS-style
            }
          }
        }
      }
      // getters
      frMIdx x() const {
        return xIdx;
      }
      frMIdx y() const {
        return yIdx;
      }
      frMIdx z() const {
        return zIdx;
      }
      frCost getPathCost() const {
        return pathCost;
      }
      frCost getCost() const {
        return cost;
      }
      auto getBackTraceBuffer() const {
        return backTraceBuffer;
      }
      frCoord getLayerPathArea() const {
        return layerPathArea;
      }
      frCoord getLength() const {
        return vLengthX;
      }
      void getVLength(frCoord &vLengthXIn, frCoord &vLengthYIn) const {
        vLengthXIn = vLengthX;
        vLengthYIn = vLengthY;
      }
      bool isPrevViaUp() const {
        return prevViaUp;
      }
      frCoord getTLength() const {
        return tLength;
      }
      // setters
      void addLayerPathArea(frCoord in) {
        layerPathArea += in;
      }
      void resetLayerPathArea() {
        layerPathArea = 0;
      }

      void resetLength() {
        vLengthX = 0;
        vLengthY = 0;
      }
      void setPrevViaUp(bool in) {
        prevViaUp = in;
      }
      frDirEnum getLastDir() const {
        auto currDirVal = backTraceBuffer & 0b111u;
        return static_cast<frDirEnum>(currDirVal);
      }
      bool isBufferFull() const {
        unsigned int mask = cuWAVEFRONTBUFFERHIGHMASK;
        return any(mask & backTraceBuffer);
      }
      frDirEnum shiftAddBuffer(const frDirEnum &dir) {
        auto retBS = static_cast<frDirEnum>((backTraceBuffer >> (cuWAVEFRONTBITSIZE - cuDIRBITSIZE)));
        backTraceBuffer <<= cuDIRBITSIZE;
        auto newBS = (unsigned)dir;
        backTraceBuffer |= newBS;
        return retBS;
      }
    protected:
      frMIdx xIdx, yIdx, zIdx;
      frCost pathCost; // path cost
      frCost cost; // path + est cost
      frCoord layerPathArea;
      frCoord vLengthX;
      frCoord vLengthY;
      frCoord dist; // to maze center
      bool    prevViaUp;
      frCoord tLength; // length since last turn
      unsigned int backTraceBuffer;
    private:
      bool any(unsigned int num) const{
        auto size = sizeof(unsigned int);
        auto maxPow = 1<<(size*8-1);
        int i=0;
        for(;i<size;++i){
          for(;i<size*8;++i){
            // print last bit and shift left.
            auto bit = num & maxPow;
            if (bit) {
              return true;
            }
            num = num<<1;
          }       
        }
        return false;
      }
  };
}

#endif
