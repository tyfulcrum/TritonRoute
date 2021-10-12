#ifndef _FR_CORE_LANG_TYPES_H_
#define _FR_CORE_LANG_TYPES_H_
namespace coret {
  using frLayerNum = int;
  using frCoord = int;
  using frUInt4 = unsigned int;
  using frDist  = double;
  using frCost = unsigned int;
  using frMIdx = int; // negative value expected 
  enum class frDirEnum { UNKNOWN = 0, D = 1, S = 2, W = 3, E = 4, N = 5, U = 6 };
}

#endif
