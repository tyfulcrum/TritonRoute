/* Authors: Lutong Wang and Bangqi Xu */
/*
 * Copyright (c) 2019, The Regents of the University of California
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *     * Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *     * Neither the name of the University nor the
 *       names of its contributors may be used to endorse or promote products
 *       derived from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE REGENTS BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <memory>
#include <utility>
#include "dr/FlexGridGraph.h"
#include "dr/FlexDR.h"
#include "pathway/GPU-solver.hpp"
#include <fmt/core.h>

using namespace std;
using namespace fr;


void FlexGridGraph::expand(FlexWavefrontGrid &currGrid, const frDirEnum &dir, 
                                      const FlexMazeIdx &dstMazeIdx1, const FlexMazeIdx &dstMazeIdx2,
                                      const frPoint &centerPt) {
  bool enableOutput = false;
  //bool enableOutput = true;
  frCost nextEstCost, nextPathCost;
  int gridX = currGrid.x();
  int gridY = currGrid.y();
  int gridZ = currGrid.z();

  getNextGrid(gridX, gridY, gridZ, dir);
  
  FlexMazeIdx nextIdx(gridX, gridY, gridZ);
  // get cost
  nextEstCost = getEstCost(nextIdx, dstMazeIdx1, dstMazeIdx2, dir);
  nextPathCost = getNextPathCost(currGrid, dir, false);  
  if (enableOutput) {
    std::cout << "  expanding from (" << currGrid.x() << ", " << currGrid.y() << ", " << currGrid.z() 
              << ") [pathCost / totalCost = " << currGrid.getPathCost() << " / " << currGrid.getCost() << "] to "
              << "(" << gridX << ", " << gridY << ", " << gridZ << ") [pathCost / totalCost = " 
              << nextPathCost << " / " << nextPathCost + nextEstCost << "]\n";
  }
  auto lNum = getLayerNum(currGrid.z());
  auto pathWidth = getDesign()->getTech()->getLayer(lNum)->getWidth();
  frPoint currPt;
  getPoint(currPt, gridX, gridY);
  frCoord currDist = abs(currPt.x() - centerPt.x()) + abs(currPt.y() - centerPt.y());

  // vlength calculation
  frCoord currVLengthX = 0;
  frCoord currVLengthY = 0;
  currGrid.getVLength(currVLengthX, currVLengthY);
  auto nextVLengthX = currVLengthX;
  auto nextVLengthY = currVLengthY;
  bool nextIsPrevViaUp = currGrid.isPrevViaUp();
  if (dir == frDirEnum::U || dir == frDirEnum::D) {
    nextVLengthX = 0;
    nextVLengthY = 0;
    nextIsPrevViaUp = (dir == frDirEnum::D); // up via if current path goes down
  } else {
    if (currVLengthX != std::numeric_limits<frCoord>::max() &&
        currVLengthY != std::numeric_limits<frCoord>::max()) {
      if (dir == frDirEnum::W || dir == frDirEnum::E) {
        nextVLengthX += getEdgeLength(currGrid.x(), currGrid.y(), currGrid.z(), dir);
      } else { 
        nextVLengthY += getEdgeLength(currGrid.x(), currGrid.y(), currGrid.z(), dir);
      }
    }
  }
  
  // tlength calculation
  auto currTLength = currGrid.getTLength();
  auto nextTLength = currTLength;
  // if there was a turn, then add tlength
  if (currTLength != std::numeric_limits<frCoord>::max()) {
    nextTLength += getEdgeLength(currGrid.x(), currGrid.y(), currGrid.z(), dir);
  }
  // if current is a turn, then reset tlength
  if (currGrid.getLastDir() != frDirEnum::UNKNOWN && currGrid.getLastDir() != dir) {
    nextTLength = getEdgeLength(currGrid.x(), currGrid.y(), currGrid.z(), dir);
  }
  // if current is a via, then reset tlength
  if (dir == frDirEnum::U || dir == frDirEnum::D) {
    nextTLength = std::numeric_limits<frCoord>::max();
  }

  FlexWavefrontGrid nextWavefrontGrid(gridX, gridY, gridZ, 
                                      currGrid.getLayerPathArea() + getEdgeLength(currGrid.x(), currGrid.y(), currGrid.z(), dir) * pathWidth, 
                                      nextVLengthX, nextVLengthY, nextIsPrevViaUp,
                                      nextTLength,
                                      currDist,
                                      nextPathCost, nextPathCost + nextEstCost, currGrid.getBackTraceBuffer());
  if (dir == frDirEnum::U || dir == frDirEnum::D) {
    nextWavefrontGrid.resetLayerPathArea();
    nextWavefrontGrid.resetLength();
    if (dir == frDirEnum::U) {
      nextWavefrontGrid.setPrevViaUp(false);
    } else {
      nextWavefrontGrid.setPrevViaUp(true);
    }
    nextWavefrontGrid.addLayerPathArea((dir == frDirEnum::U) ? getHalfViaEncArea(currGrid.z(), false) : getHalfViaEncArea(gridZ, true));
  }
  // update wavefront buffer
  auto tailDir = nextWavefrontGrid.shiftAddBuffer(dir);
  // non-buffer enablement is faster for ripup all
  // commit grid prev direction if needed
  auto tailIdx = getTailIdx(nextIdx, nextWavefrontGrid);
  if (tailDir != frDirEnum::UNKNOWN) {
    if (getPrevAstarNodeDir(tailIdx.x(), tailIdx.y(), tailIdx.z()) == frDirEnum::UNKNOWN ||
        getPrevAstarNodeDir(tailIdx.x(), tailIdx.y(), tailIdx.z()) == tailDir) {
      setPrevAstarNodeDir(tailIdx.x(), tailIdx.y(), tailIdx.z(), tailDir);
      wavefront.push(nextWavefrontGrid);
      if (enableOutput) {
        std::cout << "    commit (" << tailIdx.x() << ", " << tailIdx.y() << ", " << tailIdx.z() << ") prev accessing dir = " << (int)tailDir << "\n";
      }
    }
  } else {  
    // add to wavefront
    wavefront.push(nextWavefrontGrid);
  }

  return;
}

void FlexGridGraph::expandWavefront(FlexWavefrontGrid &currGrid, const FlexMazeIdx &dstMazeIdx1, 
                                               const FlexMazeIdx &dstMazeIdx2, const frPoint &centerPt) {
  bool enableOutput = false;
  //bool enableOutput = true;
  if (enableOutput) {
    cout << "start expand from (" << currGrid.x() << ", " << currGrid.y() << ", " << currGrid.z() << ")\n";
  }
  //if (currGrid.y() == 19 && currGrid.z() == 0) {
  //  cout <<"is expandable (" <<currGrid.x() <<", " <<currGrid.y() <<", " <<currGrid.z() <<") NESWUD "
  //       <<isExpandable(currGrid, frDirEnum::N)
  //       <<isExpandable(currGrid, frDirEnum::E)
  //       <<isExpandable(currGrid, frDirEnum::S)
  //       <<isExpandable(currGrid, frDirEnum::W)
  //       <<isExpandable(currGrid, frDirEnum::U)
  //       <<isExpandable(currGrid, frDirEnum::D)
  //       <<endl;
  //  cout <<"has edge " 
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::N)
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::E)
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::S)
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::W)
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::U)
  //       <<gridGraph.hasEdge(currGrid.x(), currGrid.y(), currGrid.z(), frDirEnum::D)
  //       <<endl;
  //  int gridX = currGrid.x();
  //  int gridY = currGrid.y();
  //  int gridZ = currGrid.z();
  //  if (!gridGraph.hasEdge(gridX, gridY, gridZ, frDirEnum::E)) {
  //    ;
  //  } else {
  //    getNextGrid(gridX, gridY, gridZ, frDirEnum::E);
  //    if (gridGraph.isBlocked(gridX, gridY, gridZ)) {
  //      cout <<"blocked" <<endl;
  //    } else if (gridGraph.isSrc(gridX, gridY, gridZ)) {
  //      cout <<"src" <<endl;
  //    } else if (gridGraph.getPrevAstarNodeDir(gridX, gridY, gridZ) != frDirEnum::UNKNOWN) {
  //      cout <<"visited" <<endl;
  //    } else {
  //      ;
  //    }
  //  }
  //}
  //auto tmpGrid = currWavefrontGrid;
  //// commit grid prev direction if needed
  //auto tailIdx = getTailIdx(nextIdx, nextWavefrontGrid);
  //if (tailDir != frDirEnum::UNKNOWN) {
  //  if (getPrevAstarNodeDir(tailIdx.x(), tailIdx.y(), tailIdx.z()) == frDirEnum::UNKNOWN) {
  //    ;
  //  }
  //}
  // N
  // auto gpuSolver = GPUPathwaySolver(&currGrid, this);
  // gpuSolver.initialize(bits);
  // fmt::print("\nCPU knows bits[0]: {}\n", bits[0]);
  
  /*
  bool dirn = isExpandable(currGrid, frDirEnum::N);
  bool cuDirn = gpuSolver.cuIsExpanable(frDirEnum::N);
  */










  vector<frUInt4> path_widths;
  auto btm = getDesign()->getTech()->getBottomLayerNum();
  auto top = getDesign()->getTech()->getTopLayerNum();
  for (int i = btm; i <= top; ++i) {
    auto pathWidth = getDesign()->getTech()->getLayer(i)->getWidth();
    path_widths.push_back(pathWidth);
  }

  bool drWorker_ava = false;
  int drIter = -1;
  int ripupMode = -1;

  if (drWorker) {
    drIter = drWorker->getDRIter();
    ripupMode = drWorker->getRipupMode();
  }

  // int *a = nullptr;
  // a[32] = 3243;
  frMIdx gridX = currGrid.x();
  frMIdx gridY = currGrid.y();
  frMIdx gridZ = currGrid.z();
  auto dir = frDirEnum::E;
  getNextGrid(gridX, gridY, gridZ, dir);
  
  FlexMazeIdx nextIdx(gridX, gridY, gridZ);
  //  [9][8] 
  vector<int> via2ViaForbOverlapLen;
  viaData(via2ViaForbOverlapLen, getTech()->getVia2ViaForbiddenOverlapLen());
  // [9][8] 
  vector<int> via2viaForbLen;
  viaData(via2viaForbLen, getTech()->getVia2ViaForbiddenLen());
  // [9][4] 
  vector<int> viaForbiTurnLen;
  viaData(viaForbiTurnLen, getTech()->getViaForbiddenTurnLen());

  frCoord cuGrid_vLengthX, cuGrid_vLengthY;
  currGrid.getVLength(cuGrid_vLengthX, cuGrid_vLengthY);

  frPoint currPt;
  getPoint(currPt, gridX, gridY);
  frCoord currDist = abs(currPt.x() - centerPt.x()) + abs(currPt.y() - centerPt.y());

  auto estc = getNextPathCost(currGrid, dir, true);
  auto cuGrid = cuWavefrontGrid(currGrid.x(), currGrid.y(), currGrid.z(), 
      currGrid.getLayerPathArea(), cuGrid_vLengthX, cuGrid_vLengthY, 
      currGrid.isPrevViaUp(), currGrid.getTLength(), currDist, currGrid.getPathCost(), 
      currGrid.getCost(),
      currGrid.getBackTraceBuffer().to_ulong());
  auto ddir =  currGrid.getLastDir();
  auto hdir =  cuGrid.getLastDir();

    int p_viaFOLen_size = via2ViaForbOverlapLen.size();
    int p_viaFLen_size = via2viaForbLen.size();

    int p_viaFTLen_size = (getTech()->getViaForbiddenTurnLen())[0].size();
    // fmt::print("viaForbiTurnLen size: {}\n", p_viaFTLen_size);
  //if (hdir != frDirEnum::N && hdir != frDirEnum::S) {
  // fmt::print("CPU dir: {}, GPU dir: {}\n", hdir, ddir);
  auto gpu = GPUPathwaySolver();
  gpu.initialize(bits, prevDirs, srcs, guides, zDirs, xCoords, yCoords, zCoords,
      zHeights, path_widths, ggDRCCost, ggMarkerCost, 
      via2ViaForbOverlapLen, via2viaForbLen, viaForbiTurnLen, 
      drWorker, drIter, ripupMode, 
    p_viaFOLen_size, p_viaFLen_size, p_viaFTLen_size);
  /*
  auto d_estc = gpu.test_npCost(dir, currGrid.x(), currGrid.y(), currGrid.z(), 
      currGrid.getLayerPathArea(), cuGrid_vLengthX, cuGrid_vLengthY, 
      currGrid.isPrevViaUp(), currGrid.getTLength(), currDist, currGrid.getPathCost(), 
      currGrid.getCost(), 
      currGrid.getBackTraceBuffer().to_ulong());
      */
  // auto d_estc = gpu.test_expand(frDirEnum::D, grid


  // if (d_estc == estc) {
  //   ++gpu_pass;
  //   // fmt::print("SUCCESS in ({}, {}, {}) CPU: {} GPU EstCost: {}\n", tx, ty, tz, estc, d_estc);
  // } else {
  //   // ++gpu_failed;
  //   // fmt::print("GPU FAIL!\n");
  // }
  /*
  else if (d_estc > estc) {
    fmt::print("FAILED GREATER in ({}, {}, {}) CPU: {} GPU EstCost: {}\n", tx, ty, tz, estc, d_estc);
  } else {
    fmt::print("FAILED LESS in ({}, {}, {}) CPU: {} GPU EstCost: {}\n", tx, ty, tz, estc, d_estc);
  }
  auto d_estc = gpu.test_estcost(nextIdx, dstMazeIdx1, dstMazeIdx2, dir);
  auto estc = getEstCost(nextIdx, dstMazeIdx1, dstMazeIdx2, dir);
  if (d_estc != estc) {
    fmt::print("GPU EstCost: {} in ({}, {}, {}) FAILED!\n", d_estc,  gridX, gridY, gridZ);
  }
  if (gridX == 36 && gridY == 31 && gridZ == 1) {
    auto gpu = GPUPathwaySolver();
    gpu.initialize(bits, prevDirs, srcs, guides, zDirs, xCoords, yCoords, zCoords,
        zHeights, xCoords.size(), yCoords.size(), zCoords.size());
    auto d_estc = gpu.test_estcost(nextIdx, dstMazeIdx1, dstMazeIdx2, dir);
    auto estc = getEstCost(nextIdx, dstMazeIdx1, dstMazeIdx2, dir);
    if (d_estc == estc) {
      fmt::print("GPU EstCost in ({}, {}, {}) SUCCESS!\n", gridX, gridY, gridZ);
    }
    fmt::print("GPU EstCost result: {}\n", d_estc);
    fmt::print("CPU EstCost result: {}\n", estc);
    auto dir = currGrid.getLastDir();
    auto gpures = gpu.isEx(gridX, gridY, gridZ, dir);
    auto cpures = isExpandable(currGrid, frDirEnum::S);
    fmt::print("GPU isExpandable ({}, {}, {}) S result: {}\n", gridX, gridY, gridZ, gpures);
    fmt::print("CPU isExpandable ({}, {}, {}) S result: {}\n", gridX, gridY, gridZ,  cpures);
  }
  */
  dir = frDirEnum::N;
  /*
  auto cpu_isex = isExpandable(currGrid, dir);
  auto gpu_isex = gpu.isEx(gridX, gridY, gridZ, dir);
  */
  auto const all_dirs = {frDirEnum::D, frDirEnum::E, frDirEnum::N, frDirEnum::S, 
    frDirEnum::U, frDirEnum::W};
  auto tx = currGrid.x();
  auto ty = currGrid.y();
  auto tz = currGrid.z();
  auto const lastdir = currGrid.getLastDir();
  for (auto dir: all_dirs) {
    auto const gpu_value = gpu.isEx(tx, ty, tz, dir, lastdir);
    auto const cpu_value = isExpandable(currGrid, dir);
    std::map<frDirEnum, string> const dir_str = {
      {frDirEnum::E, "East"}, 
        {frDirEnum::W, "West"}, 
        {frDirEnum::N, "North"}, 
        {frDirEnum::S, "South"}, 
        {frDirEnum::U, "Upper"}, 
        {frDirEnum::D, "Down"}, 
        {frDirEnum::UNKNOWN, "Unknow"}}; 
    if (cpu_value == gpu_value) {
      ++gpu_pass;
    } else {
      ++gpu_failed;
      fmt::print("==========================\n");
      fmt::print("At: {}, {}, {}, {}\n", tx, ty, tz, dir_str.at(dir));
      // fmt::print("GPU result in ({}, {}, {}, {}) FAILED!\n", gridX, gridY, gridZ, dir);
    }
    /*
       auto cpu_value = isExpandable(currGrid, dir);
       auto gpu_value = gpu.isEx(gridX, gridY, gridZ, dir, currGrid.getLastDir());
       if (cpu_value == gpu_value) {
       ++gpu_pass;
       } else {
       ++gpu_failed;
       fmt::print("GPU result in ({}, {}, {}, {}) FAILED!\n", gridX, gridY, gridZ, dir);
       }
       */
  }
  if (isExpandable(currGrid, frDirEnum::N)) {
    //cout << "(" << currGrid.x() << "," << currGrid.y() + 1 << "," << currGrid.z() << ") expendable! " << endl;
    expand(currGrid, frDirEnum::N, dstMazeIdx1, dstMazeIdx2, centerPt);
    // gpu.test_expand(currGrid, frDirEnum::N, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no N" <<endl;
  // }
  // E
  // dir = frDirEnum::E;
  // cpures = isExpandable(currGrid, dir);
  // gpures = gpu.isEx(gridX, gridY, gridZ, dir);
  // if (cpures == gpures) {
  //   fmt::print("GPU result in ({}, {}, {}) FAILED!\n", gridX, gridY, gridZ);
  // }
  if (isExpandable(currGrid, frDirEnum::E)) {
   // cout << "(" << currGrid.x() + 1 << "," << currGrid.y() << "," << currGrid.z() << ") expendable! " << endl;
    expand(currGrid, frDirEnum::E, dstMazeIdx1, dstMazeIdx2, centerPt);
    // gpu.test_expand(currGrid, frDirEnum::E, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no E" <<endl;
  // }
  // S
  // dir = frDirEnum::S;
  // cpures = isExpandable(currGrid, dir);
  // gpures = gpu.isEx(gridX, gridY, gridZ, dir);
  // if (cpures == gpures) {
  //   fmt::print("GPU result in ({}, {}, {}) FAILED!\n", gridX, gridY, gridZ);
  // }
  if (isExpandable(currGrid, frDirEnum::S)) {
    //cout << "(" << currGrid.x() << "," << currGrid.y() - 1 << "," << currGrid.z() << ") expendable! " << endl;
    expand(currGrid, frDirEnum::S, dstMazeIdx1, dstMazeIdx2, centerPt);
    // gpu.test_expand(currGrid, frDirEnum::S, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no S" <<endl;
  // }
  // W
  // dir = frDirEnum::W;
  // cpures = isExpandable(currGrid, dir);
  // gpures = gpu.isEx(gridX, gridY, gridZ, dir);
  // if (cpures == gpures) {
  //   fmt::print("GPU result in ({}, {}, {}) FAILED!\n", gridX, gridY, gridZ);
  // }
  if (isExpandable(currGrid, frDirEnum::W)) {
    //cout << "(" << currGrid.x() - 1 << "," << currGrid.y() << "," << currGrid.z() << ") expendable! " << endl;
    expand(currGrid, frDirEnum::W, dstMazeIdx1, dstMazeIdx2, centerPt);
    // cuExpand(currGrid, frDirEnum::W, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no W" <<endl;
  // }
  // U
  // dir = frDirEnum::U;
  // cpures = isExpandable(currGrid, dir);
  // gpures = gpu.isEx(gridX, gridY, gridZ, dir);
  // if (cpures == gpures) {
  //   fmt::print("GPU result in ({}, {}, {}) FAILED!\n", gridX, gridY, gridZ);
  // }
  if (isExpandable(currGrid, frDirEnum::U)) {
    //cout << "(" << currGrid.x() << "," << currGrid.y() << "," << currGrid.z() + 1 << ") expendable! " << endl;
    expand(currGrid, frDirEnum::U, dstMazeIdx1, dstMazeIdx2, centerPt);
    // gpu.test_expand(currGrid, frDirEnum::U, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no U" <<endl;
  // }
  // D
  // dir = frDirEnum::D;
  // cpures = isExpandable(currGrid, dir);
  // gpures = gpu.isEx(gridX, gridY, gridZ, dir);
  // if (cpures == gpures) {
  //   fmt::print("GPU result in ({}, {}, {}) FAILED!\n", gridX, gridY, gridZ);
  // }
  if (isExpandable(currGrid, frDirEnum::D)) {
    //cout << "(" << currGrid.x() << "," << currGrid.y() << "," << currGrid.z() - 1 << ") expendable! " << endl;
    expand(currGrid, frDirEnum::D, dstMazeIdx1, dstMazeIdx2, centerPt);
    // gpu.test_expand(currGrid, frDirEnum::D, dstMazeIdx1, dstMazeIdx2, centerPt);
  }
  // else {
  //   std::cout <<"no D" <<endl;
  // }
}

frCost FlexGridGraph::getEstCost(const FlexMazeIdx &src, const FlexMazeIdx &dstMazeIdx1,
                                const FlexMazeIdx &dstMazeIdx2, const frDirEnum &dir) const {
  //bool enableOutput = true;
  bool enableOutput = false;
  if (enableOutput) {
    cout <<"est from (" <<src.x() <<", " <<src.y() <<", " <<src.z() <<") "
         <<"to ("       <<dstMazeIdx1.x() <<", " <<dstMazeIdx1.y() <<", " <<dstMazeIdx1.z() <<") ("
                        <<dstMazeIdx2.x() <<", " <<dstMazeIdx2.y() <<", " <<dstMazeIdx2.z() <<")";
  }
  // bend cost
  int bendCnt = 0;
  int forbiddenPenalty = 0;
  frPoint srcPoint, dstPoint1, dstPoint2;
  getPoint(srcPoint, src.x(), src.y());
  getPoint(dstPoint1, dstMazeIdx1.x(), dstMazeIdx1.y());
  getPoint(dstPoint2, dstMazeIdx2.x(), dstMazeIdx2.y());
  frCoord minCostX = max(max(dstPoint1.x() - srcPoint.x(), srcPoint.x() - dstPoint2.x()), 0) * 1;
  frCoord minCostY = max(max(dstPoint1.y() - srcPoint.y(), srcPoint.y() - dstPoint2.y()), 0) * 1;
  frCoord minCostZ = max(max(getZHeight(dstMazeIdx1.z()) - getZHeight(src.z()), 
                       getZHeight(src.z()) - getZHeight(dstMazeIdx2.z())), 0) * 1;
  if (enableOutput) {
    cout <<" x/y/z min cost = (" <<minCostX <<", " <<minCostY <<", " <<minCostZ <<") " <<endl;
  }

  bendCnt += (minCostX && dir != frDirEnum::UNKNOWN && dir != frDirEnum::E && dir != frDirEnum::W) ? 1 : 0;
  bendCnt += (minCostY && dir != frDirEnum::UNKNOWN && dir != frDirEnum::S && dir != frDirEnum::N) ? 1 : 0;
  bendCnt += (minCostZ && dir != frDirEnum::UNKNOWN && dir != frDirEnum::U && dir != frDirEnum::D) ? 1 : 0;
  if (enableOutput) {
    cout << "  est cost = " << minCostX + minCostY + minCostZ + bendCnt << endl;
  }

  int gridX = src.x();
  int gridY = src.y();
  int gridZ = src.z();
  getNextGrid(gridX, gridY, gridZ, dir);
  frPoint nextPoint;
  getPoint(nextPoint, gridX, gridY);

  // avoid propagating to location that will cause fobidden via spacing to boundary pin
  if (DBPROCESSNODE == "GF14_13M_3Mx_2Cx_4Kx_2Hx_2Gx_LB") {
    if (drWorker && drWorker->getDRIter() >= 30 && drWorker->getRipupMode() == 0) {
      if (dstMazeIdx1 == dstMazeIdx2 && gridZ == dstMazeIdx1.z()) {
        auto layerNum = (gridZ + 1) * 2;
        auto layer = getDesign()->getTech()->getLayer(layerNum);
        bool isH = (layer->getDir() == frPrefRoutingDirEnum::frcHorzPrefRoutingDir);
        if (isH) {
          auto gap = abs(nextPoint.y() - dstPoint1.y());
          if (gap &&
              (getDesign()->getTech()->isVia2ViaForbiddenLen(gridZ, false, false, false, gap, false) || layerNum - 2 < BOTTOM_ROUTING_LAYER) &&
              (getDesign()->getTech()->isVia2ViaForbiddenLen(gridZ, true, true, false, gap, false) || layerNum + 2 > getDesign()->getTech()->getTopLayerNum())) {
            forbiddenPenalty = layer->getPitch() * ggDRCCost * 20;
          }
        } else {
          auto gap = abs(nextPoint.x() - dstPoint1.x());
          if (gap &&
              (getDesign()->getTech()->isVia2ViaForbiddenLen(gridZ, false, false, true, gap, false) || layerNum - 2 < BOTTOM_ROUTING_LAYER) &&
              (getDesign()->getTech()->isVia2ViaForbiddenLen(gridZ, true, true, true, gap, false) || layerNum + 2 > getDesign()->getTech()->getTopLayerNum())) {
            forbiddenPenalty = layer->getPitch() * ggDRCCost * 20;
          }
        }
      }
    }
  }

  return (minCostX + minCostY + minCostZ + bendCnt + forbiddenPenalty);

}

frDirEnum FlexGridGraph::getLastDir(const std::bitset<WAVEFRONTBITSIZE> &buffer) const {
  auto currDirVal = buffer.to_ulong() & 0b111u;
  return static_cast<frDirEnum>(currDirVal);
}

void FlexGridGraph::getNextGrid(frMIdx &gridX, frMIdx &gridY, frMIdx &gridZ, const frDirEnum dir) const {
  switch(dir) {
    case frDirEnum::E:
      ++gridX;
      break;
    case frDirEnum::S:
      --gridY;
      break;
    case frDirEnum::W:
      --gridX;
      break;
    case frDirEnum::N:
      ++gridY;
      break;
    case frDirEnum::U:
      ++gridZ;
      break;
    case frDirEnum::D:
      --gridZ;
      break;
    default:
      ;
  }
  return;
}

void FlexGridGraph::getPrevGrid(frMIdx &gridX, frMIdx &gridY, frMIdx &gridZ, const frDirEnum dir) const {
  switch(dir) {
    case frDirEnum::E:
      --gridX;
      break;
    case frDirEnum::S:
      ++gridY;
      break;
    case frDirEnum::W:
      ++gridX;
      break;
    case frDirEnum::N:
      --gridY;
      break;
    case frDirEnum::U:
      --gridZ;
      break;
    case frDirEnum::D:
      ++gridZ;
      break;
    default:
      ;
  }
  return;
}

/*inline*/ frCost FlexGridGraph::getNextPathCost(const FlexWavefrontGrid &currGrid, const frDirEnum &dir, bool test) const {
  // bool enableOutput = true;
  bool enableOutput = false;
  frMIdx gridX = currGrid.x();
  frMIdx gridY = currGrid.y();
  frMIdx gridZ = currGrid.z();
  frCost nextPathCost = currGrid.getPathCost();
  // bending cost
  auto currDir = currGrid.getLastDir();
  auto lNum = getLayerNum(currGrid.z());
  auto pathWidth = getDesign()->getTech()->getLayer(lNum)->getWidth();
  // printf("CPU currDir: %d\n", currDir);

  if (currDir != dir && currDir != frDirEnum::UNKNOWN) {
    // original
    ++nextPathCost;
    // printf("CPU Old cost add 1\n");
  }
  auto oldcost = nextPathCost;

  // via2viaForbiddenLen enablement
  if (dir == frDirEnum::U || dir == frDirEnum::D) {
    frCoord currVLengthX = 0;
    frCoord currVLengthY = 0;
    currGrid.getVLength(currVLengthX, currVLengthY);
    bool isCurrViaUp = (dir == frDirEnum::U);
    bool isForbiddenVia2Via = false;
    // check only y
    if (currVLengthX == 0 && currVLengthY > 0 && getTech()->isVia2ViaForbiddenLen(gridZ, !(currGrid.isPrevViaUp()), !isCurrViaUp, false, currVLengthY, false)) {
      isForbiddenVia2Via = true;
    // check only x
    } else if (currVLengthX > 0 && currVLengthY == 0 && getTech()->isVia2ViaForbiddenLen(gridZ, !(currGrid.isPrevViaUp()), !isCurrViaUp, true, currVLengthX, false)) {
      isForbiddenVia2Via = true;
    // check both x and y
    } else if (currVLengthX > 0 && currVLengthY > 0 && 
               (getTech()->isVia2ViaForbiddenLen(gridZ, !(currGrid.isPrevViaUp()), !isCurrViaUp, false, currVLengthY) &&
                getTech()->isVia2ViaForbiddenLen(gridZ, !(currGrid.isPrevViaUp()), !isCurrViaUp, true, currVLengthX))) {
      isForbiddenVia2Via = true;
    }

    if (isForbiddenVia2Via) {
      if (drWorker && drWorker->getDRIter() >= 3) {
        nextPathCost += ggMarkerCost * getEdgeLength(gridX, gridY, gridZ, dir);
      } else {
        nextPathCost += ggDRCCost * getEdgeLength(gridX, gridY, gridZ, dir);
      }
    }
  }

  // via2turn forbidden len enablement
  frCoord tLength    = std::numeric_limits<frCoord>::max();
  frCoord tLengthDummy = 0;
  bool    isTLengthViaUp = false;
  bool    isForbiddenTLen = false;
  if (currDir != frDirEnum::UNKNOWN && currDir != dir) {
    // next dir is a via
    if (dir == frDirEnum::U || dir == frDirEnum::D) {
      isTLengthViaUp = (dir == frDirEnum::U);
      // if there was a turn before
      if (tLength != std::numeric_limits<frCoord>::max()) {
        if (currDir == frDirEnum::W || currDir == frDirEnum::E) {
          tLength = currGrid.getTLength();
          if (getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, true, tLength)) {
            isForbiddenTLen = true;
          }
        } else if (currDir == frDirEnum::S || currDir == frDirEnum::N) {
          tLength = currGrid.getTLength();
          if (getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, false, tLength)) {
            isForbiddenTLen = true;
          }
        }
      }
    // curr is a planar turn
    } else {
      isTLengthViaUp = currGrid.isPrevViaUp();
      if (currDir == frDirEnum::W || currDir == frDirEnum::E) {
        currGrid.getVLength(tLength, tLengthDummy);
        /*
        auto vvv = getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, true, tLength);
        printf("CPU isViaForbiddenTurnLen: %d\n", vvv);
        */
        if (getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, true, tLength)) {
          isForbiddenTLen = true;
        }
      } else if (currDir == frDirEnum::S || currDir == frDirEnum::N) {
        currGrid.getVLength(tLengthDummy, tLength);
        /*
          auto vvv = getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, false, tLength);
          printf("CPU isViaForbiddenTurnLen: %d\n", vvv);
          */
        if (getTech()->isViaForbiddenTurnLen(gridZ, !isTLengthViaUp, false, tLength)) {
          isForbiddenTLen = true;
        }
      }
    }
    if (isForbiddenTLen) {
      if (drWorker && drWorker->getDRIter() >= 3) {
        nextPathCost += ggDRCCost * getEdgeLength(gridX, gridY, gridZ, dir);
      } else {
        nextPathCost += ggMarkerCost * getEdgeLength(gridX, gridY, gridZ, dir);
      }
    }
  }

  bool gridCost   = hasGridCost(gridX, gridY, gridZ, dir);
  bool drcCost    = hasDRCCost(gridX, gridY, gridZ, dir);
  bool markerCost = hasMarkerCost(gridX, gridY, gridZ, dir);
  bool shapeCost  = hasShapeCost(gridX, gridY, gridZ, dir);
  bool blockCost  = isBlocked(gridX, gridY, gridZ, dir);
  bool guideCost  = hasGuide(gridX, gridY, gridZ, dir);

  auto gridCostv = gridCost    ? GRIDCOST         * getEdgeLength(gridX, gridY, gridZ, dir) : 0;
  auto drcCostv = drcCost     ? ggDRCCost        * getEdgeLength(gridX, gridY, gridZ, dir) : 0;
  auto markerCostv = markerCost  ? ggMarkerCost     * getEdgeLength(gridX, gridY, gridZ, dir) : 0;
  auto shapeCostv = shapeCost   ? SHAPECOST        * getEdgeLength(gridX, gridY, gridZ, dir) : 0;
  auto blockCostv = blockCost   ? BLOCKCOST        * pathWidth * 20                          : 0;
  auto guideCostv =!guideCost ? GUIDECOST        * getEdgeLength(gridX, gridY, gridZ, dir) : 0;
  // temporarily disable guideCost
  nextPathCost += getEdgeLength(gridX, gridY, gridZ, dir)
                  + (gridCost   ? GRIDCOST         * getEdgeLength(gridX, gridY, gridZ, dir) : 0)
                  + (drcCost    ? ggDRCCost        * getEdgeLength(gridX, gridY, gridZ, dir) : 0)
                  + (markerCost ? ggMarkerCost     * getEdgeLength(gridX, gridY, gridZ, dir) : 0)
                  + (shapeCost  ? SHAPECOST        * getEdgeLength(gridX, gridY, gridZ, dir) : 0)
                  + (blockCost  ? BLOCKCOST        * pathWidth * 20                          : 0)
                  + (!guideCost ? GUIDECOST        * getEdgeLength(gridX, gridY, gridZ, dir) : 0);
  if (enableOutput) {
    cout <<"edge grid/shape/drc/marker/blk/length = " 
         <<hasGridCost(gridX, gridY, gridZ, dir)   <<"/"
         <<hasShapeCost(gridX, gridY, gridZ, dir)  <<"/"
         <<hasDRCCost(gridX, gridY, gridZ, dir)    <<"/"
         <<hasMarkerCost(gridX, gridY, gridZ, dir) <<"/"
         <<isBlocked(gridX, gridY, gridZ, dir) <<"/"
         <<getEdgeLength(gridX, gridY, gridZ, dir) <<endl;
  }
  /*
  printf("CPU Costs %d: \n%d %d %d %d %d %d\n", nextPathCost, 
      getEdgeLength(gridX, gridY, gridZ, dir), gridCost, drcCost, markerCost, 
      shapeCost, blockCost, guideCost);
      */
  // if (test) {
  // printf("CPU Costs %u: old cost: %u\n", nextPathCost, oldcost);
  // }
  return nextPathCost;

}

/*inline*/ FlexMazeIdx FlexGridGraph::getTailIdx(const FlexMazeIdx &currIdx, const FlexWavefrontGrid &currGrid) const {
  int gridX = currIdx.x();
  int gridY = currIdx.y();
  int gridZ = currIdx.z();
  auto backTraceBuffer = currGrid.getBackTraceBuffer();
  for (int i = 0; i < WAVEFRONTBUFFERSIZE; ++i) {
    int currDirVal = backTraceBuffer.to_ulong() - ((backTraceBuffer.to_ulong() >> DIRBITSIZE) << DIRBITSIZE);
    frDirEnum currDir = static_cast<frDirEnum>(currDirVal);
    backTraceBuffer >>= DIRBITSIZE;
    getPrevGrid(gridX, gridY, gridZ, currDir);
  }
  return FlexMazeIdx(gridX, gridY, gridZ);
}

/*inline*/ bool FlexGridGraph::isExpandable(const FlexWavefrontGrid &currGrid, frDirEnum dir) const {
  //bool enableOutput = true;
  bool enableOutput = false;
  frMIdx gridX = currGrid.x();
  frMIdx gridY = currGrid.y();
  frMIdx gridZ = currGrid.z();
  bool hg = hasEdge(gridX, gridY, gridZ, dir);
  // fmt::print("CPU hasEdge: {}\n", hg);
  if (enableOutput) {
    if (!hasEdge(gridX, gridY, gridZ, dir)) {
      cout <<"no edge@(" <<gridX <<", " <<gridY <<", " <<gridZ <<") " <<(int)dir <<endl;
    }
    if (hasEdge(gridX, gridY, gridZ, dir) && !hasGuide(gridX, gridY, gridZ, dir)) {
      cout <<"no guide@(" <<gridX <<", " <<gridY <<", " <<gridZ <<") " <<(int)dir <<endl;
    }
  }
  reverse(gridX, gridY, gridZ, dir);
  if (!hg || 
      isSrc(gridX, gridY, gridZ) || 
      (getPrevAstarNodeDir(gridX, gridY, gridZ) != frDirEnum::UNKNOWN) || // comment out for non-buffer enablement
      currGrid.getLastDir() == dir) {
    return false;
  } else {
    return true;
  }
}

void FlexGridGraph::traceBackPath(const FlexWavefrontGrid &currGrid, vector<FlexMazeIdx> &path, vector<FlexMazeIdx> &root,
                                  FlexMazeIdx &ccMazeIdx1, FlexMazeIdx &ccMazeIdx2) const {
  bool enableOutput = true;
  // bool enableOutput = false;
  if (enableOutput) {
    cout << "    start traceBackPath...\n";
  }
  frDirEnum prevDir = frDirEnum::UNKNOWN, currDir = frDirEnum::UNKNOWN;
  int currX = currGrid.x(), currY = currGrid.y(), currZ = currGrid.z();
  // pop content in buffer
  auto backTraceBuffer = currGrid.getBackTraceBuffer();
  for (int i = 0; i < WAVEFRONTBUFFERSIZE; ++i) {
    // current grid is src
    if (isSrc(currX, currY, currZ)) {
      break;
    }
    // get last direction
    currDir = getLastDir(backTraceBuffer);
    backTraceBuffer >>= DIRBITSIZE;
    if (currDir == frDirEnum::UNKNOWN) {
      cout << "Warning: unexpected direction in tracBackPath\n";
      break;
    }
    root.push_back(FlexMazeIdx(currX, currY, currZ));
    // push point to path
    if (currDir != prevDir) {
      path.push_back(FlexMazeIdx(currX, currY, currZ));
      if (enableOutput) {
        cout <<" -- (" <<currX <<", " <<currY <<", " <<currZ <<")";
      }
    }
    getPrevGrid(currX, currY, currZ, currDir);
    prevDir = currDir;
  }
  // trace back according to grid prev dir
  while (isSrc(currX, currY, currZ) == false) {
    // get last direction
    currDir = getPrevAstarNodeDir(currX, currY, currZ);
    root.push_back(FlexMazeIdx(currX, currY, currZ));
    if (currDir == frDirEnum::UNKNOWN) {
      cout << "Warning: unexpected direction in tracBackPath\n";
      break;
    }
    if (currDir != prevDir) {
      path.push_back(FlexMazeIdx(currX, currY, currZ));
      if (enableOutput) {
        cout <<" -- (" <<currX <<", " <<currY <<", " <<currZ <<")";
      }
    }
    getPrevGrid(currX, currY, currZ, currDir);
    prevDir = currDir;
  }
  // add final path to src, only add when path exists; no path exists (src = dst)
  if (!path.empty()) {
    path.push_back(FlexMazeIdx(currX, currY, currZ));
    if (enableOutput) {
      cout <<" -- (" <<currX <<", " <<currY <<", " <<currZ <<")";
    }
  }
  for (auto &mi: path) {
    ccMazeIdx1.set(min(ccMazeIdx1.x(), mi.x()),
                   min(ccMazeIdx1.y(), mi.y()),
                   min(ccMazeIdx1.z(), mi.z()));
    ccMazeIdx2.set(max(ccMazeIdx2.x(), mi.x()),
                   max(ccMazeIdx2.y(), mi.y()),
                   max(ccMazeIdx2.z(), mi.z()));
  }
  if (enableOutput) {
    cout <<endl;
  }

}

bool FlexGridGraph::search(vector<FlexMazeIdx> &connComps, drPin* nextPin, vector<FlexMazeIdx> &path, 
                           FlexMazeIdx &ccMazeIdx1, FlexMazeIdx &ccMazeIdx2, const frPoint &centerPt) {
  bool enableOutput = true;
  cout << "ccMazeIdx1: (" << ccMazeIdx1.x() << "," << ccMazeIdx1.y() << "," << ccMazeIdx1.z() << ")" << endl;
  cout << "ccMazeIdx2: (" << ccMazeIdx2.x() << "," << ccMazeIdx2.y() << "," << ccMazeIdx2.z() << ")" << endl;
  cout << "centerPt: (" << centerPt.x()  << "," << centerPt.y() << ")" << endl;
  //bool enableOutput = false;
  int stepCnt = 0;

  // prep nextPinBox
  frMIdx xDim, yDim, zDim;
  getDim(xDim, yDim, zDim);


  FlexMazeIdx dstMazeIdx1(xDim - 1, yDim - 1, zDim - 1);
  FlexMazeIdx dstMazeIdx2(0, 0, 0);
  FlexMazeIdx mi;
  for (auto &ap: nextPin->getAccessPatterns()) {
    ap->getMazeIdx(mi);
    dstMazeIdx1.set(min(dstMazeIdx1.x(), mi.x()),
                    min(dstMazeIdx1.y(), mi.y()),
                    min(dstMazeIdx1.z(), mi.z()));
    dstMazeIdx2.set(max(dstMazeIdx2.x(), mi.x()),
                    max(dstMazeIdx2.y(), mi.y()),
                    max(dstMazeIdx2.z(), mi.z()));
  }
  cout << "dstMazeIdx1: (" << dstMazeIdx1.x() << "," << dstMazeIdx1.y() << "," << dstMazeIdx1.z() << ")" << endl;
  cout << "dstMazeIdx2: (" << dstMazeIdx2.x() << "," << dstMazeIdx2.y() << "," << dstMazeIdx2.z() << ")" << endl;

  wavefront.cleanup();
  // init wavefront
  frPoint currPt;
  cout << "cnncomps size: " << connComps.size() << endl;
  for (auto &idx: connComps) {
    if (isDst(idx.x(), idx.y(), idx.z())) {
      if (enableOutput) {
        cout <<"message: astarSearch dst covered (" <<idx.x() <<", " <<idx.y() <<", " <<idx.z() <<")" <<endl;
      }
      path.push_back(FlexMazeIdx(idx.x(), idx.y(), idx.z()));
      return true;
    }
    // get min area min length
    auto lNum = getLayerNum(idx.z());
    auto minAreaConstraint = getDesign()->getTech()->getLayer(lNum)->getAreaConstraint();
    frCoord fakeArea = minAreaConstraint ? minAreaConstraint->getMinArea() : 0;
    getPoint(currPt, idx.x(), idx.y());
    frCoord currDist = abs(currPt.x() - centerPt.x()) + abs(currPt.y() - centerPt.y());
    FlexWavefrontGrid currGrid(idx.x(), idx.y(), idx.z(), fakeArea, 
                               std::numeric_limits<frCoord>::max(), std::numeric_limits<frCoord>::max(), true, 
                               std::numeric_limits<frCoord>::max(),
                               currDist, 0, getEstCost(idx, dstMazeIdx1, dstMazeIdx2, frDirEnum::UNKNOWN));
    wavefront.push(currGrid);
    if (enableOutput) {
      cout <<"src add to wavefront (" <<idx.x() <<", " <<idx.y() <<", " <<idx.z() <<")" <<endl;
    }
  }
  while(!wavefront.empty()) {
    auto currGrid = wavefront.top();
    wavefront.pop();
    if (getPrevAstarNodeDir(currGrid.x(), currGrid.y(), currGrid.z()) != frDirEnum::UNKNOWN) {
      continue;
    }

    // test
    if (enableOutput) {
      ++stepCnt;
      // cout << " (" <<currGrid.x() << "," << currGrid.y() << "," << currGrid.z() << ") -> ";
    }
    // if (stepCnt % 100000 == 0) {
    //   std::cout << "wavefront size = " << wavefront.size() << " at step = " << stepCnt << "\n";
    // }
    if (isDst(currGrid.x(), currGrid.y(), currGrid.z())) {
      traceBackPath(currGrid, path, connComps, ccMazeIdx1, ccMazeIdx2);
      cout << "ccMazeIdx1: (" << ccMazeIdx1.x() << "," << ccMazeIdx1.y() << "," << ccMazeIdx1.z() << ")" << endl;
      cout << "ccMazeIdx2: (" << ccMazeIdx2.x() << "," << ccMazeIdx2.y() << "," << ccMazeIdx2.z() << ")" << endl;
      if (enableOutput) {
        cout << "path found. stepCnt = " << stepCnt << "\n";
        cout << "CUDA free space: " << xDim << " * " << yDim << " * " << zDim << endl;
      }
      return true;
    } else {
      // expand and update wavefront
      expandWavefront(currGrid, dstMazeIdx1, dstMazeIdx2, centerPt);
    }
    
  }
  return false;
}

bool FlexGridGraph::cuSearch(vector<FlexMazeIdx> &connComps, drPin* nextPin, vector<FlexMazeIdx> &path, 
                           FlexMazeIdx &ccMazeIdx1, FlexMazeIdx &ccMazeIdx2, const frPoint &centerPt) {
  bool enableOutput = true;
  cout << "ccMazeIdx1: (" << ccMazeIdx1.x() << "," << ccMazeIdx1.y() << "," << ccMazeIdx1.z() << ")" << endl;
  cout << "ccMazeIdx2: (" << ccMazeIdx2.x() << "," << ccMazeIdx2.y() << "," << ccMazeIdx2.z() << ")" << endl;
  cout << "centerPt: (" << centerPt.x()  << "," << centerPt.y() << ")" << endl;
  //bool enableOutput = false;
  int stepCnt = 0;

  // prep nextPinBox
  frMIdx xDim, yDim, zDim;
  getDim(xDim, yDim, zDim);



  FlexMazeIdx dstMazeIdx1(xDim - 1, yDim - 1, zDim - 1);
  FlexMazeIdx dstMazeIdx2(0, 0, 0);
  FlexMazeIdx mi;
  for (auto &ap: nextPin->getAccessPatterns()) {
    ap->getMazeIdx(mi);
    dstMazeIdx1.set(min(dstMazeIdx1.x(), mi.x()),
                    min(dstMazeIdx1.y(), mi.y()),
                    min(dstMazeIdx1.z(), mi.z()));
    dstMazeIdx2.set(max(dstMazeIdx2.x(), mi.x()),
                    max(dstMazeIdx2.y(), mi.y()),
                    max(dstMazeIdx2.z(), mi.z()));
  }
  cout << "dstMazeIdx1: (" << dstMazeIdx1.x() << "," << dstMazeIdx1.y() << "," << dstMazeIdx1.z() << ")" << endl;
  cout << "dstMazeIdx2: (" << dstMazeIdx2.x() << "," << dstMazeIdx2.y() << "," << dstMazeIdx2.z() << ")" << endl;

  wavefront.cleanup();
  // init wavefront
  frPoint currPt;
  cout << "cnncomps size: " << connComps.size() << endl;
  for (auto &idx: connComps) {
    if (isDst(idx.x(), idx.y(), idx.z())) {
      if (enableOutput) {
        cout <<"message: astarSearch dst covered (" <<idx.x() <<", " <<idx.y() <<", " <<idx.z() <<")" <<endl;
      }
      path.push_back(FlexMazeIdx(idx.x(), idx.y(), idx.z()));
      return true;
    }
    // get min area min length
    auto lNum = getLayerNum(idx.z());
    auto minAreaConstraint = getDesign()->getTech()->getLayer(lNum)->getAreaConstraint();
    frCoord fakeArea = minAreaConstraint ? minAreaConstraint->getMinArea() : 0;
    getPoint(currPt, idx.x(), idx.y());
    frCoord currDist = abs(currPt.x() - centerPt.x()) + abs(currPt.y() - centerPt.y());
    FlexWavefrontGrid currGrid(idx.x(), idx.y(), idx.z(), fakeArea, 
                               std::numeric_limits<frCoord>::max(), std::numeric_limits<frCoord>::max(), true, 
                               std::numeric_limits<frCoord>::max(),
                               currDist, 0, getEstCost(idx, dstMazeIdx1, dstMazeIdx2, frDirEnum::UNKNOWN));
    wavefront.push(currGrid);
    if (enableOutput) {
      cout <<"src add to wavefront (" <<idx.x() <<", " <<idx.y() <<", " <<idx.z() <<")" <<endl;
    }
  }
    bool goon = true;
  while(!wavefront.empty()) {
    auto currGrid = wavefront.top();
    wavefront.pop();
    if (getPrevAstarNodeDir(currGrid.x(), currGrid.y(), currGrid.z()) != frDirEnum::UNKNOWN) {
      continue;
    }

    /*
    if (goon) {
      auto gpu = GPUPathwaySolver();
      gpu.initialize(bits, prevDirs, srcs, guides, zDirs, 
          xCoords.size(), yCoords.size(), zCoords.size());
      auto x = 2;
      auto z = 2;
      auto y = 2;
      auto dir = frDirEnum::N;
      auto gpures = gpu.isEx(2, 2, 2, dir);
      auto cpures = getPrevAstarNodeDir(2, 2, 2);
      fmt::print("GPU prevD result: {}\n", gpures);
      fmt::print("CPU prevD result: {}\n",  cpures);
      goon = false;
    }
    */
    // test
    if (enableOutput) {
      ++stepCnt;
      // cout << " (" <<currGrid.x() << "," << currGrid.y() << "," << currGrid.z() << ") -> ";
    }
    // if (stepCnt % 100000 == 0) {
    //   std::cout << "wavefront size = " << wavefront.size() << " at step = " << stepCnt << "\n";
    // }
    if (isDst(currGrid.x(), currGrid.y(), currGrid.z())) {
      traceBackPath(currGrid, path, connComps, ccMazeIdx1, ccMazeIdx2);
      cout << "ccMazeIdx1: (" << ccMazeIdx1.x() << "," << ccMazeIdx1.y() << "," << ccMazeIdx1.z() << ")" << endl;
      cout << "ccMazeIdx2: (" << ccMazeIdx2.x() << "," << ccMazeIdx2.y() << "," << ccMazeIdx2.z() << ")" << endl;
      if (enableOutput) {
        cout << "path found. stepCnt = " << stepCnt << "\n";
      }
      fmt::print("=========isEx GPU========\n");
      fmt::print("GPU Passed: {}, Failed: {}\n", gpu_pass, gpu_failed);
      fmt::print("===============================\n");
      return true;
    } else {
      // expand and update wavefront
      expandWavefront(currGrid, dstMazeIdx1, dstMazeIdx2, centerPt);
    }
    
  }
  return false;
}

void FlexGridGraph::viaData(vector<int> &dest, const vector<vector<vector<std::pair<frCoord, frCoord>>>> &data) {
  int total_size = 0;
  int sizex = data.size();
  int sizey = data[0].size();
  // fmt::print("sizex: {}, sizey: {}\n", sizex, sizey);
  for (auto &two_d_vec: data){
    total_size += two_d_vec.size();
  }
  if (sizex * sizey != total_size) {
    // fmt::print("Total size: {}, sizex: {}, sizey: {}\n", total_size, sizex, sizey);
  }

  auto flatvec = dest;
  for (int i = 0; i < sizex; ++i) {
    for (int j = 0; j < sizey; ++j) {
      auto &p = data[i][j];
      if (p.size() > 1) {
      // fmt::print("size of data[{}][{}]: {}\n", i, j, p.size());
      }
      if (p.size() > 0) {
        flatvec.push_back(p[0].first);
        flatvec.push_back(p[0].second);
      } else {
        flatvec.push_back(-114514);
        flatvec.push_back(-114514);
      }
    }
  }
}
