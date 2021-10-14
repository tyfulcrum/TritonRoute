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

#ifndef _FLEX_MAZE_TYPES_H_
#define _FLEX_MAZE_TYPES_H_

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 
#include <frCoreLangTypes.h>
using coret::frMIdx;

namespace fr {
  class FlexMazeIdx {
  public:
    CUDA_CALLABLE_MEMBER FlexMazeIdx(): xIdx(-1), yIdx(-1), zIdx(-1) {}
    CUDA_CALLABLE_MEMBER FlexMazeIdx(frMIdx xIn, frMIdx yIn, frMIdx zIn): xIdx(xIn), yIdx(yIn), zIdx(zIn) {}
    // getters
    CUDA_CALLABLE_MEMBER frMIdx x() const {
      return xIdx;
    }
    CUDA_CALLABLE_MEMBER frMIdx y() const {
      return yIdx;
    }
    CUDA_CALLABLE_MEMBER frMIdx z() const {
      return zIdx;
    }
    CUDA_CALLABLE_MEMBER bool empty() const {
      return (xIdx == -1 && yIdx == -1 && zIdx == -1);
    }
    // setters
    CUDA_CALLABLE_MEMBER void set(frMIdx xIn, frMIdx yIn, frMIdx zIn) {
      xIdx = xIn;
      yIdx = yIn;
      zIdx = zIn;
    }
    CUDA_CALLABLE_MEMBER void set(const FlexMazeIdx &in) {
      xIdx = in.x();
      yIdx = in.y();
      zIdx = in.z();
    }
    // others
    CUDA_CALLABLE_MEMBER bool operator<(const FlexMazeIdx &rhs) const {
      if (this->xIdx != rhs.x()) {
        return this->xIdx < rhs.x();
      } else if (this->yIdx != rhs.y()) {
        return this->yIdx < rhs.y();
      } else {
        return this->zIdx < rhs.z();
      }
    }
    CUDA_CALLABLE_MEMBER bool operator==(const FlexMazeIdx &rhs) const {
      return (xIdx == rhs.xIdx && yIdx == rhs.yIdx && zIdx == rhs.zIdx);
    }

    CUDA_CALLABLE_MEMBER friend std::ostream& operator<<(std::ostream& os, const FlexMazeIdx &mIdx) {
      os <<"(" <<mIdx.x() <<", " <<mIdx.y() <<", " <<mIdx.z() <<")";
      return os;
    }

  protected:
    frMIdx xIdx;
    frMIdx yIdx;
    frMIdx zIdx;
  };
}

#endif
