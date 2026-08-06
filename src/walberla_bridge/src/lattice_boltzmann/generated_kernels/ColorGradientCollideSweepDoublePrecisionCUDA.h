//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file ColorGradientCollideSweepDoublePrecisionCUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"

#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"

#include <functional>
#include <unordered_map>



#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {


class ColorGradientCollideSweepDoublePrecisionCUDA
{
public:
   ColorGradientCollideSweepDoublePrecisionCUDA( BlockDataID force_aID_, BlockDataID force_bID_, BlockDataID pdfs_aID_, BlockDataID pdfs_bID_, BlockDataID phasefieldID_, BlockDataID rho_aID_, BlockDataID rho_bID_, BlockDataID velocityID_, double beta, double kT, double omega_shear_a, double omega_shear_b, uint32_t seed, double sigma, uint32_t time_step )
     : force_aID(force_aID_), force_bID(force_bID_), pdfs_aID(pdfs_aID_), pdfs_bID(pdfs_bID_), phasefieldID(phasefieldID_), rho_aID(rho_aID_), rho_bID(rho_bID_), velocityID(velocityID_), beta_(beta), kT_(kT), omega_shear_a_(omega_shear_a), omega_shear_b_(omega_shear_b), seed_(seed), sigma_(sigma), time_step_(time_step), block_offset_0_(uint32_t(0)) , block_offset_1_(uint32_t(0)) , block_offset_2_(uint32_t(0)) , configured_(false)
   {}

   

   void run(IBlock * block, gpuStream_t stream = nullptr);

   void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream = nullptr);

   
   void operator() (IBlock * block, gpuStream_t stream = nullptr)
   {
     run(block, stream);
   }
   

   static std::function<void (IBlock *)> getSweep(const shared_ptr<ColorGradientCollideSweepDoublePrecisionCUDA> & kernel)
   {
     return [kernel]
            (IBlock * b)
            { kernel->run(b); };
   }

   static std::function<void (IBlock*, gpuStream_t )> getSweepOnCellInterval(const shared_ptr<ColorGradientCollideSweepDoublePrecisionCUDA> & kernel, const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
   {
     return [kernel, blocks, globalCellInterval, ghostLayers]
            (IBlock * b, gpuStream_t stream = nullptr)
            { kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b, stream); };
   }

   std::function<void (IBlock *)> getSweep(gpuStream_t stream = nullptr)
   {
     return [this, stream]
            (IBlock * b)
            { this->run(b, stream); };
   }

   std::function<void (IBlock *)> getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1, gpuStream_t stream = nullptr)
   {
     return [this, blocks, globalCellInterval, ghostLayers, stream]
            (IBlock * b)
            { this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b, stream); };
   }

   
   void configure( const shared_ptr<StructuredBlockStorage> & blocks, IBlock * block )
   {
   Cell BlockCellBB = blocks->getBlockCellBB( *block).min();
   block_offset_0_ = uint32_t(BlockCellBB[0]);
   block_offset_1_ = uint32_t(BlockCellBB[1]);
   block_offset_2_ = uint32_t(BlockCellBB[2]);
   configured_ = true;
   }
   

   

   inline double getBeta() const { return beta_; }
   inline uint32_t getBlock_offset_0() const { return block_offset_0_; }
   inline uint32_t getBlock_offset_1() const { return block_offset_1_; }
   inline uint32_t getBlock_offset_2() const { return block_offset_2_; }
   inline double getKt() const { return kT_; }
   inline double getOmega_shear_a() const { return omega_shear_a_; }
   inline double getOmega_shear_b() const { return omega_shear_b_; }
   inline uint32_t getSeed() const { return seed_; }
   inline double getSigma() const { return sigma_; }
   inline uint32_t getTime_step() const { return time_step_; }
   inline void setBeta(const double value) { beta_ = value; }
   inline void setBlock_offset_0(const uint32_t value) { block_offset_0_ = value; }
   inline void setBlock_offset_1(const uint32_t value) { block_offset_1_ = value; }
   inline void setBlock_offset_2(const uint32_t value) { block_offset_2_ = value; }
   inline void setKt(const double value) { kT_ = value; }
   inline void setOmega_shear_a(const double value) { omega_shear_a_ = value; }
   inline void setOmega_shear_b(const double value) { omega_shear_b_ = value; }
   inline void setSeed(const uint32_t value) { seed_ = value; }
   inline void setSigma(const double value) { sigma_ = value; }
   inline void setTime_step(const uint32_t value) { time_step_ = value; }

private:
   
   BlockDataID force_aID;
   BlockDataID force_bID;
   BlockDataID pdfs_aID;
   BlockDataID pdfs_bID;
   BlockDataID phasefieldID;
   BlockDataID rho_aID;
   BlockDataID rho_bID;
   BlockDataID velocityID;
   double beta_;
   uint32_t block_offset_0_;
   uint32_t block_offset_1_;
   uint32_t block_offset_2_;
   double kT_;
   double omega_shear_a_;
   double omega_shear_b_;
   uint32_t seed_;
   double sigma_;
   uint32_t time_step_;

   

   bool configured_;
   
};


} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
