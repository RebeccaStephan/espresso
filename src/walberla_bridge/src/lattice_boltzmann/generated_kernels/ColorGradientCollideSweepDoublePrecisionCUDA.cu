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
//! \\file ColorGradientCollideSweepDoublePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepDoublePrecisionCUDA.h"




#define FUNC_PREFIX __global__

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wfloat-equal"
#   pragma GCC diagnostic ignored "-Wshadow"
#   pragma GCC diagnostic ignored "-Wconversion"
#   pragma GCC diagnostic ignored "-Wunused-variable"
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning push
#pragma warning( disable :  1599 )
#endif

using namespace std;

namespace walberla {
namespace pystencils {


namespace internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double beta, double omega_shear_a, double omega_shear_b, double sigma)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_force_a_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_force_a_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_force_a_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      const double xi_2 = ((1.0) / (omega_shear_a + omega_shear_b));
      const double xi_3 = omega_shear_a*omega_shear_b;
      const double xi_4 = xi_2*xi_3;
      const double xi_5 = xi_4*2.0;
      const double xi_6 = omega_shear_a*-4.0;
      const double xi_7 = xi_4*8.0;
      const double xi_8 = xi_6 + xi_7;
      const double xi_9 = -xi_8;
      const double xi_12 = omega_shear_b*-4.0;
      const double xi_13 = xi_12 + xi_7;
      const double xi_206 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_207 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_208 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_209 = -0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_210 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_211 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_212 = -xi_208 + xi_209 + xi_210 + xi_211;
      const double xi_213 = -0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_214 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_215 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_216 = 0.013888888888888888*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_217 = xi_213 - xi_214 + xi_215 + xi_216;
      const double xi_218 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_219 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_220 = -xi_218 + xi_219;
      const double xi_221 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_222 = -xi_221;
      const double xi_223 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_224 = xi_222 + xi_223;
      const double xi_225 = -0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_226 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_227 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_228 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_229 = -xi_227 + xi_228;
      const double xi_230 = xi_225 + xi_226 + xi_229;
      const double xi_231 = xi_224 + xi_230;
      const double xi_232 = -xi_206 + xi_207 + xi_212 + xi_217 + xi_220 + xi_231;
      const double xi_233 = xi_232*xi_232;
      const double xi_234 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_235 = -xi_234;
      const double xi_236 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_237 = xi_222 + xi_223 + xi_235 + xi_236;
      const double xi_238 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_239 = 0.055555555555555552*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_240 = -xi_238 + xi_239;
      const double xi_241 = xi_206 - xi_207;
      const double xi_242 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_243 = -xi_242;
      const double xi_244 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_245 = xi_227 - xi_228 + xi_243 + xi_244;
      const double xi_246 = xi_212 + xi_214 - xi_216 + xi_237 + xi_240 + xi_241 + xi_245;
      const double xi_247 = xi_246*xi_246;
      const double xi_248 = xi_238 - xi_239;
      const double xi_249 = xi_235 + xi_236;
      const double xi_250 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_251 = -xi_250;
      const double xi_252 = 0.22222222222222221*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_253 = xi_218 - xi_219 + xi_251 + xi_252;
      const double xi_254 = xi_208 + xi_209 + xi_210 - xi_211 + xi_217 + xi_241 + xi_248 + xi_249 + xi_253;
      const double xi_255 = xi_254*xi_254;
      const double xi_256 = xi_233 + xi_247 + xi_255;
      const double xi_257 = pow(xi_256, 0.5);
      const double xi_258 = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const double xi_259 = sigma*xi_257*((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.5) ? (omega_shear_a): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] < -0.5) ? (omega_shear_b): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.0) ? (xi_258*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]): (xi_13*xi_258 + xi_13*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] + xi_5))));
      const bool xi_260 = xi_257 > 0.0;
      const double xi_261 = ((xi_260) ? (xi_259*0.25): (0.0));
      const double xi_279 = ((1.0) / (xi_256));
      const double xi_280 = xi_259*1.125;
      const double xi_281 = ((xi_260) ? (xi_280*(xi_247*xi_279*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_295 = ((xi_260) ? (xi_280*(xi_279*0.055555555555555552*(xi_246*xi_246) - 0.018518518518518517)): (0.0));
      const double xi_302 = ((xi_260) ? (xi_280*(xi_279*0.055555555555555552*(xi_254*xi_254) - 0.018518518518518517)): (0.0));
      const double xi_305 = ((xi_260) ? (xi_280*(xi_255*xi_279*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_306 = -((xi_260) ? (xi_280*(xi_279*0.055555555555555552*(xi_232*xi_232) - 0.018518518518518517)): (0.0));
      const double xi_313 = ((xi_260) ? (xi_280*(xi_233*xi_279*0.055555555555555552 - 0.018518518518518517)): (0.0));
      const double xi_318 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_319 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_320 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_321 = xi_213 + xi_215;
      const double xi_322 = xi_253 + xi_321;
      const double xi_323 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_324 = xi_242 - xi_244 + xi_323;
      const double xi_325 = xi_221 - xi_223 + xi_229 - xi_318 + xi_319 + xi_320 + xi_322 + xi_324;
      const double xi_326 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
      const double xi_332 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_333 = 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const double xi_334 = -0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_335 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const double xi_336 = xi_224 + xi_245 + xi_322 - xi_332 + xi_333 + xi_334 + xi_335;
      const double xi_337 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_336*xi_336) - 0.037037037037037035)): (0.0));
      const double xi_338 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_336*xi_336) - 0.037037037037037035)): (0.0));
      const double xi_340 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_325*xi_325) - 0.037037037037037035)): (0.0));
      const double xi_341 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_342 = xi_220 + xi_225 + xi_226 + xi_321;
      const double xi_343 = xi_234 - xi_236 + xi_332 - xi_333;
      const double xi_344 = xi_248 + xi_324 + xi_341 + xi_342 + xi_343;
      const double xi_345 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_344*xi_344) - 0.037037037037037035)): (0.0));
      const double xi_359 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_360 = xi_240 + xi_318 - xi_319;
      const double xi_361 = xi_243 + xi_244 + xi_249 + xi_334 + xi_342 + xi_359 + xi_360;
      const double xi_362 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_361*xi_361) - 0.037037037037037035)): (0.0));
      const double xi_367 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_368 = xi_230 + xi_237 + xi_248 + xi_251 + xi_252 + xi_323 + xi_334 + xi_367;
      const double xi_369 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_368*xi_368) - 0.037037037037037035)): (0.0));
      const double xi_376 = -0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const double xi_377 = xi_231 + xi_250 - xi_252 + xi_343 + xi_360 + xi_376;
      const double xi_378 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_377*xi_377) - 0.037037037037037035)): (0.0));
      const double xi_380 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_361*xi_361) - 0.037037037037037035)): (0.0));
      const double xi_382 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_344*xi_344) - 0.037037037037037035)): (0.0));
      const double xi_383 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_377*xi_377) - 0.037037037037037035)): (0.0));
      const double xi_385 = ((xi_260) ? (xi_280*(xi_279*0.027777777777777776*(xi_368*xi_368) - 0.037037037037037035)): (0.0));
      const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
      const double xi_272 = -xia_3_collide;
      const double xia_4_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
      const double xia_5_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
      const double xia_6_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const double xia_7_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
      const double xi_50 = xia_7_collide*0.083333333333333329;
      const double xia_8_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_0 = xia_8_collide > 0.5;
      const bool xi_1 = xia_8_collide < -0.5;
      const double xi_10 = xia_8_collide*xia_8_collide;
      const bool xi_11 = xia_8_collide > 0.0;
      const double xia_9_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
      const double xi_41 = xia_9_collide*0.083333333333333329;
      const double xia_10_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
      const double xi_203 = xia_10_collide*2.0 + xia_4_collide*2.0;
      const double xi_296 = xia_10_collide - xia_4_collide;
      const double xia_11_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
      const double xi_289 = -xia_11_collide;
      const double xia_12_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const double xia_13_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
      const double xia_14_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const double xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
      const double xia_16_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
      const double xi_270 = -xia_16_collide;
      const double xia_17_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
      const double xi_194 = xia_16_collide + xia_17_collide;
      const double xi_263 = -xia_17_collide;
      const double xi_264 = xi_263 + xia_16_collide;
      const double xia_18_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
      const double xi_198 = xia_15_collide + xia_18_collide;
      const double xi_274 = xia_15_collide*2.0 + xia_18_collide*-2.0;
      const double xia_19_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
      const double xi_271 = xia_19_collide*2.0;
      const double xia_20_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
      const double xi_199 = xi_198 + xia_19_collide + xia_20_collide;
      const double xi_283 = -xia_19_collide - xia_20_collide;
      const double xia_21_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
      const double xi_195 = xia_21_collide + xia_3_collide;
      const double xi_262 = -xia_21_collide;
      const double xi_273 = xi_272 + xia_21_collide;
      const double xia_22_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const double xia_23_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
      const double xi_26 = xia_23_collide*0.083333333333333329;
      const double xi_60 = xia_23_collide*0.25;
      const double xia_24_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
      const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
      const double xi_201 = xia_25_collide*2.0;
      const double xia_26_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
      const double xi_189 = -xia_26_collide;
      const double xi_190 = xi_189 + xia_13_collide*5.0 + xia_16_collide*2.0 + xia_17_collide*2.0 + xia_21_collide*2.0 + xia_24_collide*5.0 + xia_3_collide*2.0;
      const double xia_27_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
      const double xi_288 = -xia_27_collide;
      const double xi_308 = xi_288 + xia_11_collide;
      const double xi_309 = xi_194 + xi_262 + xi_272 + xi_308;
      const double xia_28_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
      const double xi_265 = -xia_28_collide;
      const double xi_266 = xi_265 + xia_5_collide;
      const double xi_267 = -xi_262 - xi_264 - xi_266 - xia_3_collide;
      const double xi_275 = -xi_266 - xi_270 - xi_271 - xi_273 - xi_274 - xia_17_collide + xia_20_collide*2.0;
      const double xia_29_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
      const double xi_202 = xia_29_collide*2.0;
      const double xi_204 = xi_190 + xi_201 + xi_202 + xi_203 + xia_11_collide*-4.0 + xia_15_collide*-7.0 + xia_18_collide*-7.0 + xia_19_collide*-7.0 + xia_20_collide*-7.0 + xia_27_collide*-4.0 + xia_28_collide*5.0 + xia_5_collide*5.0;
      const double xi_297 = xi_296 + xia_25_collide - xia_29_collide;
      const double xi_298 = -xi_297 - xia_13_collide + xia_24_collide;
      const double xi_300 = -xi_271 + xi_274 + xi_297 - xia_13_collide + xia_20_collide*2.0 + xia_24_collide;
      const double xi_310 = xi_195 - xi_201 - xi_202 + xi_203 + xi_263 + xi_270 + xi_308;
      const double xia_0_collide = ((1.0) / (-0.25*((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_8_collide): (xi_10*xi_13 + xi_13*xia_8_collide + xi_5)))) + 2.0));
      const double xia_1_collide = xi_2;
      const double xi_14 = xi_3*xia_1_collide;
      const double xi_15 = xi_14*2.0;
      const double xi_16 = xi_14*8.0;
      const double xi_17 = xi_16 + xi_6;
      const double xi_18 = xi_12 + xi_16;
      const double xi_19 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_17 + xi_15 - xi_17*xia_8_collide): (xi_10*xi_18 + xi_15 + xi_18*xia_8_collide))));
      const double rr_0_a_collide = 0.0;
      const double xi_27 = rr_0_a_collide*xi_26 + xia_23_collide*-0.16666666666666666;
      const double xi_42 = rr_0_a_collide*xi_41 + xia_9_collide*-0.16666666666666666;
      const double xi_51 = rr_0_a_collide*xi_50 + xia_7_collide*-0.16666666666666666;
      const double xi_54 = rr_0_a_collide*0.041666666666666664;
      const double xi_55 = -xi_26 + xi_54*xia_23_collide;
      const double xi_56 = xi_54*xia_9_collide;
      const double xi_57 = xi_41 - xi_56;
      const double xi_58 = xi_55 + xi_57;
      const double xi_70 = -xi_41 + xi_56;
      const double xi_71 = xi_55 + xi_70;
      const double xi_75 = xi_54*xia_7_collide;
      const double xi_76 = -xi_50 + xi_75;
      const double xi_77 = xi_55 + xi_76;
      const double xi_84 = xi_50 + xi_55 - xi_75;
      const double xi_86 = xi_57 + xi_76;
      const double xi_95 = xi_70 + xi_76;
      const double rr_1_a_collide = xi_19;
      const double xi_23 = rr_1_a_collide*0.5;
      const double xi_64 = rr_1_a_collide*0.125;
      const double xi_65 = xi_64*xia_23_collide;
      const double xi_73 = rr_1_a_collide*0.041666666666666664;
      const double xi_188 = rr_1_a_collide*(xia_11_collide*3.0 + xia_16_collide*-3.0 + xia_17_collide*-3.0 + xia_21_collide*-3.0 - xia_26_collide + xia_27_collide*3.0 + xia_28_collide*3.0 + xia_3_collide*-3.0 + xia_5_collide*3.0);
      const double xi_191 = rr_1_a_collide*(xi_190 + xia_10_collide*-5.0 + xia_11_collide*3.0 + xia_25_collide*-5.0 + xia_27_collide*3.0 + xia_28_collide*-2.0 + xia_29_collide*-5.0 + xia_4_collide*-5.0 + xia_5_collide*-2.0);
      const double xi_205 = rr_1_a_collide*xi_204;
      const double xi_282 = xi_205*-0.01984126984126984;
      const double xi_287 = xi_188*0.050000000000000003;
      const double xi_307 = xi_191*0.021428571428571429;
      const double xi_315 = rr_1_a_collide*0.25;
      const double xi_328 = xi_205*0.013888888888888888;
      const double xi_346 = xi_205*-0.003968253968253968;
      const double xi_347 = xi_191*-0.0071428571428571426;
      const double xi_349 = xi_188*0.025000000000000001;
      const double xi_371 = xi_191*0.017857142857142856;
      const double rr_2_a_collide = xia_0_collide*(xi_19*-2.0 + 4.0);
      const double xi_268 = rr_2_a_collide*0.25;
      const double xi_269 = xi_267*xi_268;
      const double xi_276 = rr_2_a_collide*0.083333333333333329;
      const double xi_277 = xi_275*xi_276;
      const double xi_278 = -xi_277;
      const double xi_299 = xi_268*xi_298;
      const double xi_301 = xi_276*xi_300;
      const double xi_304 = -xi_301;
      const double xi_311 = xi_276*xi_310;
      const double xi_312 = xi_268*xi_309 + xi_311;
      const double xi_351 = rr_2_a_collide*0.041666666666666664;
      const double xi_352 = xi_310*xi_351;
      const double xi_353 = rr_2_a_collide*0.125;
      const double xi_354 = xi_309*xi_353;
      const double xi_355 = -xi_352 + xi_354;
      const double xi_356 = xi_275*xi_351;
      const double xi_357 = xi_267*xi_353;
      const double xi_358 = -xi_356 + xi_357;
      const double xi_363 = xi_356 - xi_357;
      const double xi_373 = xi_298*xi_353;
      const double xi_374 = xi_300*xi_351;
      const double xi_375 = -xi_373 + xi_374;
      const double xi_379 = xi_373 - xi_374;
      const double xi_381 = xi_352 - xi_354;
      const double rho_a_collide = xia_6_collide;
      const double u_0_a_collide = xia_22_collide;
      const double xi_20 = u_0_a_collide*xia_9_collide;
      const double xi_32 = xi_20*0.16666666666666666;
      const double xi_33 = xi_20*0.083333333333333329;
      const double xi_34 = rr_1_a_collide*xi_33;
      const double xi_35 = xi_32 - xi_34;
      const double xi_38 = -xi_32 + xi_34;
      const double xi_40 = rr_1_a_collide*xi_32;
      const double xi_61 = u_0_a_collide*xi_60;
      const double xi_66 = u_0_a_collide*xi_65;
      const double xi_88 = u_0_a_collide*xia_7_collide;
      const double xi_89 = xi_88*0.25;
      const double xi_92 = xi_64*xi_88;
      const double xi_192 = u_0_a_collide*u_0_a_collide;
      const double u_1_a_collide = xia_14_collide;
      const double xi_21 = u_1_a_collide*xia_23_collide;
      const double xi_24 = xi_21*0.16666666666666666;
      const double xi_25 = rr_1_a_collide*xi_24;
      const double xi_43 = xi_21*0.083333333333333329;
      const double xi_44 = rr_1_a_collide*xi_43;
      const double xi_45 = -xi_24 + xi_44;
      const double xi_47 = xi_24 - xi_44;
      const double xi_52 = xi_35 + xi_47;
      const double xi_53 = xi_38 + xi_45;
      const double xi_62 = u_1_a_collide*0.25;
      const double xi_63 = xi_62*xia_9_collide;
      const double xi_67 = u_1_a_collide*xi_64;
      const double xi_68 = xi_67*xia_9_collide;
      const double xi_69 = xi_61 + xi_63 - xi_66 - xi_68;
      const double xi_72 = -xi_61 - xi_63 + xi_66 + xi_68;
      const double xi_79 = xi_62*xia_7_collide;
      const double xi_81 = xi_67*xia_7_collide;
      const double xi_193 = rho_a_collide*(u_1_a_collide*u_1_a_collide);
      const double xi_284 = xi_193 + xi_265 - xia_5_collide;
      const double xi_314 = rho_a_collide*u_1_a_collide;
      const double xi_316 = xi_315*(u_0_a_collide*xi_314 + xi_198 + xi_283);
      const double xi_317 = -xi_316;
      const double u_2_a_collide = xia_12_collide;
      const double xi_22 = u_2_a_collide*xia_7_collide;
      const double xi_28 = xi_22*0.16666666666666666;
      const double xi_29 = xi_22*0.083333333333333329;
      const double xi_30 = rr_1_a_collide*xi_29;
      const double xi_31 = xi_28 - xi_30;
      const double xi_36 = xi_31 + xi_35;
      const double xi_37 = -xi_28 + xi_30;
      const double xi_39 = xi_37 + xi_38;
      const double xi_46 = xi_37 + xi_45;
      const double xi_48 = xi_31 + xi_47;
      const double xi_49 = rr_1_a_collide*xi_28;
      const double xi_59 = rr_1_a_collide*u_2_a_collide*xia_7_collide*-0.041666666666666664 + xi_29 + xi_53;
      const double xi_74 = xi_22*xi_73 - xi_29 + xi_52;
      const double xi_78 = rr_1_a_collide*u_0_a_collide*xia_9_collide*-0.041666666666666664 + xi_33 + xi_46;
      const double xi_80 = u_2_a_collide*xi_60;
      const double xi_82 = u_2_a_collide*xi_65;
      const double xi_83 = -xi_79 - xi_80 + xi_81 + xi_82;
      const double xi_85 = xi_20*xi_73 - xi_33 + xi_48;
      const double xi_87 = rr_1_a_collide*u_1_a_collide*xia_23_collide*-0.041666666666666664 + xi_39 + xi_43;
      const double xi_90 = u_2_a_collide*xia_9_collide;
      const double xi_91 = xi_90*0.25;
      const double xi_93 = xi_64*xi_90;
      const double xi_94 = xi_89 + xi_91 - xi_92 - xi_93;
      const double xi_96 = -xi_89 - xi_91 + xi_92 + xi_93;
      const double xi_97 = xi_79 + xi_80 - xi_81 - xi_82;
      const double xi_98 = xi_21*xi_73 + xi_36 - xi_43;
      const double xi_196 = rho_a_collide*(u_2_a_collide*u_2_a_collide);
      const double xi_197 = -xi_196 + xia_10_collide + xia_25_collide + xia_29_collide + xia_4_collide;
      const double xi_200 = rr_1_a_collide*(rho_a_collide*xi_192 - xi_189 + xi_193 - xi_194 - xi_195 - xi_197 - xi_199);
      const double xi_285 = rr_1_a_collide*(xi_197 + xi_283 + xi_284 + xia_11_collide - xia_15_collide - xia_18_collide + xia_27_collide);
      const double xi_286 = xi_285*0.125;
      const double xi_290 = rr_1_a_collide*(rho_a_collide*xi_192*2.0 - xi_196 - xi_199 - xi_284 - xi_288 - xi_289 - xia_10_collide + xia_13_collide*-2.0 + xia_16_collide*2.0 + xia_17_collide*2.0 + xia_21_collide*2.0 + xia_24_collide*-2.0 - xia_25_collide - xia_29_collide + xia_3_collide*2.0 - xia_4_collide);
      const double xi_291 = xi_290*0.041666666666666664;
      const double xi_292 = -xi_291;
      const double xi_293 = -xi_287 + xi_292;
      const double xi_294 = xi_191*0.014285714285714285 + xi_282 + xi_286 + xi_293;
      const double xi_303 = xi_191*-0.035714285714285712 + xi_282 + xi_290*0.083333333333333329;
      const double xi_327 = xi_285*0.0625;
      const double xi_329 = xi_200*0.041666666666666664;
      const double xi_330 = xi_290*0.020833333333333332 + xi_329;
      const double xi_331 = xi_277 + xi_327 + xi_328 + xi_330;
      const double xi_339 = xi_278 + xi_327 + xi_328 + xi_330;
      const double xi_348 = xi_315*(u_2_a_collide*xi_314 + xi_264 + xi_273);
      const double xi_350 = xi_292 + xi_329 + xi_346 + xi_347 + xi_348 + xi_349;
      const double xi_364 = xi_292 + xi_329 + xi_346 + xi_347 - xi_348 + xi_349;
      const double xi_365 = xi_315*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_296 - xia_25_collide + xia_29_collide);
      const double xi_366 = -xi_365;
      const double xi_370 = -xi_327;
      const double xi_372 = xi_311 + xi_330 + xi_346 + xi_370 + xi_371;
      const double xi_384 = -xi_311 + xi_330 + xi_346 + xi_370 + xi_371;
      const double forceTerm_0_a_collide = xi_20*xi_23 - xi_20 + xi_21*xi_23 - xi_21 + xi_22*xi_23 - xi_22;
      const double forceTerm_1_a_collide = u_1_a_collide*xia_23_collide*0.33333333333333331 - xi_25 - xi_27 - xi_36;
      const double forceTerm_2_a_collide = xi_21*0.33333333333333331 - xi_25 + xi_27 + xi_39;
      const double forceTerm_3_a_collide = xi_20*0.33333333333333331 - xi_40 + xi_42 + xi_46;
      const double forceTerm_4_a_collide = u_0_a_collide*xia_9_collide*0.33333333333333331 - xi_40 - xi_42 - xi_48;
      const double forceTerm_5_a_collide = u_2_a_collide*xia_7_collide*0.33333333333333331 - xi_49 - xi_51 - xi_52;
      const double forceTerm_6_a_collide = xi_22*0.33333333333333331 - xi_49 + xi_51 + xi_53;
      const double forceTerm_7_a_collide = -xi_58 - xi_59 - xi_69;
      const double forceTerm_8_a_collide = -xi_59 - xi_71 - xi_72;
      const double forceTerm_9_a_collide = xi_69 + xi_71 + xi_74;
      const double forceTerm_10_a_collide = xi_58 + xi_72 + xi_74;
      const double forceTerm_11_a_collide = -xi_77 - xi_78 - xi_83;
      const double forceTerm_12_a_collide = xi_83 + xi_84 + xi_85;
      const double forceTerm_13_a_collide = -xi_86 - xi_87 - xi_94;
      const double forceTerm_14_a_collide = -xi_87 - xi_95 - xi_96;
      const double forceTerm_15_a_collide = -xi_78 - xi_84 - xi_97;
      const double forceTerm_16_a_collide = xi_77 + xi_85 + xi_97;
      const double forceTerm_17_a_collide = xi_94 + xi_95 + xi_98;
      const double forceTerm_18_a_collide = xi_86 + xi_96 + xi_98;
      const double xib_3_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
      const double xib_4_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
      const double xi_402 = -xib_4_collide;
      const double xib_5_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
      const double xib_6_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
      const double xib_8_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
      const double xi_425 = -xib_8_collide;
      const double xib_9_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
      const double xib_10_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
      const double xib_11_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_99 = xib_11_collide > 0.5;
      const bool xi_100 = xib_11_collide < -0.5;
      const double xi_101 = xib_11_collide*xib_11_collide;
      const bool xi_102 = xib_11_collide > 0.0;
      const double xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
      const double xi_440 = xi_425 + xib_12_collide;
      const double xib_13_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const double xib_14_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const double xib_15_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
      const double xi_399 = xib_15_collide*2.0 + xib_3_collide*2.0;
      const double xi_431 = xib_15_collide - xib_3_collide;
      const double xib_16_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
      const double xi_393 = xib_16_collide + xib_4_collide;
      const double xi_412 = -xib_16_collide;
      const double xi_413 = xi_412 + xib_4_collide;
      const double xib_17_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
      const double xi_410 = -xib_17_collide;
      const double xib_18_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
      const double xi_139 = xib_18_collide*0.083333333333333329;
      const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
      const double xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
      const double xi_387 = -xib_20_collide;
      const double xib_21_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
      const double xi_405 = -xib_21_collide;
      const double xi_406 = xi_405 + xib_7_collide;
      const double xib_22_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const double xib_23_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
      const double xi_396 = xib_10_collide + xib_23_collide;
      const double xi_414 = xib_10_collide*2.0 + xib_23_collide*-2.0;
      const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
      const double xi_388 = xi_387 + xib_16_collide*2.0 + xib_17_collide*2.0 + xib_19_collide*5.0 + xib_24_collide*2.0 + xib_4_collide*2.0 + xib_9_collide*5.0;
      const double xi_392 = xib_17_collide + xib_24_collide;
      const double xi_403 = -xib_24_collide;
      const double xi_404 = xi_403 + xib_17_collide;
      const double xi_407 = xi_402 + xi_404 + xi_406 + xib_16_collide;
      const double xi_441 = -xi_392 - xi_402 - xi_412 - xi_440;
      const double xib_25_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
      const double xi_115 = xib_25_collide*0.083333333333333329;
      const double xi_150 = xib_25_collide*0.25;
      const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
      const double xib_27_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
      const double xi_397 = xi_396 + xib_26_collide + xib_27_collide;
      const double xi_411 = xib_27_collide*2.0;
      const double xi_415 = xi_406 + xi_410 - xi_411 + xi_413 + xi_414 + xib_24_collide + xib_26_collide*2.0;
      const double xi_420 = -xib_26_collide - xib_27_collide;
      const double xib_28_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
      const double xi_400 = xi_388 + xi_399 + xib_10_collide*-7.0 + xib_12_collide*-4.0 + xib_21_collide*5.0 + xib_23_collide*-7.0 + xib_26_collide*-7.0 + xib_27_collide*-7.0 + xib_28_collide*2.0 + xib_5_collide*2.0 + xib_7_collide*5.0 + xib_8_collide*-4.0;
      const double xi_432 = xi_431 + xib_28_collide - xib_5_collide;
      const double xi_433 = xi_432 - xib_19_collide + xib_9_collide;
      const double xi_435 = -xi_411 - xi_414 - xi_432 - xib_19_collide + xib_26_collide*2.0 + xib_9_collide;
      const double xi_442 = -xi_393 - xi_399 - xi_403 - xi_410 - xi_440 + xib_28_collide*2.0 + xib_5_collide*2.0;
      const double xib_29_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
      const double xi_130 = xib_29_collide*0.083333333333333329;
      const double xib_0_collide = ((1.0) / (-0.25*((xi_99) ? (omega_shear_a): ((xi_100) ? (omega_shear_b): ((xi_102) ? (xi_101*xi_8 + xi_5 + xi_9*xib_11_collide): (xi_101*xi_13 + xi_13*xib_11_collide + xi_5)))) + 2.0));
      const double xib_1_collide = xi_2;
      const double xi_103 = xi_3*xib_1_collide;
      const double xi_104 = xi_103*2.0;
      const double xi_105 = xi_103*8.0;
      const double xi_106 = xi_105 + xi_6;
      const double xi_107 = xi_105 + xi_12;
      const double xi_108 = ((xi_99) ? (omega_shear_a): ((xi_100) ? (omega_shear_b): ((xi_102) ? (xi_101*xi_106 + xi_104 - xi_106*xib_11_collide): (xi_101*xi_107 + xi_104 + xi_107*xib_11_collide))));
      const double rr_0_b_collide = 0.0;
      const double xi_116 = rr_0_b_collide*xi_115 + xib_25_collide*-0.16666666666666666;
      const double xi_131 = rr_0_b_collide*xi_130 + xib_29_collide*-0.16666666666666666;
      const double xi_140 = rr_0_b_collide*xi_139 + xib_18_collide*-0.16666666666666666;
      const double xi_143 = rr_0_b_collide*0.041666666666666664;
      const double xi_144 = xi_143*xib_25_collide;
      const double xi_145 = -xi_115 + xi_144;
      const double xi_146 = xi_143*xib_29_collide;
      const double xi_147 = xi_130 - xi_146;
      const double xi_148 = xi_145 + xi_147;
      const double xi_160 = -xi_130 + xi_146;
      const double xi_161 = xi_145 + xi_160;
      const double xi_165 = -xi_139 + xi_143*xib_18_collide;
      const double xi_166 = xi_145 + xi_165;
      const double xi_173 = xi_115 - xi_144 + xi_165;
      const double xi_175 = xi_147 + xi_165;
      const double xi_184 = xi_160 + xi_165;
      const double rr_1_b_collide = xi_108;
      const double xi_112 = rr_1_b_collide*0.5;
      const double xi_154 = rr_1_b_collide*0.125;
      const double xi_155 = xi_154*xib_25_collide;
      const double xi_163 = rr_1_b_collide*0.041666666666666664;
      const double xi_386 = rr_1_b_collide*(xib_12_collide*3.0 + xib_16_collide*-3.0 + xib_17_collide*-3.0 - xib_20_collide + xib_21_collide*3.0 + xib_24_collide*-3.0 + xib_4_collide*-3.0 + xib_7_collide*3.0 + xib_8_collide*3.0);
      const double xi_389 = rr_1_b_collide*(xi_388 + xib_12_collide*3.0 + xib_15_collide*-5.0 + xib_21_collide*-2.0 + xib_28_collide*-5.0 + xib_3_collide*-5.0 + xib_5_collide*-5.0 + xib_7_collide*-2.0 + xib_8_collide*3.0);
      const double xi_401 = rr_1_b_collide*xi_400;
      const double xi_419 = xi_401*-0.01984126984126984;
      const double xi_424 = xi_386*0.050000000000000003;
      const double xi_439 = xi_389*0.021428571428571429;
      const double xi_446 = rr_1_b_collide*0.25;
      const double xi_450 = xi_401*0.013888888888888888;
      const double xi_455 = xi_401*-0.003968253968253968;
      const double xi_456 = xi_389*-0.0071428571428571426;
      const double xi_458 = xi_386*0.025000000000000001;
      const double xi_473 = xi_389*0.017857142857142856;
      const double rr_2_b_collide = xib_0_collide*(xi_108*-2.0 + 4.0);
      const double xi_408 = rr_2_b_collide*0.25;
      const double xi_409 = xi_407*xi_408;
      const double xi_416 = rr_2_b_collide*0.083333333333333329;
      const double xi_417 = xi_415*xi_416;
      const double xi_418 = -xi_417;
      const double xi_434 = xi_408*xi_433;
      const double xi_436 = xi_416*xi_435;
      const double xi_438 = -xi_436;
      const double xi_443 = xi_416*xi_442;
      const double xi_444 = xi_408*xi_441 + xi_443;
      const double xi_460 = rr_2_b_collide*0.041666666666666664;
      const double xi_461 = xi_442*xi_460;
      const double xi_462 = rr_2_b_collide*0.125;
      const double xi_463 = xi_441*xi_462;
      const double xi_464 = -xi_461 + xi_463;
      const double xi_465 = xi_415*xi_460;
      const double xi_466 = xi_407*xi_462;
      const double xi_467 = -xi_465 + xi_466;
      const double xi_468 = xi_465 - xi_466;
      const double xi_475 = xi_433*xi_462;
      const double xi_476 = xi_435*xi_460;
      const double xi_477 = -xi_475 + xi_476;
      const double xi_478 = xi_475 - xi_476;
      const double xi_479 = xi_461 - xi_463;
      const double rho_b_collide = xib_6_collide;
      const double u_0_b_collide = xib_22_collide;
      const double xi_109 = u_0_b_collide*xib_29_collide;
      const double xi_121 = xi_109*0.16666666666666666;
      const double xi_122 = xi_109*0.083333333333333329;
      const double xi_123 = rr_1_b_collide*xi_122;
      const double xi_124 = xi_121 - xi_123;
      const double xi_127 = -xi_121 + xi_123;
      const double xi_129 = rr_1_b_collide*xi_121;
      const double xi_151 = u_0_b_collide*xi_150;
      const double xi_156 = u_0_b_collide*xi_155;
      const double xi_177 = u_0_b_collide*xib_18_collide;
      const double xi_178 = xi_177*0.25;
      const double xi_181 = xi_154*xi_177;
      const double xi_390 = u_0_b_collide*u_0_b_collide;
      const double u_1_b_collide = xib_14_collide;
      const double xi_110 = u_1_b_collide*xib_25_collide;
      const double xi_113 = xi_110*0.16666666666666666;
      const double xi_114 = rr_1_b_collide*xi_113;
      const double xi_132 = xi_110*0.083333333333333329;
      const double xi_133 = rr_1_b_collide*xi_132;
      const double xi_134 = -xi_113 + xi_133;
      const double xi_136 = xi_113 - xi_133;
      const double xi_141 = xi_124 + xi_136;
      const double xi_142 = xi_127 + xi_134;
      const double xi_152 = u_1_b_collide*0.25;
      const double xi_153 = xi_152*xib_29_collide;
      const double xi_157 = u_1_b_collide*xi_154;
      const double xi_158 = xi_157*xib_29_collide;
      const double xi_159 = xi_151 + xi_153 - xi_156 - xi_158;
      const double xi_162 = -xi_151 - xi_153 + xi_156 + xi_158;
      const double xi_168 = xi_152*xib_18_collide;
      const double xi_170 = xi_157*xib_18_collide;
      const double xi_391 = rho_b_collide*(u_1_b_collide*u_1_b_collide);
      const double xi_421 = xi_391 + xi_405 - xib_7_collide;
      const double xi_445 = rho_b_collide*u_1_b_collide;
      const double xi_447 = xi_446*(u_0_b_collide*xi_445 + xi_396 + xi_420);
      const double xi_448 = -xi_447;
      const double u_2_b_collide = xib_13_collide;
      const double xi_111 = u_2_b_collide*xib_18_collide;
      const double xi_117 = xi_111*0.16666666666666666;
      const double xi_118 = xi_111*0.083333333333333329;
      const double xi_119 = rr_1_b_collide*xi_118;
      const double xi_120 = xi_117 - xi_119;
      const double xi_125 = xi_120 + xi_124;
      const double xi_126 = -xi_117 + xi_119;
      const double xi_128 = xi_126 + xi_127;
      const double xi_135 = xi_126 + xi_134;
      const double xi_137 = xi_120 + xi_136;
      const double xi_138 = rr_1_b_collide*xi_117;
      const double xi_149 = rr_1_b_collide*u_2_b_collide*xib_18_collide*-0.041666666666666664 + xi_118 + xi_142;
      const double xi_164 = xi_111*xi_163 - xi_118 + xi_141;
      const double xi_167 = rr_1_b_collide*u_0_b_collide*xib_29_collide*-0.041666666666666664 + xi_122 + xi_135;
      const double xi_169 = u_2_b_collide*xi_150;
      const double xi_171 = u_2_b_collide*xi_155;
      const double xi_172 = -xi_168 - xi_169 + xi_170 + xi_171;
      const double xi_174 = xi_168 + xi_169 - xi_170 - xi_171;
      const double xi_176 = rr_1_b_collide*u_1_b_collide*xib_25_collide*-0.041666666666666664 + xi_128 + xi_132;
      const double xi_179 = u_2_b_collide*xib_29_collide;
      const double xi_180 = xi_179*0.25;
      const double xi_182 = xi_154*xi_179;
      const double xi_183 = xi_178 + xi_180 - xi_181 - xi_182;
      const double xi_185 = -xi_178 - xi_180 + xi_181 + xi_182;
      const double xi_186 = xi_109*xi_163 - xi_122 + xi_137;
      const double xi_187 = xi_110*xi_163 + xi_125 - xi_132;
      const double xi_394 = rho_b_collide*(u_2_b_collide*u_2_b_collide);
      const double xi_395 = -xi_394 + xib_15_collide + xib_28_collide + xib_3_collide + xib_5_collide;
      const double xi_398 = rr_1_b_collide*(rho_b_collide*xi_390 - xi_387 + xi_391 - xi_392 - xi_393 - xi_395 - xi_397);
      const double xi_422 = rr_1_b_collide*(xi_395 + xi_420 + xi_421 - xib_10_collide + xib_12_collide - xib_23_collide + xib_8_collide);
      const double xi_423 = xi_422*0.125;
      const double xi_426 = rr_1_b_collide*(rho_b_collide*xi_390*2.0 - xi_394 - xi_397 - xi_421 - xi_425 + xib_12_collide - xib_15_collide + xib_16_collide*2.0 + xib_17_collide*2.0 + xib_19_collide*-2.0 + xib_24_collide*2.0 - xib_28_collide - xib_3_collide + xib_4_collide*2.0 - xib_5_collide + xib_9_collide*-2.0);
      const double xi_427 = xi_426*0.041666666666666664;
      const double xi_428 = -xi_427;
      const double xi_429 = -xi_424 + xi_428;
      const double xi_430 = xi_389*0.014285714285714285 + xi_419 + xi_423 + xi_429;
      const double xi_437 = xi_389*-0.035714285714285712 + xi_419 + xi_426*0.083333333333333329;
      const double xi_449 = xi_422*0.0625;
      const double xi_451 = xi_398*0.041666666666666664;
      const double xi_452 = xi_426*0.020833333333333332 + xi_451;
      const double xi_453 = xi_417 + xi_449 + xi_450 + xi_452;
      const double xi_454 = xi_418 + xi_449 + xi_450 + xi_452;
      const double xi_457 = xi_446*(u_2_b_collide*xi_445 + xi_404 + xi_413);
      const double xi_459 = xi_428 + xi_451 + xi_455 + xi_456 + xi_457 + xi_458;
      const double xi_469 = xi_428 + xi_451 + xi_455 + xi_456 - xi_457 + xi_458;
      const double xi_470 = xi_446*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_431 - xib_28_collide + xib_5_collide);
      const double xi_471 = -xi_470;
      const double xi_472 = -xi_449;
      const double xi_474 = xi_443 + xi_452 + xi_455 + xi_472 + xi_473;
      const double xi_480 = -xi_443 + xi_452 + xi_455 + xi_472 + xi_473;
      const double forceTerm_0_b_collide = xi_109*xi_112 - xi_109 + xi_110*xi_112 - xi_110 + xi_111*xi_112 - xi_111;
      const double forceTerm_1_b_collide = u_1_b_collide*xib_25_collide*0.33333333333333331 - xi_114 - xi_116 - xi_125;
      const double forceTerm_2_b_collide = xi_110*0.33333333333333331 - xi_114 + xi_116 + xi_128;
      const double forceTerm_3_b_collide = xi_109*0.33333333333333331 - xi_129 + xi_131 + xi_135;
      const double forceTerm_4_b_collide = u_0_b_collide*xib_29_collide*0.33333333333333331 - xi_129 - xi_131 - xi_137;
      const double forceTerm_5_b_collide = u_2_b_collide*xib_18_collide*0.33333333333333331 - xi_138 - xi_140 - xi_141;
      const double forceTerm_6_b_collide = xi_111*0.33333333333333331 - xi_138 + xi_140 + xi_142;
      const double forceTerm_7_b_collide = -xi_148 - xi_149 - xi_159;
      const double forceTerm_8_b_collide = -xi_149 - xi_161 - xi_162;
      const double forceTerm_9_b_collide = xi_159 + xi_161 + xi_164;
      const double forceTerm_10_b_collide = xi_148 + xi_162 + xi_164;
      const double forceTerm_11_b_collide = -xi_166 - xi_167 - xi_172;
      const double forceTerm_12_b_collide = -xi_167 - xi_173 - xi_174;
      const double forceTerm_13_b_collide = -xi_175 - xi_176 - xi_183;
      const double forceTerm_14_b_collide = -xi_176 - xi_184 - xi_185;
      const double forceTerm_15_b_collide = xi_172 + xi_173 + xi_186;
      const double forceTerm_16_b_collide = xi_166 + xi_174 + xi_186;
      const double forceTerm_17_b_collide = xi_183 + xi_184 + xi_187;
      const double forceTerm_18_b_collide = xi_175 + xi_185 + xi_187;
      const double tmp_a0 = forceTerm_0_a_collide + xi_188*0.10000000000000001 + xi_191*0.042857142857142858 + xi_200*-0.5 + xi_205*0.023809523809523808 + xi_261 + xia_26_collide;
      const double tmp_a1 = forceTerm_1_a_collide - xi_269 + xi_278 + xi_281 + xi_294 + xia_28_collide;
      const double tmp_a2 = forceTerm_2_a_collide + xi_269 + xi_277 + xi_294 + xi_295 + xia_5_collide;
      const double tmp_a3 = forceTerm_3_a_collide + xi_299 + xi_301 + xi_302 + xi_303 + xia_13_collide;
      const double tmp_a4 = forceTerm_4_a_collide - xi_299 + xi_303 + xi_304 + xi_305 + xia_24_collide;
      const double tmp_a5 = forceTerm_5_a_collide + rr_1_a_collide*xi_204*0.015873015873015872 - xi_286 - xi_287 - xi_289 - xi_291 - xi_306 - xi_307 - xi_312;
      const double tmp_a6 = forceTerm_6_a_collide + xi_205*0.015873015873015872 - xi_286 + xi_293 - xi_307 + xi_312 + xi_313 + xia_27_collide;
      const double tmp_a7 = forceTerm_7_a_collide + xi_304 + xi_317 + xi_326 + xi_331 + xia_15_collide;
      const double tmp_a8 = forceTerm_8_a_collide + xi_301 + xi_316 + xi_331 + xi_337 + xia_19_collide;
      const double tmp_a9 = forceTerm_9_a_collide + xi_304 + xi_316 + xi_338 + xi_339 + xia_20_collide;
      const double tmp_a10 = forceTerm_10_a_collide + xi_301 + xi_317 + xi_339 + xi_340 + xia_18_collide;
      const double tmp_a11 = forceTerm_11_a_collide + xi_345 + xi_350 + xi_355 + xi_358 + xia_3_collide;
      const double tmp_a12 = forceTerm_12_a_collide + xi_355 + xi_362 + xi_363 + xi_364 + xia_21_collide;
      const double tmp_a13 = forceTerm_13_a_collide + xi_366 + xi_369 + xi_372 + xi_375 + xia_29_collide;
      const double tmp_a14 = forceTerm_14_a_collide + xi_365 + xi_372 + xi_378 + xi_379 + xia_25_collide;
      const double tmp_a15 = forceTerm_15_a_collide + xi_358 + xi_364 + xi_380 + xi_381 + xia_16_collide;
      const double tmp_a16 = forceTerm_16_a_collide + xi_350 + xi_363 + xi_381 + xi_382 + xia_17_collide;
      const double tmp_a17 = forceTerm_17_a_collide + xi_365 + xi_375 + xi_383 + xi_384 + xia_4_collide;
      const double tmp_a18 = forceTerm_18_a_collide + xi_366 + xi_379 + xi_384 + xi_385 + xia_10_collide;
      const double tmp_b0 = forceTerm_0_b_collide + xi_261 + xi_386*0.10000000000000001 + xi_389*0.042857142857142858 + xi_398*-0.5 + xi_401*0.023809523809523808 + xib_20_collide;
      const double tmp_b1 = forceTerm_1_b_collide + xi_281 - xi_409 + xi_418 + xi_430 + xib_7_collide;
      const double tmp_b2 = forceTerm_2_b_collide + xi_295 + xi_409 + xi_417 + xi_430 + xib_21_collide;
      const double tmp_b3 = forceTerm_3_b_collide + xi_302 + xi_434 + xi_436 + xi_437 + xib_19_collide;
      const double tmp_b4 = forceTerm_4_b_collide + xi_305 - xi_434 + xi_437 + xi_438 + xib_9_collide;
      const double tmp_b5 = forceTerm_5_b_collide + rr_1_b_collide*xi_400*0.015873015873015872 - xi_306 - xi_423 - xi_424 - xi_425 - xi_427 - xi_439 - xi_444;
      const double tmp_b6 = forceTerm_6_b_collide + xi_313 + xi_401*0.015873015873015872 - xi_423 + xi_429 - xi_439 + xi_444 + xib_12_collide;
      const double tmp_b7 = forceTerm_7_b_collide + xi_326 + xi_438 + xi_448 + xi_453 + xib_23_collide;
      const double tmp_b8 = forceTerm_8_b_collide + xi_337 + xi_436 + xi_447 + xi_453 + xib_27_collide;
      const double tmp_b9 = forceTerm_9_b_collide + xi_338 + xi_438 + xi_447 + xi_454 + xib_26_collide;
      const double tmp_b10 = forceTerm_10_b_collide + xi_340 + xi_436 + xi_448 + xi_454 + xib_10_collide;
      const double tmp_b11 = forceTerm_11_b_collide + xi_345 + xi_459 + xi_464 + xi_467 + xib_24_collide;
      const double tmp_b12 = forceTerm_12_b_collide + xi_362 + xi_464 + xi_468 + xi_469 + xib_17_collide;
      const double tmp_b13 = forceTerm_13_b_collide + xi_369 + xi_471 + xi_474 + xi_477 + xib_15_collide;
      const double tmp_b14 = forceTerm_14_b_collide + xi_378 + xi_470 + xi_474 + xi_478 + xib_3_collide;
      const double tmp_b15 = forceTerm_15_b_collide + xi_380 + xi_467 + xi_469 + xi_479 + xib_4_collide;
      const double tmp_b16 = forceTerm_16_b_collide + xi_382 + xi_459 + xi_468 + xi_479 + xib_16_collide;
      const double tmp_b17 = forceTerm_17_b_collide + xi_383 + xi_470 + xi_477 + xi_480 + xib_28_collide;
      const double tmp_b18 = forceTerm_18_b_collide + xi_385 + xi_471 + xi_478 + xi_480 + xib_5_collide;
      const double xirecolor_0 = tmp_a0 + tmp_b0;
      const double xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_2 = ((1.0) / (xirecolor_1));
      const double xi_502 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const double xirecolor_4 = tmp_a1 + tmp_b1;
      const double xirecolor_5 = xi_216;
      const double xirecolor_6 = xi_214;
      const double xirecolor_7 = xi_234;
      const double xi_481 = -xirecolor_7;
      const double xirecolor_8 = xi_481;
      const double xirecolor_9 = xi_221;
      const double xirecolor_10 = -xirecolor_9;
      const double xirecolor_11 = xi_223;
      const double xi_482 = xirecolor_10 + xirecolor_11;
      const double xirecolor_12 = xi_236;
      const double xi_483 = xirecolor_12 + xirecolor_8;
      const double xirecolor_13 = xi_482 + xi_483;
      const double xirecolor_14 = xi_238;
      const double xirecolor_15 = xi_239;
      const double xi_484 = xirecolor_14 - xirecolor_15;
      const double xirecolor_16 = -xi_484;
      const double xirecolor_17 = xi_208;
      const double xirecolor_18 = xi_209;
      const double xirecolor_19 = xi_210;
      const double xi_485 = xirecolor_18 + xirecolor_19;
      const double xirecolor_20 = xi_211;
      const double xirecolor_21 = xi_485 - xirecolor_17 + xirecolor_20;
      const double xirecolor_22 = xi_207;
      const double xirecolor_23 = xi_206;
      const double xi_486 = xirecolor_22 - xirecolor_23;
      const double xirecolor_24 = -xi_486;
      const double xirecolor_25 = xi_242;
      const double xirecolor_26 = -xirecolor_25;
      const double xirecolor_27 = xi_228;
      const double xirecolor_28 = xi_227;
      const double xirecolor_29 = xi_244;
      const double xi_487 = xirecolor_26 + xirecolor_29;
      const double xirecolor_30 = xi_487 - xirecolor_27 + xirecolor_28;
      const double xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
      const double xirecolor_32 = xi_213;
      const double xirecolor_33 = xi_215;
      const double xi_488 = xirecolor_32 + xirecolor_33;
      const double xirecolor_34 = xi_488 + xirecolor_5 - xirecolor_6;
      const double xirecolor_35 = xi_218;
      const double xirecolor_36 = xi_219;
      const double xi_489 = xirecolor_35 - xirecolor_36;
      const double xirecolor_37 = -xi_489;
      const double xirecolor_38 = xi_482;
      const double xirecolor_39 = xi_225;
      const double xirecolor_40 = xi_226;
      const double xi_490 = xirecolor_39 + xirecolor_40;
      const double xirecolor_41 = xirecolor_27 - xirecolor_28;
      const double xirecolor_42 = xi_490 + xirecolor_41;
      const double xirecolor_43 = xirecolor_38 + xirecolor_42;
      const double xirecolor_44 = xi_486 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
      const double xirecolor_45 = xi_484;
      const double xirecolor_46 = xi_483;
      const double xirecolor_47 = xi_250;
      const double xirecolor_48 = -xirecolor_47;
      const double xirecolor_49 = xi_252;
      const double xi_491 = xirecolor_48 + xirecolor_49;
      const double xirecolor_50 = xi_489 + xi_491;
      const double xirecolor_51 = xi_485 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
      const double xirecolor_52 = pow(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5);
      const double xirecolor_53 = ((1.0) / (xirecolor_52));
      const double xi_492 = xirecolor_31*xirecolor_53;
      const double xi_493 = xirecolor_51*xirecolor_53;
      const double xi_494 = xirecolor_44*xirecolor_53;
      const bool xirecolor_54 = xirecolor_52 > 0.0;
      const double xirecolor_55 = beta*((1.0) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const double xirecolor_56 = xirecolor_55*(0.055555555555555552*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.055555555555555552*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const double xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_492): (0.0));
      const double xirecolor_58 = tmp_a2 + tmp_b2;
      const double xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_492): (0.0));
      const double xirecolor_60 = tmp_a3 + tmp_b3;
      const double xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_493): (0.0));
      const double xirecolor_62 = tmp_a4 + tmp_b4;
      const double xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_493): (0.0));
      const double xirecolor_64 = tmp_a5 + tmp_b5;
      const double xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_494): (0.0));
      const double xirecolor_66 = tmp_a6 + tmp_b6;
      const double xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_494): (0.0));
      const double xirecolor_68 = tmp_a7 + tmp_b7;
      const double xirecolor_69 = xi_318;
      const double xirecolor_70 = xi_319;
      const double xirecolor_71 = xi_488;
      const double xirecolor_72 = xirecolor_50 + xirecolor_71;
      const double xirecolor_73 = xi_323;
      const double xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
      const double xirecolor_75 = xi_320 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
      const double xirecolor_76 = xirecolor_53*0.70710678118654757;
      const double xi_495 = xirecolor_75*xirecolor_76;
      const double xirecolor_77 = xirecolor_55*(0.027777777777777776*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.027777777777777776*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const double xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_495): (0.0));
      const double xirecolor_79 = tmp_a8 + tmp_b8;
      const double xirecolor_80 = xi_332;
      const double xirecolor_81 = xi_333;
      const double xi_496 = -xirecolor_80 + xirecolor_81;
      const double xirecolor_82 = xi_334;
      const double xirecolor_83 = xi_335 + xi_496 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
      const double xi_497 = xirecolor_76*xirecolor_83;
      const double xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_497): (0.0));
      const double xirecolor_85 = tmp_a9 + tmp_b9;
      const double xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_497): (0.0));
      const double xirecolor_87 = tmp_a10 + tmp_b10;
      const double xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_495): (0.0));
      const double xirecolor_89 = tmp_a11 + tmp_b11;
      const double xirecolor_90 = xi_490 + xirecolor_37 + xirecolor_71;
      const double xirecolor_91 = -xi_481 - xi_496 - xirecolor_12;
      const double xirecolor_92 = xi_341 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
      const double xi_498 = xirecolor_76*xirecolor_92;
      const double xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_498): (0.0));
      const double xirecolor_94 = tmp_a12 + tmp_b12;
      const double xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
      const double xirecolor_96 = xi_359 + xi_487 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
      const double xi_499 = xirecolor_76*xirecolor_96;
      const double xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_499): (0.0));
      const double xirecolor_98 = tmp_a13 + tmp_b13;
      const double xirecolor_99 = xi_367 + xi_491 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
      const double xi_500 = xirecolor_76*xirecolor_99;
      const double xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_500): (0.0));
      const double xirecolor_101 = tmp_a14 + tmp_b14;
      const double xirecolor_102 = xi_376 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
      const double xi_501 = xirecolor_102*xirecolor_76;
      const double xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_501): (0.0));
      const double xirecolor_104 = tmp_a15 + tmp_b15;
      const double xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_499): (0.0));
      const double xirecolor_106 = tmp_a16 + tmp_b16;
      const double xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_498): (0.0));
      const double xirecolor_108 = tmp_a17 + tmp_b17;
      const double xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_501): (0.0));
      const double xirecolor_110 = tmp_a18 + tmp_b18;
      const double xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_500): (0.0));
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2] = xirecolor_0*xirecolor_3;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3] = xirecolor_3*xirecolor_4 + xirecolor_57;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3] = xirecolor_3*xirecolor_58 + xirecolor_59;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3] = xirecolor_3*xirecolor_60 + xirecolor_61;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3] = xirecolor_3*xirecolor_62 + xirecolor_63;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3] = xirecolor_3*xirecolor_64 + xirecolor_65;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3] = xirecolor_3*xirecolor_66 + xirecolor_67;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3] = xirecolor_3*xirecolor_68 + xirecolor_78;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3] = xirecolor_3*xirecolor_79 + xirecolor_84;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3] = xirecolor_3*xirecolor_85 + xirecolor_86;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3] = xirecolor_3*xirecolor_87 + xirecolor_88;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3] = xirecolor_3*xirecolor_89 + xirecolor_93;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3] = xirecolor_3*xirecolor_94 + xirecolor_97;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3] = xirecolor_100 + xirecolor_3*xirecolor_98;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3] = xirecolor_101*xirecolor_3 + xirecolor_103;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3] = xirecolor_104*xirecolor_3 + xirecolor_105;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3] = xirecolor_106*xirecolor_3 + xirecolor_107;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3] = xirecolor_108*xirecolor_3 + xirecolor_109;
      _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3] = xirecolor_110*xirecolor_3 + xirecolor_111;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = xi_502*xirecolor_0;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = xi_502*xirecolor_4 - xirecolor_57;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = xi_502*xirecolor_58 - xirecolor_59;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = xi_502*xirecolor_60 - xirecolor_61;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = xi_502*xirecolor_62 - xirecolor_63;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = xi_502*xirecolor_64 - xirecolor_65;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = xi_502*xirecolor_66 - xirecolor_67;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = xi_502*xirecolor_68 - xirecolor_78;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = xi_502*xirecolor_79 - xirecolor_84;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = xi_502*xirecolor_85 - xirecolor_86;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = xi_502*xirecolor_87 - xirecolor_88;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = xi_502*xirecolor_89 - xirecolor_93;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = xi_502*xirecolor_94 - xirecolor_97;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = xi_502*xirecolor_98 - xirecolor_100;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = xi_502*xirecolor_101 - xirecolor_103;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = xi_502*xirecolor_104 - xirecolor_105;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = xi_502*xirecolor_106 - xirecolor_107;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = xi_502*xirecolor_108 - xirecolor_109;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = xi_502*xirecolor_110 - xirecolor_111;
   } 
}
}


void ColorGradientCollideSweepDoublePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   
    auto pdfs_a = block->getData< gpu::GPUField<double> >(pdfs_aID);
    auto rho_a = block->getData< gpu::GPUField<double> >(rho_aID);
    auto velocity = block->getData< gpu::GPUField<double> >(velocityID);
    auto rho_b = block->getData< gpu::GPUField<double> >(rho_bID);
    auto pdfs_b = block->getData< gpu::GPUField<double> >(pdfs_bID);
    auto force_a = block->getData< gpu::GPUField<double> >(force_aID);
    auto force_b = block->getData< gpu::GPUField<double> >(force_bID);
    auto phasefield = block->getData< gpu::GPUField<double> >(phasefieldID);

    auto & beta = this->beta_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & sigma = this->sigma_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda::colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepDoublePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto pdfs_a = block->getData< gpu::GPUField<double> >(pdfs_aID);
    auto rho_a = block->getData< gpu::GPUField<double> >(rho_aID);
    auto velocity = block->getData< gpu::GPUField<double> >(velocityID);
    auto rho_b = block->getData< gpu::GPUField<double> >(rho_bID);
    auto pdfs_b = block->getData< gpu::GPUField<double> >(pdfs_bID);
    auto force_a = block->getData< gpu::GPUField<double> >(force_aID);
    auto force_b = block->getData< gpu::GPUField<double> >(force_bID);
    auto phasefield = block->getData< gpu::GPUField<double> >(phasefieldID);

    auto & beta = this->beta_;
    auto & omega_shear_a = this->omega_shear_a_;
    auto & omega_shear_b = this->omega_shear_b_;
    auto & sigma = this->sigma_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    const int64_t _stride_force_a_0 = int64_t(force_a->xStride());
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_0 = int64_t(force_b->xStride());
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_0 = int64_t(pdfs_a->xStride());
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_0 = int64_t(pdfs_b->xStride());
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_0 = int64_t(phasefield->xStride());
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_0 = int64_t(rho_a->xStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_0 = int64_t(rho_b->xStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_0 = int64_t(velocity->xStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    dim3 _block(uint32_c(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)), uint32_c(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))), uint32_c(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))));
    dim3 _grid(uint32_c(( (_size_force_a_0 - 2) % (((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) == 0 ? (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) : ( (int64_t)(_size_force_a_0 - 2) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)) ) +1 )), uint32_c(( (_size_force_a_1 - 2) % (((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) == 0 ? (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) : ( (int64_t)(_size_force_a_1 - 2) / (int64_t)(((1024 < ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))) ? 1024 : ((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))) ) +1 )), uint32_c(( (_size_force_a_2 - 2) % (((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) == 0 ? (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) : ( (int64_t)(_size_force_a_2 - 2) / (int64_t)(((64 < ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))))))) ? 64 : ((_size_force_a_2 - 2 < ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))) ? _size_force_a_2 - 2 : ((int64_t)(256) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)*((_size_force_a_1 - 2 < 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2)))) ? _size_force_a_1 - 2 : 2*((int64_t)(128) / (int64_t)(((128 < _size_force_a_0 - 2) ? 128 : _size_force_a_0 - 2))))))))) ) +1 )));
    internal_colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda::colorgradientcollidesweepdoubleprecisioncuda_colorgradientcollidesweepdoubleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
