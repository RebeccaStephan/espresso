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
//! \\file ColorGradientCollideSweepSinglePrecisionCUDA.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepSinglePrecisionCUDA.h"




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


namespace internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda {
static FUNC_PREFIX __launch_bounds__(256) void colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda(float * RESTRICT const _data_force_a, float * RESTRICT const _data_force_b, float * RESTRICT  _data_pdfs_a, float * RESTRICT  _data_pdfs_b, float * RESTRICT const _data_phasefield, float * RESTRICT const _data_rho_a, float * RESTRICT const _data_rho_b, float * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_0, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_0, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_0, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_0, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_0, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_0, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_0, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_0, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, float beta, float omega_shear_a, float omega_shear_b, float sigma)
{
   if (blockDim.x*blockIdx.x + threadIdx.x + 1 < _size_force_a_0 - 1 && blockDim.y*blockIdx.y + threadIdx.y + 1 < _size_force_a_1 - 1 && blockDim.z*blockIdx.z + threadIdx.z + 1 < _size_force_a_2 - 1)
   {
      const int64_t ctr_0 = blockDim.x*blockIdx.x + threadIdx.x + 1;
      const int64_t ctr_1 = blockDim.y*blockIdx.y + threadIdx.y + 1;
      const int64_t ctr_2 = blockDim.z*blockIdx.z + threadIdx.z + 1;
      const float xi_2 = ((1.0f) / (omega_shear_a + omega_shear_b));
      const float xi_3 = omega_shear_a*omega_shear_b;
      const float xi_4 = xi_2*xi_3;
      const float xi_5 = xi_4*2.0f;
      const float xi_6 = omega_shear_a*-4.0f;
      const float xi_7 = xi_4*8.0f;
      const float xi_8 = xi_6 + xi_7;
      const float xi_9 = -xi_8;
      const float xi_12 = omega_shear_b*-4.0f;
      const float xi_13 = xi_12 + xi_7;
      const float xi_213 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_214 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_215 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_216 = -0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_217 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_218 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_219 = -xi_215 + xi_216 + xi_217 + xi_218;
      const float xi_220 = -0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_221 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_222 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_223 = 0.013888888888888888f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_224 = xi_220 - xi_221 + xi_222 + xi_223;
      const float xi_225 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_226 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_227 = -xi_225 + xi_226;
      const float xi_228 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_229 = -xi_228;
      const float xi_230 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_231 = xi_229 + xi_230;
      const float xi_232 = -0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_233 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_234 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_235 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_236 = -xi_234 + xi_235;
      const float xi_237 = xi_232 + xi_233 + xi_236;
      const float xi_238 = xi_231 + xi_237;
      const float xi_239 = -xi_213 + xi_214 + xi_219 + xi_224 + xi_227 + xi_238;
      const float xi_240 = xi_239*xi_239;
      const float xi_241 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_242 = -xi_241;
      const float xi_243 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_244 = xi_229 + xi_230 + xi_242 + xi_243;
      const float xi_245 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_246 = 0.055555555555555552f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_247 = -xi_245 + xi_246;
      const float xi_248 = xi_213 - xi_214;
      const float xi_249 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_250 = -xi_249;
      const float xi_251 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_252 = xi_234 - xi_235 + xi_250 + xi_251;
      const float xi_253 = xi_219 + xi_221 - xi_223 + xi_244 + xi_247 + xi_248 + xi_252;
      const float xi_254 = xi_253*xi_253;
      const float xi_255 = xi_245 - xi_246;
      const float xi_256 = xi_242 + xi_243;
      const float xi_257 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_258 = -xi_257;
      const float xi_259 = 0.22222222222222221f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_260 = xi_225 - xi_226 + xi_258 + xi_259;
      const float xi_261 = xi_215 + xi_216 + xi_217 - xi_218 + xi_224 + xi_248 + xi_255 + xi_256 + xi_260;
      const float xi_262 = xi_261*xi_261;
      const float xi_263 = xi_240 + xi_254 + xi_262;
      const float xi_264 = powf(xi_263, 0.5f);
      const float xi_265 = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const float xi_266 = sigma*xi_264*((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.5f) ? (omega_shear_a): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] < -0.5f) ? (omega_shear_b): ((_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] > 0.0f) ? (xi_265*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2]): (xi_13*xi_265 + xi_13*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2] + xi_5))));
      const bool xi_267 = xi_264 > 0.0f;
      const float xi_268 = ((xi_267) ? (xi_266*0.25f): (0.0f));
      const float xi_280 = ((1.0f) / (xi_263));
      const float xi_281 = xi_266*1.125f;
      const float xi_282 = ((xi_267) ? (xi_281*(xi_254*xi_280*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_286 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552f*(xi_253*xi_253) - 0.018518518518518517f)): (0.0f));
      const float xi_294 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552f*(xi_261*xi_261) - 0.018518518518518517f)): (0.0f));
      const float xi_297 = ((xi_267) ? (xi_281*(xi_262*xi_280*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_303 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552f*(xi_239*xi_239) - 0.018518518518518517f)): (0.0f));
      const float xi_305 = ((xi_267) ? (xi_281*(xi_240*xi_280*0.055555555555555552f - 0.018518518518518517f)): (0.0f));
      const float xi_313 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_314 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_315 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_316 = xi_220 + xi_222;
      const float xi_317 = xi_260 + xi_316;
      const float xi_318 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_319 = xi_249 - xi_251 + xi_318;
      const float xi_320 = xi_228 - xi_230 + xi_236 - xi_313 + xi_314 + xi_315 + xi_317 + xi_319;
      const float xi_321 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_320*xi_320) - 0.037037037037037035f)): (0.0f));
      const float xi_325 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_326 = 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2];
      const float xi_327 = -0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.027777777777777776f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_328 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2];
      const float xi_329 = xi_231 + xi_252 + xi_317 - xi_325 + xi_326 + xi_327 + xi_328;
      const float xi_330 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_329*xi_329) - 0.037037037037037035f)): (0.0f));
      const float xi_331 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_329*xi_329) - 0.037037037037037035f)): (0.0f));
      const float xi_333 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_320*xi_320) - 0.037037037037037035f)): (0.0f));
      const float xi_336 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_337 = xi_227 + xi_232 + xi_233 + xi_316;
      const float xi_338 = xi_241 - xi_243 + xi_325 - xi_326;
      const float xi_339 = xi_255 + xi_319 + xi_336 + xi_337 + xi_338;
      const float xi_340 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_339*xi_339) - 0.037037037037037035f)): (0.0f));
      const float xi_346 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_347 = xi_247 + xi_313 - xi_314;
      const float xi_348 = xi_250 + xi_251 + xi_256 + xi_327 + xi_337 + xi_346 + xi_347;
      const float xi_349 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_348*xi_348) - 0.037037037037037035f)): (0.0f));
      const float xi_354 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_355 = xi_237 + xi_244 + xi_255 + xi_258 + xi_259 + xi_318 + xi_327 + xi_354;
      const float xi_356 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_355*xi_355) - 0.037037037037037035f)): (0.0f));
      const float xi_360 = -0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2] + 0.1111111111111111f*_data_phasefield[_stride_phasefield_0*ctr_0 - _stride_phasefield_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2];
      const float xi_361 = xi_238 + xi_257 - xi_259 + xi_338 + xi_347 + xi_360;
      const float xi_362 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_361*xi_361) - 0.037037037037037035f)): (0.0f));
      const float xi_363 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_348*xi_348) - 0.037037037037037035f)): (0.0f));
      const float xi_365 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_339*xi_339) - 0.037037037037037035f)): (0.0f));
      const float xi_366 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_361*xi_361) - 0.037037037037037035f)): (0.0f));
      const float xi_368 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776f*(xi_355*xi_355) - 0.037037037037037035f)): (0.0f));
      const float xia_3_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3];
      const float xia_4_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_0 = xia_4_collide > 0.5f;
      const bool xi_1 = xia_4_collide < -0.5f;
      const float xi_10 = xia_4_collide*xia_4_collide;
      const bool xi_11 = xia_4_collide > 0.0f;
      const float xia_5_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3];
      const float xia_6_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3];
      const float xia_7_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3];
      const float xia_8_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3];
      const float xia_9_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3];
      const float xia_10_collide = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xia_11_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3];
      const float xia_12_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xia_13_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3];
      const float xi_276 = -xia_13_collide;
      const float xi_277 = xi_276 + xia_5_collide;
      const float xia_14_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3];
      const float xia_15_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3];
      const float xia_16_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3];
      const float xia_17_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xia_18_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3];
      const float xia_19_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3];
      const float xia_20_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3];
      const float xi_51 = xia_20_collide*0.16666666666666666f;
      const float xi_77 = xia_20_collide*0.083333333333333329f;
      const float xia_21_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3];
      const float xi_271 = -xia_21_collide;
      const float xi_272 = xi_271 + xia_9_collide;
      const float xia_22_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2];
      const float xia_23_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3];
      const float xi_26 = xia_23_collide*0.16666666666666666f;
      const float xi_54 = xia_23_collide*0.083333333333333329f;
      const float xi_58 = xia_23_collide*0.25f;
      const float xia_24_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3];
      const float xi_290 = -xia_24_collide;
      const float xi_291 = xi_290 + xia_6_collide;
      const float xia_25_collide = _data_force_a[_stride_force_a_0*ctr_0 + _stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2];
      const float xi_42 = xia_25_collide*0.16666666666666666f;
      const float xi_67 = xia_25_collide*0.083333333333333329f;
      const float xia_26_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3];
      const float xi_205 = xia_16_collide + xia_26_collide;
      const float xi_206 = xi_205 + xia_21_collide + xia_9_collide;
      const float xia_27_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3];
      const float xi_198 = xia_18_collide + xia_27_collide;
      const float xi_199 = xi_198 + xia_13_collide + xia_5_collide;
      const float xia_28_collide = _data_pdfs_a[_stride_pdfs_a_0*ctr_0 + _stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3];
      const float xi_202 = xia_19_collide + xia_28_collide;
      const float xi_203 = xi_202 + xia_24_collide + xia_6_collide;
      const float xia_29_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xia_0_collide = ((1.0f) / (-0.25f*((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_4_collide): (xi_10*xi_13 + xi_13*xia_4_collide + xi_5)))) + 2.0f));
      const float xia_1_collide = xi_2;
      const float xi_14 = xi_3*xia_1_collide;
      const float xi_15 = xi_14*2.0f;
      const float xi_16 = xi_14*8.0f;
      const float xi_17 = xi_16 + xi_6;
      const float xi_18 = xi_12 + xi_16;
      const float xi_19 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_17 + xi_15 - xi_17*xia_4_collide): (xi_10*xi_18 + xi_15 + xi_18*xia_4_collide))));
      const float rr_0_a_collide = 0.0f;
      const float xi_27 = rr_0_a_collide*0.25f;
      const float rr_1_a_collide = xi_19;
      const float xi_23 = rr_1_a_collide*0.5f;
      const float xi_62 = rr_1_a_collide*0.125f;
      const float xi_63 = xi_62*xia_23_collide;
      const float xi_71 = rr_1_a_collide*0.041666666666666664f;
      const float xi_307 = rr_1_a_collide*0.25f;
      const float rr_2_a_collide = xia_0_collide*(xi_19*-2.0f + 4.0f);
      const float xi_28 = rr_2_a_collide*xi_26 + xi_26 - xi_27*xia_23_collide;
      const float xi_43 = rr_2_a_collide*xi_42 - xi_27*xia_25_collide + xi_42;
      const float xi_52 = rr_2_a_collide*xi_51 - xi_27*xia_20_collide + xi_51;
      const float xi_55 = rr_2_a_collide*0.041666666666666664f;
      const float xi_56 = xi_55*xia_23_collide;
      const float xi_57 = xi_54 - xi_56;
      const float xi_68 = xi_55*xia_25_collide;
      const float xi_69 = -xi_67 + xi_68;
      const float xi_73 = -xi_54 + xi_56;
      const float xi_75 = xi_67 - xi_68;
      const float xi_78 = xi_55*xia_20_collide;
      const float xi_79 = xi_77 - xi_78;
      const float xi_98 = -xi_77 + xi_78;
      const float xi_274 = rr_2_a_collide*0.5f;
      const float xi_310 = rr_2_a_collide*0.25f;
      const float rho_a_collide = xia_10_collide;
      const float xi_195 = rho_a_collide*-0.1111111111111111f;
      const float xi_208 = rho_a_collide*-0.33333333333333331f;
      const float xi_210 = xi_199 + xi_208;
      const float xi_269 = rho_a_collide*0.33333333333333331f;
      const float u_0_a_collide = xia_17_collide;
      const float xi_20 = u_0_a_collide*xia_25_collide;
      const float xi_33 = xi_20*0.16666666666666666f;
      const float xi_34 = xi_20*0.083333333333333329f;
      const float xi_35 = rr_1_a_collide*xi_34;
      const float xi_36 = -xi_33 + xi_35;
      const float xi_39 = xi_33 - xi_35;
      const float xi_41 = rr_1_a_collide*xi_33;
      const float xi_59 = -u_0_a_collide*xi_58;
      const float xi_64 = u_0_a_collide*xi_63;
      const float xi_80 = -xi_34;
      const float xi_81 = xi_20*xi_71;
      const float xi_90 = u_0_a_collide*xia_20_collide;
      const float xi_91 = xi_90*-0.25f;
      const float xi_94 = xi_62*xi_90;
      const float xi_194 = u_0_a_collide*u_0_a_collide;
      const float xi_201 = rho_a_collide*xi_194*-0.33333333333333331f;
      const float xi_209 = rr_1_a_collide*(rho_a_collide*xi_194 - xi_203 - xi_206 - xi_208 - xia_14_collide - xia_15_collide);
      const float xi_287 = u_0_a_collide*xi_269;
      const float xi_288 = xi_272 + xi_287 + xia_16_collide - xia_26_collide;
      const float xi_289 = xi_274*xi_288;
      const float xi_292 = xi_287 + xi_291 + xia_19_collide - xia_28_collide;
      const float xi_293 = xi_274*xi_292;
      const float xi_311 = xi_288*xi_310;
      const float xi_312 = -xi_311;
      const float xi_352 = xi_292*xi_310;
      const float xi_353 = -xi_352;
      const float u_1_a_collide = xia_29_collide;
      const float xi_21 = u_1_a_collide*xia_23_collide;
      const float xi_24 = xi_21*0.16666666666666666f;
      const float xi_25 = rr_1_a_collide*xi_24;
      const float xi_44 = xi_21*0.083333333333333329f;
      const float xi_45 = rr_1_a_collide*xi_44;
      const float xi_46 = xi_24 - xi_45;
      const float xi_48 = -xi_24 + xi_45;
      const float xi_50 = xi_36 + xi_48;
      const float xi_53 = xi_39 + xi_46;
      const float xi_60 = u_1_a_collide*0.25f;
      const float xi_61 = -xi_60*xia_25_collide;
      const float xi_65 = u_1_a_collide*xi_62;
      const float xi_66 = xi_65*xia_25_collide;
      const float xi_70 = xi_59 + xi_61 + xi_64 + xi_66 + xi_69;
      const float xi_76 = xi_59 + xi_61 + xi_64 + xi_66 + xi_75;
      const float xi_83 = xi_60*xia_20_collide;
      const float xi_85 = xi_65*xia_20_collide;
      const float xi_196 = u_1_a_collide*u_1_a_collide;
      const float xi_197 = rho_a_collide*xi_196*-0.33333333333333331f + xi_195;
      const float xi_211 = rr_1_a_collide*(rho_a_collide*xi_196 - xi_206 - xi_210 - xia_11_collide - xia_7_collide);
      const float xi_270 = u_1_a_collide*xi_269;
      const float xi_273 = xi_270 + xi_272 - xia_16_collide + xia_26_collide;
      const float xi_275 = xi_273*xi_274;
      const float xi_278 = xi_270 + xi_277 - xia_18_collide + xia_27_collide;
      const float xi_279 = xi_274*xi_278;
      const float xi_306 = rho_a_collide*u_1_a_collide;
      const float xi_308 = xi_307*(u_0_a_collide*xi_306 + xi_205 + xi_271 - xia_9_collide);
      const float xi_309 = -xi_308;
      const float xi_323 = xi_273*xi_310;
      const float xi_335 = xi_278*xi_310;
      const float xi_345 = -xi_335;
      const float u_2_a_collide = xia_12_collide;
      const float xi_22 = u_2_a_collide*xia_20_collide;
      const float xi_29 = xi_22*0.16666666666666666f;
      const float xi_30 = xi_22*0.083333333333333329f;
      const float xi_31 = rr_1_a_collide*xi_30;
      const float xi_32 = -xi_29 + xi_31;
      const float xi_37 = xi_32 + xi_36;
      const float xi_38 = xi_29 - xi_31;
      const float xi_40 = xi_38 + xi_39;
      const float xi_47 = xi_38 + xi_46;
      const float xi_49 = rr_1_a_collide*xi_29;
      const float xi_72 = xi_22*xi_71 - xi_30 + xi_53;
      const float xi_74 = rr_1_a_collide*u_2_a_collide*xia_20_collide*-0.041666666666666664f + xi_30 + xi_50;
      const float xi_82 = xi_47 + xi_57 + xi_80 + xi_81;
      const float xi_84 = u_2_a_collide*xi_58;
      const float xi_86 = u_2_a_collide*xi_63;
      const float xi_87 = xi_83 + xi_84 - xi_85 - xi_86;
      const float xi_88 = -xi_83 - xi_84 + xi_85 + xi_86;
      const float xi_89 = xi_47 + xi_73 + xi_80 + xi_81;
      const float xi_92 = u_2_a_collide*xia_25_collide;
      const float xi_93 = xi_92*-0.25f;
      const float xi_95 = xi_62*xi_92;
      const float xi_96 = xi_69 + xi_91 + xi_93 + xi_94 + xi_95;
      const float xi_97 = xi_21*xi_71 + xi_40 - xi_44;
      const float xi_99 = rr_1_a_collide*u_1_a_collide*xia_23_collide*-0.041666666666666664f + xi_37 + xi_44;
      const float xi_100 = xi_75 + xi_91 + xi_93 + xi_94 + xi_95;
      const float xi_192 = u_2_a_collide*u_2_a_collide;
      const float xi_193 = rho_a_collide*xi_192*-0.33333333333333331f;
      const float xi_200 = rr_1_a_collide*(rho_a_collide*xi_194*-0.16666666666666666f - xi_193 - xi_197 - xi_199);
      const float xi_204 = rr_1_a_collide*(rho_a_collide*xi_196*-0.16666666666666666f - xi_193 - xi_195 - xi_201 - xi_203);
      const float xi_207 = rr_1_a_collide*(rho_a_collide*xi_192*-0.16666666666666666f - xi_197 - xi_201 - xi_206);
      const float xi_212 = rr_1_a_collide*(rho_a_collide*xi_192 - xi_203 - xi_210 - xia_3_collide - xia_8_collide);
      const float xi_283 = xi_200*-0.5f;
      const float xi_284 = xi_207*-0.5f;
      const float xi_285 = xi_211*0.5f + xi_283 + xi_284;
      const float xi_295 = xi_204*-0.5f;
      const float xi_296 = xi_209*0.5f + xi_284 + xi_295;
      const float xi_298 = u_2_a_collide*xi_269;
      const float xi_299 = xi_277 + xi_298 + xia_18_collide - xia_27_collide;
      const float xi_300 = xi_274*xi_299;
      const float xi_301 = xi_291 + xi_298 - xia_19_collide + xia_28_collide;
      const float xi_302 = xi_274*xi_301;
      const float xi_304 = xi_212*0.5f + xi_283 + xi_295;
      const float xi_322 = xi_207*0.25f;
      const float xi_324 = xi_322 + xi_323;
      const float xi_332 = xi_322 - xi_323;
      const float xi_334 = xi_307*(u_2_a_collide*xi_306 + xi_198 + xi_276 - xia_5_collide);
      const float xi_341 = xi_200*0.25f;
      const float xi_342 = xi_299*xi_310;
      const float xi_343 = xi_341 + xi_342;
      const float xi_344 = -xi_334;
      const float xi_350 = xi_307*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_202 + xi_290 - xia_6_collide);
      const float xi_351 = -xi_350;
      const float xi_357 = xi_204*0.25f;
      const float xi_358 = xi_301*xi_310;
      const float xi_359 = xi_357 + xi_358;
      const float xi_364 = xi_341 - xi_342;
      const float xi_367 = xi_357 - xi_358;
      const float forceTerm_0_a_collide = xi_20*xi_23 - xi_20 + xi_21*xi_23 - xi_21 + xi_22*xi_23 - xi_22;
      const float forceTerm_1_a_collide = xi_21*0.33333333333333331f - xi_25 + xi_28 + xi_37;
      const float forceTerm_2_a_collide = u_1_a_collide*xia_23_collide*0.33333333333333331f - xi_25 - xi_28 - xi_40;
      const float forceTerm_3_a_collide = u_0_a_collide*xia_25_collide*0.33333333333333331f - xi_41 - xi_43 - xi_47;
      const float forceTerm_4_a_collide = xi_20*0.33333333333333331f + xi_32 - xi_41 + xi_43 + xi_48;
      const float forceTerm_5_a_collide = xi_22*0.33333333333333331f - xi_49 + xi_50 + xi_52;
      const float forceTerm_6_a_collide = u_2_a_collide*xia_20_collide*0.33333333333333331f - xi_49 - xi_52 - xi_53;
      const float forceTerm_7_a_collide = xi_57 + xi_70 + xi_72;
      const float forceTerm_8_a_collide = -xi_70 - xi_73 - xi_74;
      const float forceTerm_9_a_collide = -xi_57 - xi_74 - xi_76;
      const float forceTerm_10_a_collide = xi_72 + xi_73 + xi_76;
      const float forceTerm_11_a_collide = xi_79 + xi_82 + xi_87;
      const float forceTerm_12_a_collide = xi_79 + xi_88 + xi_89;
      const float forceTerm_13_a_collide = xi_79 + xi_96 + xi_97;
      const float forceTerm_14_a_collide = -xi_96 - xi_98 - xi_99;
      const float forceTerm_15_a_collide = xi_82 + xi_88 + xi_98;
      const float forceTerm_16_a_collide = xi_87 + xi_89 + xi_98;
      const float forceTerm_17_a_collide = -xi_100 - xi_79 - xi_99;
      const float forceTerm_18_a_collide = xi_100 + xi_97 + xi_98;
      const float xib_3_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3];
      const float xi_407 = -xib_3_collide;
      const float xib_4_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2];
      const float xi_132 = xib_4_collide*0.16666666666666666f;
      const float xi_157 = xib_4_collide*0.083333333333333329f;
      const float xib_5_collide = _data_phasefield[_stride_phasefield_0*ctr_0 + _stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2];
      const bool xi_101 = xib_5_collide > 0.5f;
      const bool xi_102 = xib_5_collide < -0.5f;
      const float xi_103 = xib_5_collide*xib_5_collide;
      const bool xi_104 = xib_5_collide > 0.0f;
      const float xib_6_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3];
      const float xib_7_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3];
      const float xi_141 = xib_7_collide*0.16666666666666666f;
      const float xi_167 = xib_7_collide*0.083333333333333329f;
      const float xib_8_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3];
      const float xib_9_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3];
      const float xib_10_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3];
      const float xib_11_collide = _data_force_b[_stride_force_b_0*ctr_0 + _stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3];
      const float xi_117 = xib_11_collide*0.16666666666666666f;
      const float xi_144 = xib_11_collide*0.083333333333333329f;
      const float xi_148 = xib_11_collide*0.25f;
      const float xib_12_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3];
      const float xib_13_collide = _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xib_14_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3];
      const float xib_15_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3];
      const float xib_16_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2];
      const float xib_17_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2];
      const float xib_18_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3];
      const float xi_397 = -xib_18_collide;
      const float xib_19_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3];
      const float xi_398 = xi_397 + xib_19_collide;
      const float xib_20_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3];
      const float xib_21_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3];
      const float xi_392 = -xib_21_collide;
      const float xi_393 = xi_392 + xib_20_collide;
      const float xib_22_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3];
      const float xi_382 = xib_14_collide + xib_22_collide;
      const float xi_383 = xi_382 + xib_20_collide + xib_21_collide;
      const float xib_23_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3];
      const float xi_408 = xi_407 + xib_23_collide;
      const float xib_24_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3];
      const float xib_25_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3];
      const float xib_26_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3];
      const float xi_379 = xib_15_collide + xib_26_collide;
      const float xi_380 = xi_379 + xib_23_collide + xib_3_collide;
      const float xib_27_collide = _data_velocity[_stride_velocity_0*ctr_0 + _stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3];
      const float xib_28_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3];
      const float xib_29_collide = _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3];
      const float xi_375 = xib_29_collide + xib_6_collide;
      const float xi_376 = xi_375 + xib_18_collide + xib_19_collide;
      const float xib_0_collide = ((1.0f) / (-0.25f*((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_8 + xi_5 + xi_9*xib_5_collide): (xi_103*xi_13 + xi_13*xib_5_collide + xi_5)))) + 2.0f));
      const float xib_1_collide = xi_2;
      const float xi_105 = xi_3*xib_1_collide;
      const float xi_106 = xi_105*2.0f;
      const float xi_107 = xi_105*8.0f;
      const float xi_108 = xi_107 + xi_6;
      const float xi_109 = xi_107 + xi_12;
      const float xi_110 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_108 + xi_106 - xi_108*xib_5_collide): (xi_103*xi_109 + xi_106 + xi_109*xib_5_collide))));
      const float rr_0_b_collide = 0.0f;
      const float xi_118 = rr_0_b_collide*0.25f;
      const float rr_1_b_collide = xi_110;
      const float xi_114 = rr_1_b_collide*0.5f;
      const float xi_152 = rr_1_b_collide*0.125f;
      const float xi_153 = xi_152*xib_11_collide;
      const float xi_161 = rr_1_b_collide*0.041666666666666664f;
      const float xi_420 = rr_1_b_collide*0.25f;
      const float rr_2_b_collide = xib_0_collide*(xi_110*-2.0f + 4.0f);
      const float xi_119 = rr_2_b_collide*xi_117 + xi_117 - xi_118*xib_11_collide;
      const float xi_133 = rr_2_b_collide*xi_132 - xi_118*xib_4_collide + xi_132;
      const float xi_142 = rr_2_b_collide*xi_141 - xi_118*xib_7_collide + xi_141;
      const float xi_145 = rr_2_b_collide*0.041666666666666664f;
      const float xi_146 = xi_145*xib_11_collide;
      const float xi_147 = xi_144 - xi_146;
      const float xi_158 = xi_145*xib_4_collide;
      const float xi_159 = -xi_157 + xi_158;
      const float xi_163 = -xi_144 + xi_146;
      const float xi_165 = xi_157 - xi_158;
      const float xi_168 = xi_145*xib_7_collide;
      const float xi_169 = xi_167 - xi_168;
      const float xi_180 = -xi_167 + xi_168;
      const float xi_395 = rr_2_b_collide*0.5f;
      const float xi_423 = rr_2_b_collide*0.25f;
      const float rho_b_collide = xib_13_collide;
      const float xi_369 = rho_b_collide*-0.1111111111111111f;
      const float xi_385 = rho_b_collide*-0.33333333333333331f;
      const float xi_386 = xi_380 + xi_385;
      const float xi_390 = rho_b_collide*0.33333333333333331f;
      const float u_0_b_collide = xib_17_collide;
      const float xi_111 = u_0_b_collide*xib_4_collide;
      const float xi_124 = xi_111*0.16666666666666666f;
      const float xi_125 = xi_111*0.083333333333333329f;
      const float xi_126 = rr_1_b_collide*xi_125;
      const float xi_127 = -xi_124 + xi_126;
      const float xi_130 = xi_124 - xi_126;
      const float xi_131 = rr_1_b_collide*xi_124;
      const float xi_149 = -u_0_b_collide*xi_148;
      const float xi_154 = u_0_b_collide*xi_153;
      const float xi_170 = -xi_125;
      const float xi_171 = xi_111*xi_161;
      const float xi_183 = u_0_b_collide*xib_7_collide;
      const float xi_184 = xi_183*0.25f;
      const float xi_187 = xi_152*xi_183;
      const float xi_374 = u_0_b_collide*u_0_b_collide;
      const float xi_378 = rho_b_collide*xi_374*-0.33333333333333331f + xi_369;
      const float xi_387 = rr_1_b_collide*(rho_b_collide*xi_374 - xi_383 - xi_386 - xib_10_collide - xib_25_collide);
      const float xi_404 = u_0_b_collide*xi_390;
      const float xi_405 = xi_393 + xi_404 - xib_14_collide + xib_22_collide;
      const float xi_406 = xi_395*xi_405;
      const float xi_409 = xi_404 + xi_408 - xib_15_collide + xib_26_collide;
      const float xi_410 = xi_395*xi_409;
      const float xi_424 = xi_405*xi_423;
      const float xi_425 = -xi_424;
      const float xi_439 = xi_409*xi_423;
      const float xi_440 = -xi_439;
      const float u_1_b_collide = xib_27_collide;
      const float xi_112 = u_1_b_collide*xib_11_collide;
      const float xi_115 = xi_112*0.16666666666666666f;
      const float xi_116 = rr_1_b_collide*xi_115;
      const float xi_134 = xi_112*0.083333333333333329f;
      const float xi_135 = rr_1_b_collide*xi_134;
      const float xi_136 = xi_115 - xi_135;
      const float xi_138 = -xi_115 + xi_135;
      const float xi_140 = xi_127 + xi_138;
      const float xi_143 = xi_130 + xi_136;
      const float xi_150 = u_1_b_collide*0.25f;
      const float xi_151 = -xi_150*xib_4_collide;
      const float xi_155 = u_1_b_collide*xi_152;
      const float xi_156 = xi_155*xib_4_collide;
      const float xi_160 = xi_149 + xi_151 + xi_154 + xi_156 + xi_159;
      const float xi_166 = xi_149 + xi_151 + xi_154 + xi_156 + xi_165;
      const float xi_173 = xi_150*xib_7_collide;
      const float xi_175 = xi_155*xib_7_collide;
      const float xi_181 = rr_1_b_collide*u_1_b_collide*xib_11_collide*-0.041666666666666664f;
      const float xi_370 = u_1_b_collide*u_1_b_collide;
      const float xi_371 = rho_b_collide*xi_370*-0.33333333333333331f;
      const float xi_388 = rr_1_b_collide*(rho_b_collide*xi_370 - xi_376 - xi_383 - xi_385 - xib_24_collide - xib_28_collide);
      const float xi_391 = u_1_b_collide*xi_390;
      const float xi_394 = xi_391 + xi_393 + xib_14_collide - xib_22_collide;
      const float xi_396 = xi_394*xi_395;
      const float xi_399 = xi_391 + xi_398 - xib_29_collide + xib_6_collide;
      const float xi_400 = xi_395*xi_399;
      const float xi_419 = rho_b_collide*u_1_b_collide;
      const float xi_421 = xi_420*(u_0_b_collide*xi_419 + xi_382 + xi_392 - xib_20_collide);
      const float xi_422 = -xi_421;
      const float xi_427 = xi_394*xi_423;
      const float xi_431 = xi_399*xi_423;
      const float xi_436 = -xi_431;
      const float u_2_b_collide = xib_8_collide;
      const float xi_113 = u_2_b_collide*xib_7_collide;
      const float xi_120 = xi_113*0.16666666666666666f;
      const float xi_121 = xi_113*0.083333333333333329f;
      const float xi_122 = rr_1_b_collide*xi_121;
      const float xi_123 = -xi_120 + xi_122;
      const float xi_128 = xi_123 + xi_127;
      const float xi_129 = xi_120 - xi_122;
      const float xi_137 = xi_129 + xi_136;
      const float xi_139 = rr_1_b_collide*xi_120;
      const float xi_162 = xi_113*xi_161 - xi_121 + xi_143;
      const float xi_164 = rr_1_b_collide*u_2_b_collide*xib_7_collide*-0.041666666666666664f + xi_121 + xi_140;
      const float xi_172 = xi_137 + xi_147 + xi_170 + xi_171;
      const float xi_174 = u_2_b_collide*xi_148;
      const float xi_176 = u_2_b_collide*xi_153;
      const float xi_177 = xi_173 + xi_174 - xi_175 - xi_176;
      const float xi_178 = -xi_173 - xi_174 + xi_175 + xi_176;
      const float xi_179 = xi_137 + xi_163 + xi_170 + xi_171;
      const float xi_182 = xi_128 + xi_134 + xi_165 + xi_181;
      const float xi_185 = u_2_b_collide*xib_4_collide;
      const float xi_186 = xi_185*0.25f;
      const float xi_188 = xi_152*xi_185;
      const float xi_189 = xi_184 + xi_186 - xi_187 - xi_188;
      const float xi_190 = -xi_184 - xi_186 + xi_187 + xi_188;
      const float xi_191 = xi_128 + xi_134 + xi_159 + xi_181;
      const float xi_372 = u_2_b_collide*u_2_b_collide;
      const float xi_373 = rho_b_collide*xi_372*-0.33333333333333331f;
      const float xi_377 = rr_1_b_collide*(rho_b_collide*xi_374*-0.16666666666666666f - xi_369 - xi_371 - xi_373 - xi_376);
      const float xi_381 = rr_1_b_collide*(rho_b_collide*xi_370*-0.16666666666666666f - xi_373 - xi_378 - xi_380);
      const float xi_384 = rr_1_b_collide*(rho_b_collide*xi_372*-0.16666666666666666f - xi_371 - xi_378 - xi_383);
      const float xi_389 = rr_1_b_collide*(rho_b_collide*xi_372 - xi_376 - xi_386 - xib_12_collide - xib_9_collide);
      const float xi_401 = xi_384*-0.5f;
      const float xi_402 = xi_377*-0.5f;
      const float xi_403 = xi_388*0.5f + xi_401 + xi_402;
      const float xi_411 = xi_381*-0.5f;
      const float xi_412 = xi_387*0.5f + xi_401 + xi_411;
      const float xi_413 = u_2_b_collide*xi_390;
      const float xi_414 = xi_408 + xi_413 + xib_15_collide - xib_26_collide;
      const float xi_415 = xi_395*xi_414;
      const float xi_416 = xi_398 + xi_413 + xib_29_collide - xib_6_collide;
      const float xi_417 = xi_395*xi_416;
      const float xi_418 = xi_389*0.5f + xi_402 + xi_411;
      const float xi_426 = xi_384*0.25f;
      const float xi_428 = xi_426 + xi_427;
      const float xi_429 = xi_426 - xi_427;
      const float xi_430 = xi_420*(u_2_b_collide*xi_419 + xi_375 + xi_397 - xib_19_collide);
      const float xi_432 = xi_377*0.25f;
      const float xi_433 = xi_416*xi_423;
      const float xi_434 = xi_432 + xi_433;
      const float xi_435 = -xi_430;
      const float xi_437 = xi_420*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_379 + xi_407 - xib_23_collide);
      const float xi_438 = -xi_437;
      const float xi_441 = xi_381*0.25f;
      const float xi_442 = xi_414*xi_423;
      const float xi_443 = xi_441 + xi_442;
      const float xi_444 = xi_432 - xi_433;
      const float xi_445 = xi_441 - xi_442;
      const float forceTerm_0_b_collide = xi_111*xi_114 - xi_111 + xi_112*xi_114 - xi_112 + xi_113*xi_114 - xi_113;
      const float forceTerm_1_b_collide = xi_112*0.33333333333333331f - xi_116 + xi_119 + xi_128;
      const float forceTerm_2_b_collide = u_1_b_collide*xib_11_collide*0.33333333333333331f - xi_116 - xi_119 - xi_129 - xi_130;
      const float forceTerm_3_b_collide = u_0_b_collide*xib_4_collide*0.33333333333333331f - xi_131 - xi_133 - xi_137;
      const float forceTerm_4_b_collide = xi_111*0.33333333333333331f + xi_123 - xi_131 + xi_133 + xi_138;
      const float forceTerm_5_b_collide = xi_113*0.33333333333333331f - xi_139 + xi_140 + xi_142;
      const float forceTerm_6_b_collide = u_2_b_collide*xib_7_collide*0.33333333333333331f - xi_139 - xi_142 - xi_143;
      const float forceTerm_7_b_collide = xi_147 + xi_160 + xi_162;
      const float forceTerm_8_b_collide = -xi_160 - xi_163 - xi_164;
      const float forceTerm_9_b_collide = -xi_147 - xi_164 - xi_166;
      const float forceTerm_10_b_collide = xi_162 + xi_163 + xi_166;
      const float forceTerm_11_b_collide = xi_169 + xi_172 + xi_177;
      const float forceTerm_12_b_collide = xi_169 + xi_178 + xi_179;
      const float forceTerm_13_b_collide = -xi_180 - xi_182 - xi_189;
      const float forceTerm_14_b_collide = -xi_180 - xi_190 - xi_191;
      const float forceTerm_15_b_collide = xi_172 + xi_178 + xi_180;
      const float forceTerm_16_b_collide = xi_177 + xi_179 + xi_180;
      const float forceTerm_17_b_collide = -xi_169 - xi_182 - xi_190;
      const float forceTerm_18_b_collide = -xi_169 - xi_189 - xi_191;
      const float tmp_a0 = forceTerm_0_a_collide + xi_200 + xi_204 + xi_207 - xi_209 - xi_211 - xi_212 + xi_268 + xia_22_collide;
      const float tmp_a1 = forceTerm_1_a_collide - xi_275 - xi_279 + xi_282 + xi_285 + xia_7_collide;
      const float tmp_a2 = forceTerm_2_a_collide + xi_275 + xi_279 + xi_285 + xi_286 + xia_11_collide;
      const float tmp_a3 = forceTerm_3_a_collide + xi_289 + xi_293 + xi_294 + xi_296 + xia_14_collide;
      const float tmp_a4 = forceTerm_4_a_collide - xi_289 - xi_293 + xi_296 + xi_297 + xia_15_collide;
      const float tmp_a5 = forceTerm_5_a_collide - xi_300 - xi_302 + xi_303 + xi_304 + xia_8_collide;
      const float tmp_a6 = forceTerm_6_a_collide + xi_300 + xi_302 + xi_304 + xi_305 + xia_3_collide;
      const float tmp_a7 = forceTerm_7_a_collide + xi_309 + xi_312 + xi_321 + xi_324 + xia_16_collide;
      const float tmp_a8 = forceTerm_8_a_collide + xi_308 + xi_311 + xi_324 + xi_330 + xia_21_collide;
      const float tmp_a9 = forceTerm_9_a_collide + xi_308 + xi_312 + xi_331 + xi_332 + xia_9_collide;
      const float tmp_a10 = forceTerm_10_a_collide + xi_309 + xi_311 + xi_332 + xi_333 + xia_26_collide;
      const float tmp_a11 = forceTerm_11_a_collide + xi_334 + xi_335 + xi_340 + xi_343 + xia_13_collide;
      const float tmp_a12 = forceTerm_12_a_collide + xi_343 + xi_344 + xi_345 + xi_349 + xia_27_collide;
      const float tmp_a13 = forceTerm_13_a_collide + xi_351 + xi_353 + xi_356 + xi_359 + xia_19_collide;
      const float tmp_a14 = forceTerm_14_a_collide + xi_350 + xi_352 + xi_359 + xi_362 + xia_24_collide;
      const float tmp_a15 = forceTerm_15_a_collide + xi_335 + xi_344 + xi_363 + xi_364 + xia_18_collide;
      const float tmp_a16 = forceTerm_16_a_collide + xi_334 + xi_345 + xi_364 + xi_365 + xia_5_collide;
      const float tmp_a17 = forceTerm_17_a_collide + xi_350 + xi_353 + xi_366 + xi_367 + xia_6_collide;
      const float tmp_a18 = forceTerm_18_a_collide + xi_351 + xi_352 + xi_367 + xi_368 + xia_28_collide;
      const float tmp_b0 = forceTerm_0_b_collide + xi_268 + xi_377 + xi_381 + xi_384 - xi_387 - xi_388 - xi_389 + xib_16_collide;
      const float tmp_b1 = forceTerm_1_b_collide + xi_282 - xi_396 - xi_400 + xi_403 + xib_24_collide;
      const float tmp_b2 = forceTerm_2_b_collide + xi_286 + xi_396 + xi_400 + xi_403 + xib_28_collide;
      const float tmp_b3 = forceTerm_3_b_collide + xi_294 + xi_406 + xi_410 + xi_412 + xib_10_collide;
      const float tmp_b4 = forceTerm_4_b_collide + xi_297 - xi_406 - xi_410 + xi_412 + xib_25_collide;
      const float tmp_b5 = forceTerm_5_b_collide + xi_303 - xi_415 - xi_417 + xi_418 + xib_9_collide;
      const float tmp_b6 = forceTerm_6_b_collide + xi_305 + xi_415 + xi_417 + xi_418 + xib_12_collide;
      const float tmp_b7 = forceTerm_7_b_collide + xi_321 + xi_422 + xi_425 + xi_428 + xib_22_collide;
      const float tmp_b8 = forceTerm_8_b_collide + xi_330 + xi_421 + xi_424 + xi_428 + xib_21_collide;
      const float tmp_b9 = forceTerm_9_b_collide + xi_331 + xi_421 + xi_425 + xi_429 + xib_20_collide;
      const float tmp_b10 = forceTerm_10_b_collide + xi_333 + xi_422 + xi_424 + xi_429 + xib_14_collide;
      const float tmp_b11 = forceTerm_11_b_collide + xi_340 + xi_430 + xi_431 + xi_434 + xib_18_collide;
      const float tmp_b12 = forceTerm_12_b_collide + xi_349 + xi_434 + xi_435 + xi_436 + xib_6_collide;
      const float tmp_b13 = forceTerm_13_b_collide + xi_356 + xi_438 + xi_440 + xi_443 + xib_26_collide;
      const float tmp_b14 = forceTerm_14_b_collide + xi_362 + xi_437 + xi_439 + xi_443 + xib_3_collide;
      const float tmp_b15 = forceTerm_15_b_collide + xi_363 + xi_431 + xi_435 + xi_444 + xib_29_collide;
      const float tmp_b16 = forceTerm_16_b_collide + xi_365 + xi_430 + xi_436 + xi_444 + xib_19_collide;
      const float tmp_b17 = forceTerm_17_b_collide + xi_366 + xi_437 + xi_440 + xi_445 + xib_23_collide;
      const float tmp_b18 = forceTerm_18_b_collide + xi_368 + xi_438 + xi_439 + xi_445 + xib_15_collide;
      const float xirecolor_0 = tmp_a0 + tmp_b0;
      const float xirecolor_1 = _data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + _data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_2 = ((1.0f) / (xirecolor_1));
      const float xi_467 = xirecolor_2*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2];
      const float xirecolor_4 = tmp_a1 + tmp_b1;
      const float xirecolor_5 = xi_223;
      const float xirecolor_6 = xi_221;
      const float xirecolor_7 = xi_241;
      const float xi_446 = -xirecolor_7;
      const float xirecolor_8 = xi_446;
      const float xirecolor_9 = xi_228;
      const float xirecolor_10 = -xirecolor_9;
      const float xirecolor_11 = xi_230;
      const float xi_447 = xirecolor_10 + xirecolor_11;
      const float xirecolor_12 = xi_243;
      const float xi_448 = xirecolor_12 + xirecolor_8;
      const float xirecolor_13 = xi_447 + xi_448;
      const float xirecolor_14 = xi_245;
      const float xirecolor_15 = xi_246;
      const float xi_449 = xirecolor_14 - xirecolor_15;
      const float xirecolor_16 = -xi_449;
      const float xirecolor_17 = xi_215;
      const float xirecolor_18 = xi_216;
      const float xirecolor_19 = xi_217;
      const float xi_450 = xirecolor_18 + xirecolor_19;
      const float xirecolor_20 = xi_218;
      const float xirecolor_21 = xi_450 - xirecolor_17 + xirecolor_20;
      const float xirecolor_22 = xi_214;
      const float xirecolor_23 = xi_213;
      const float xi_451 = xirecolor_22 - xirecolor_23;
      const float xirecolor_24 = -xi_451;
      const float xirecolor_25 = xi_249;
      const float xirecolor_26 = -xirecolor_25;
      const float xirecolor_27 = xi_235;
      const float xirecolor_28 = xi_234;
      const float xirecolor_29 = xi_251;
      const float xi_452 = xirecolor_26 + xirecolor_29;
      const float xirecolor_30 = xi_452 - xirecolor_27 + xirecolor_28;
      const float xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
      const float xirecolor_32 = xi_220;
      const float xirecolor_33 = xi_222;
      const float xi_453 = xirecolor_32 + xirecolor_33;
      const float xirecolor_34 = xi_453 + xirecolor_5 - xirecolor_6;
      const float xirecolor_35 = xi_225;
      const float xirecolor_36 = xi_226;
      const float xi_454 = xirecolor_35 - xirecolor_36;
      const float xirecolor_37 = -xi_454;
      const float xirecolor_38 = xi_447;
      const float xirecolor_39 = xi_232;
      const float xirecolor_40 = xi_233;
      const float xi_455 = xirecolor_39 + xirecolor_40;
      const float xirecolor_41 = xirecolor_27 - xirecolor_28;
      const float xirecolor_42 = xi_455 + xirecolor_41;
      const float xirecolor_43 = xirecolor_38 + xirecolor_42;
      const float xirecolor_44 = xi_451 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
      const float xirecolor_45 = xi_449;
      const float xirecolor_46 = xi_448;
      const float xirecolor_47 = xi_257;
      const float xirecolor_48 = -xirecolor_47;
      const float xirecolor_49 = xi_259;
      const float xi_456 = xirecolor_48 + xirecolor_49;
      const float xirecolor_50 = xi_454 + xi_456;
      const float xirecolor_51 = xi_450 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
      const float xirecolor_52 = powf(xirecolor_31*xirecolor_31 + xirecolor_44*xirecolor_44 + xirecolor_51*xirecolor_51, 0.5f);
      const float xirecolor_53 = ((1.0f) / (xirecolor_52));
      const float xi_457 = xirecolor_31*xirecolor_53;
      const float xi_458 = xirecolor_51*xirecolor_53;
      const float xi_459 = xirecolor_44*xirecolor_53;
      const bool xirecolor_54 = xirecolor_52 > 0.0f;
      const float xirecolor_55 = beta*((1.0f) / (xirecolor_1*xirecolor_1))*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2]*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2];
      const float xirecolor_56 = xirecolor_55*(0.055555555555555552f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.055555555555555552f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_457): (0.0f));
      const float xirecolor_58 = tmp_a2 + tmp_b2;
      const float xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_457): (0.0f));
      const float xirecolor_60 = tmp_a3 + tmp_b3;
      const float xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_458): (0.0f));
      const float xirecolor_62 = tmp_a4 + tmp_b4;
      const float xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_458): (0.0f));
      const float xirecolor_64 = tmp_a5 + tmp_b5;
      const float xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_459): (0.0f));
      const float xirecolor_66 = tmp_a6 + tmp_b6;
      const float xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_459): (0.0f));
      const float xirecolor_68 = tmp_a7 + tmp_b7;
      const float xirecolor_69 = xi_313;
      const float xirecolor_70 = xi_314;
      const float xirecolor_71 = xi_453;
      const float xirecolor_72 = xirecolor_50 + xirecolor_71;
      const float xirecolor_73 = xi_318;
      const float xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
      const float xirecolor_75 = xi_315 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
      const float xirecolor_76 = xirecolor_53*0.70710678118654757f;
      const float xi_460 = xirecolor_75*xirecolor_76;
      const float xirecolor_77 = xirecolor_55*(0.027777777777777776f*_data_rho_a[_stride_rho_a_0*ctr_0 + _stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2] + 0.027777777777777776f*_data_rho_b[_stride_rho_b_0*ctr_0 + _stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2]);
      const float xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_460): (0.0f));
      const float xirecolor_79 = tmp_a8 + tmp_b8;
      const float xirecolor_80 = xi_325;
      const float xirecolor_81 = xi_326;
      const float xi_461 = -xirecolor_80 + xirecolor_81;
      const float xirecolor_82 = xi_327;
      const float xirecolor_83 = xi_328 + xi_461 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
      const float xi_462 = xirecolor_76*xirecolor_83;
      const float xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_462): (0.0f));
      const float xirecolor_85 = tmp_a9 + tmp_b9;
      const float xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_462): (0.0f));
      const float xirecolor_87 = tmp_a10 + tmp_b10;
      const float xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_460): (0.0f));
      const float xirecolor_89 = tmp_a11 + tmp_b11;
      const float xirecolor_90 = xi_455 + xirecolor_37 + xirecolor_71;
      const float xirecolor_91 = -xi_446 - xi_461 - xirecolor_12;
      const float xirecolor_92 = xi_336 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
      const float xi_463 = xirecolor_76*xirecolor_92;
      const float xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_463): (0.0f));
      const float xirecolor_94 = tmp_a12 + tmp_b12;
      const float xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
      const float xirecolor_96 = xi_346 + xi_452 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
      const float xi_464 = xirecolor_76*xirecolor_96;
      const float xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_464): (0.0f));
      const float xirecolor_98 = tmp_a13 + tmp_b13;
      const float xirecolor_99 = xi_354 + xi_456 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
      const float xi_465 = xirecolor_76*xirecolor_99;
      const float xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_465): (0.0f));
      const float xirecolor_101 = tmp_a14 + tmp_b14;
      const float xirecolor_102 = xi_360 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
      const float xi_466 = xirecolor_102*xirecolor_76;
      const float xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_466): (0.0f));
      const float xirecolor_104 = tmp_a15 + tmp_b15;
      const float xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_464): (0.0f));
      const float xirecolor_106 = tmp_a16 + tmp_b16;
      const float xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_463): (0.0f));
      const float xirecolor_108 = tmp_a17 + tmp_b17;
      const float xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_466): (0.0f));
      const float xirecolor_110 = tmp_a18 + tmp_b18;
      const float xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_465): (0.0f));
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
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2] = xi_467*xirecolor_0;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3] = xi_467*xirecolor_4 - xirecolor_57;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3] = xi_467*xirecolor_58 - xirecolor_59;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3] = xi_467*xirecolor_60 - xirecolor_61;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3] = xi_467*xirecolor_62 - xirecolor_63;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3] = xi_467*xirecolor_64 - xirecolor_65;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3] = xi_467*xirecolor_66 - xirecolor_67;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3] = xi_467*xirecolor_68 - xirecolor_78;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3] = xi_467*xirecolor_79 - xirecolor_84;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3] = xi_467*xirecolor_85 - xirecolor_86;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3] = xi_467*xirecolor_87 - xirecolor_88;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3] = xi_467*xirecolor_89 - xirecolor_93;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3] = xi_467*xirecolor_94 - xirecolor_97;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3] = xi_467*xirecolor_98 - xirecolor_100;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3] = xi_467*xirecolor_101 - xirecolor_103;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3] = xi_467*xirecolor_104 - xirecolor_105;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3] = xi_467*xirecolor_106 - xirecolor_107;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3] = xi_467*xirecolor_108 - xirecolor_109;
      _data_pdfs_b[_stride_pdfs_b_0*ctr_0 + _stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3] = xi_467*xirecolor_110 - xirecolor_111;
   } 
}
}


void ColorGradientCollideSweepSinglePrecisionCUDA::run(IBlock * block, gpuStream_t stream)
{
   
    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
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
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepSinglePrecisionCUDA::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block, gpuStream_t stream)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto rho_a = block->getData< gpu::GPUField<float> >(rho_aID);
    auto pdfs_a = block->getData< gpu::GPUField<float> >(pdfs_aID);
    auto pdfs_b = block->getData< gpu::GPUField<float> >(pdfs_bID);
    auto phasefield = block->getData< gpu::GPUField<float> >(phasefieldID);
    auto rho_b = block->getData< gpu::GPUField<float> >(rho_bID);
    auto velocity = block->getData< gpu::GPUField<float> >(velocityID);
    auto force_b = block->getData< gpu::GPUField<float> >(force_bID);
    auto force_a = block->getData< gpu::GPUField<float> >(force_aID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    float * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    float * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    float * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    float * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    float * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    float * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
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
    internal_colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda::colorgradientcollidesweepsingleprecisioncuda_colorgradientcollidesweepsingleprecisioncuda<<<_grid, _block, 0, stream>>>(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_0, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_0, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_0, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_0, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_0, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_0, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_0, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_0, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
