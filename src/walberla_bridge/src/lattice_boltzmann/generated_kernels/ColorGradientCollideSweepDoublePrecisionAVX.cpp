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
//! \\file ColorGradientCollideSweepDoublePrecisionAVX.cpp
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34, sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit 3247aa7395049ca5bfb69d34d55e45db19fa439c


#include <cmath>

#include "core/DataTypes.h"
#include "core/Macros.h"
#include "ColorGradientCollideSweepDoublePrecisionAVX.h"


#include <immintrin.h>



#define FUNC_PREFIX

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


namespace internal_6ed8e43fbb9656b49ba1d2a7ec0e653c {
static FUNC_PREFIX void colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(double * RESTRICT const _data_force_a, double * RESTRICT const _data_force_b, double * RESTRICT  _data_pdfs_a, double * RESTRICT  _data_pdfs_b, double * RESTRICT const _data_phasefield, double * RESTRICT const _data_rho_a, double * RESTRICT const _data_rho_b, double * RESTRICT const _data_velocity, int64_t const _size_force_a_0, int64_t const _size_force_a_1, int64_t const _size_force_a_2, int64_t const _stride_force_a_1, int64_t const _stride_force_a_2, int64_t const _stride_force_a_3, int64_t const _stride_force_b_1, int64_t const _stride_force_b_2, int64_t const _stride_force_b_3, int64_t const _stride_pdfs_a_1, int64_t const _stride_pdfs_a_2, int64_t const _stride_pdfs_a_3, int64_t const _stride_pdfs_b_1, int64_t const _stride_pdfs_b_2, int64_t const _stride_pdfs_b_3, int64_t const _stride_phasefield_1, int64_t const _stride_phasefield_2, int64_t const _stride_rho_a_1, int64_t const _stride_rho_a_2, int64_t const _stride_rho_b_1, int64_t const _stride_rho_b_2, int64_t const _stride_velocity_1, int64_t const _stride_velocity_2, int64_t const _stride_velocity_3, double beta, double omega_shear_a, double omega_shear_b, double sigma)
{
#ifdef _OPENMP
   #pragma omp parallel
#endif
   {
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
      const double xia_1_collide = xi_2;
      const double xi_14 = xi_3*xia_1_collide;
      const double xi_15 = xi_14*2.0;
      const double xi_16 = xi_14*8.0;
      const double xi_17 = xi_16 + xi_6;
      const double xi_18 = xi_12 + xi_16;
      const double rr_0_a_collide = 0.0;
      const double xi_27 = rr_0_a_collide*0.25;
      const double xib_1_collide = xi_2;
      const double xi_105 = xi_3*xib_1_collide;
      const double xi_106 = xi_105*2.0;
      const double xi_107 = xi_105*8.0;
      const double xi_108 = xi_107 + xi_6;
      const double xi_109 = xi_107 + xi_12;
      const double rr_0_b_collide = 0.0;
      const double xi_118 = rr_0_b_collide*0.25;
#ifdef _OPENMP
      #pragma omp for schedule(static)
#endif
      for (int64_t ctr_2 = 1; ctr_2 < _size_force_a_2 - 1; ctr_2 += 1)
      {
         for (int64_t ctr_1 = 1; ctr_1 < _size_force_a_1 - 1; ctr_1 += 1)
         {
            {
               for (int64_t ctr_0 = 1; ctr_0 < (int64_t)((_size_force_a_0 - 2) / (4)) * (4) + 1; ctr_0 += 4)
               {
                  const __m256d xi_213 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_214 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_215 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_216 = _mm256_mul_pd(_mm256_set_pd(-0.013888888888888888,-0.013888888888888888,-0.013888888888888888,-0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_217 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_218 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_219 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_215,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_216),xi_217),xi_218);
                  const __m256d xi_220 = _mm256_mul_pd(_mm256_set_pd(-0.055555555555555552,-0.055555555555555552,-0.055555555555555552,-0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_221 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_222 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_223 = _mm256_mul_pd(_mm256_set_pd(0.013888888888888888,0.013888888888888888,0.013888888888888888,0.013888888888888888),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_224 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_221,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_220),xi_222),xi_223);
                  const __m256d xi_225 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_226 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_227 = _mm256_add_pd(_mm256_mul_pd(xi_225,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_226);
                  const __m256d xi_228 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_229 = _mm256_mul_pd(xi_228,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_230 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_231 = _mm256_add_pd(xi_229,xi_230);
                  const __m256d xi_232 = _mm256_mul_pd(_mm256_set_pd(-0.22222222222222221,-0.22222222222222221,-0.22222222222222221,-0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_233 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_234 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_235 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0]));
                  const __m256d xi_236 = _mm256_add_pd(_mm256_mul_pd(xi_234,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_235);
                  const __m256d xi_237 = _mm256_add_pd(_mm256_add_pd(xi_232,xi_233),xi_236);
                  const __m256d xi_238 = _mm256_add_pd(xi_231,xi_237);
                  const __m256d xi_239 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_213,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_214),xi_219),xi_224),xi_227),xi_238);
                  const __m256d xi_240 = _mm256_mul_pd(xi_239,xi_239);
                  const __m256d xi_241 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_242 = _mm256_mul_pd(xi_241,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_243 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_244 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_229,xi_230),xi_242),xi_243);
                  const __m256d xi_245 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_246 = _mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_247 = _mm256_add_pd(_mm256_mul_pd(xi_245,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_246);
                  const __m256d xi_248 = _mm256_add_pd(_mm256_mul_pd(xi_214,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_213);
                  const __m256d xi_249 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_250 = _mm256_mul_pd(xi_249,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_251 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_252 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_235,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_234),xi_250),xi_251);
                  const __m256d xi_253 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_223,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_219),xi_221),xi_244),xi_247),xi_248),xi_252);
                  const __m256d xi_254 = _mm256_mul_pd(xi_253,xi_253);
                  const __m256d xi_255 = _mm256_add_pd(_mm256_mul_pd(xi_246,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_245);
                  const __m256d xi_256 = _mm256_add_pd(xi_242,xi_243);
                  const __m256d xi_257 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1]));
                  const __m256d xi_258 = _mm256_mul_pd(xi_257,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_259 = _mm256_mul_pd(_mm256_set_pd(0.22222222222222221,0.22222222222222221,0.22222222222222221,0.22222222222222221),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1]));
                  const __m256d xi_260 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_226,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_225),xi_258),xi_259);
                  const __m256d xi_261 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_218,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_215),xi_216),xi_217),xi_224),xi_248),xi_255),xi_256),xi_260);
                  const __m256d xi_262 = _mm256_mul_pd(xi_261,xi_261);
                  const __m256d xi_263 = _mm256_add_pd(_mm256_add_pd(xi_240,xi_254),xi_262);
                  const __m256d xi_264 = _mm256_sqrt_pd(xi_263);
                  const __m256d xi_265 = _mm256_mul_pd(_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]));
                  const __m256d xi_266 = _mm256_mul_pd(_mm256_mul_pd(xi_264,_mm256_set_pd(sigma,sigma,sigma,sigma)),_mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_265,_mm256_set_pd(xi_13,xi_13,xi_13,xi_13)),_mm256_mul_pd(_mm256_set_pd(xi_13,xi_13,xi_13,xi_13),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_265,_mm256_set_pd(xi_8,xi_8,xi_8,xi_8)),_mm256_mul_pd(_mm256_set_pd(xi_9,xi_9,xi_9,xi_9),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),_mm256_cmp_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),_mm256_cmp_pd(_mm256_set_pd(-0.5,-0.5,-0.5,-0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NLE_UQ)),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),_mm256_cmp_pd(_mm256_set_pd(0.5,0.5,0.5,0.5),_mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]),_CMP_NGE_UQ)));
                  const __m256d xi_267 = _mm256_cmp_pd(xi_264,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xi_268 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_266,_mm256_set_pd(0.25,0.25,0.25,0.25)),xi_267);
                  const __m256d xi_280 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xi_263);
                  const __m256d xi_281 = _mm256_mul_pd(xi_266,_mm256_set_pd(1.125,1.125,1.125,1.125));
                  const __m256d xi_282 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_254,xi_280),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_286 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_253,xi_253)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_294 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_261,xi_261)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_297 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_262,xi_280),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_303 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_mul_pd(xi_239,xi_239)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_305 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_240,xi_280),_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552)),_mm256_set_pd(-0.018518518518518517,-0.018518518518518517,-0.018518518518518517,-0.018518518518518517))),xi_267);
                  const __m256d xi_313 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_314 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_315 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256d xi_316 = _mm256_add_pd(xi_220,xi_222);
                  const __m256d xi_317 = _mm256_add_pd(xi_260,xi_316);
                  const __m256d xi_318 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.027777777777777776,-0.027777777777777776,-0.027777777777777776,-0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_319 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_251,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_249),xi_318);
                  const __m256d xi_320 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_230,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_313,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_228),xi_236),xi_314),xi_315),xi_317),xi_319);
                  const __m256d xi_321 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_320,xi_320)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_325 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1]));
                  const __m256d xi_326 = _mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1]));
                  const __m256d xi_327 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.027777777777777776,-0.027777777777777776,-0.027777777777777776,-0.027777777777777776),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_328 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1])));
                  const __m256d xi_329 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_325,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_231),xi_252),xi_317),xi_326),xi_327),xi_328);
                  const __m256d xi_330 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_329,xi_329)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_331 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_329,xi_329)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_333 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_320,xi_320)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_336 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256d xi_337 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_227,xi_232),xi_233),xi_316);
                  const __m256d xi_338 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_243,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_326,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_241),xi_325);
                  const __m256d xi_339 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_255,xi_319),xi_336),xi_337),xi_338);
                  const __m256d xi_340 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_339,xi_339)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_346 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0])));
                  const __m256d xi_347 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_314,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_247),xi_313);
                  const __m256d xi_348 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_250,xi_251),xi_256),xi_327),xi_337),xi_346),xi_347);
                  const __m256d xi_349 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_348,xi_348)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_354 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1])));
                  const __m256d xi_355 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_237,xi_244),xi_255),xi_258),xi_259),xi_318),xi_327),xi_354);
                  const __m256d xi_356 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_355,xi_355)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_360 = _mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.1111111111111111,0.1111111111111111,0.1111111111111111,0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1])),_mm256_mul_pd(_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111),_mm256_loadu_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1])));
                  const __m256d xi_361 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_259,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_238),xi_257),xi_338),xi_347),xi_360);
                  const __m256d xi_362 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_361,xi_361)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_363 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_348,xi_348)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_365 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_339,xi_339)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_366 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_361,xi_361)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xi_368 = _mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_281,_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(xi_280,_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776)),_mm256_mul_pd(xi_355,xi_355)),_mm256_set_pd(-0.037037037037037035,-0.037037037037037035,-0.037037037037037035,-0.037037037037037035))),xi_267);
                  const __m256d xia_3_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_4_collide = _mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const __m256d xi_0 = _mm256_cmp_pd(xia_4_collide,_mm256_set_pd(0.5,0.5,0.5,0.5),_CMP_NLE_UQ);
                  const __m256d xi_1 = _mm256_cmp_pd(xia_4_collide,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5),_CMP_NGE_UQ);
                  const __m256d xi_10 = _mm256_mul_pd(xia_4_collide,xia_4_collide);
                  const __m256d xi_11 = _mm256_cmp_pd(xia_4_collide,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xia_5_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_6_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_7_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_8_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_9_collide = _mm256_load_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0]);
                  const __m256d xi_42 = _mm256_mul_pd(xia_9_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_67 = _mm256_mul_pd(xia_9_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xia_10_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_271 = _mm256_mul_pd(xia_10_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xia_11_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xia_12_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_290 = _mm256_mul_pd(xia_12_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xia_13_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_198 = _mm256_add_pd(xia_13_collide,xia_7_collide);
                  const __m256d xia_14_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_199 = _mm256_add_pd(_mm256_add_pd(xi_198,xia_10_collide),xia_14_collide);
                  const __m256d xi_272 = _mm256_add_pd(xi_271,xia_14_collide);
                  const __m256d xia_15_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_16_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_291 = _mm256_add_pd(xi_290,xia_16_collide);
                  const __m256d xia_17_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0]);
                  const __m256d xi_26 = _mm256_mul_pd(xia_17_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_54 = _mm256_mul_pd(xia_17_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_58 = _mm256_mul_pd(xia_17_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xia_18_collide = _mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]);
                  const __m256d xia_19_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xia_20_collide = _mm256_loadu_pd(& _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0]);
                  const __m256d xi_51 = _mm256_mul_pd(xia_20_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_77 = _mm256_mul_pd(xia_20_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xia_21_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xia_22_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_23_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_24_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_25_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xia_26_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_202 = _mm256_add_pd(xia_26_collide,xia_5_collide);
                  const __m256d xi_203 = _mm256_add_pd(_mm256_add_pd(xi_202,xia_12_collide),xia_16_collide);
                  const __m256d xia_27_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0]);
                  const __m256d xia_28_collide = _mm256_loadu_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_205 = _mm256_add_pd(xia_28_collide,xia_6_collide);
                  const __m256d xia_29_collide = _mm256_load_pd(& _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0]);
                  const __m256d xi_206 = _mm256_add_pd(_mm256_add_pd(xi_205,xia_15_collide),xia_29_collide);
                  const __m256d xi_276 = _mm256_mul_pd(xia_29_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_277 = _mm256_add_pd(xi_276,xia_15_collide);
                  const __m256d xia_0_collide = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-0.25,-0.25,-0.25,-0.25),_mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_10,_mm256_set_pd(xi_13,xi_13,xi_13,xi_13)),_mm256_mul_pd(xia_4_collide,_mm256_set_pd(xi_13,xi_13,xi_13,xi_13))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_10,_mm256_set_pd(xi_8,xi_8,xi_8,xi_8)),_mm256_mul_pd(xia_4_collide,_mm256_set_pd(xi_9,xi_9,xi_9,xi_9))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),xi_11),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_1),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_0)),_mm256_set_pd(2.0,2.0,2.0,2.0)));
                  const __m256d xi_19 = _mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_10,_mm256_set_pd(xi_18,xi_18,xi_18,xi_18)),_mm256_mul_pd(xia_4_collide,_mm256_set_pd(xi_18,xi_18,xi_18,xi_18))),_mm256_set_pd(xi_15,xi_15,xi_15,xi_15)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_10,_mm256_set_pd(xi_17,xi_17,xi_17,xi_17)),_mm256_mul_pd(_mm256_mul_pd(xia_4_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_17,xi_17,xi_17,xi_17))),_mm256_set_pd(xi_15,xi_15,xi_15,xi_15)),xi_11),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_1),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_0);
                  const __m256d rr_1_a_collide = xi_19;
                  const __m256d xi_23 = _mm256_mul_pd(rr_1_a_collide,_mm256_set_pd(0.5,0.5,0.5,0.5));
                  const __m256d xi_62 = _mm256_mul_pd(rr_1_a_collide,_mm256_set_pd(0.125,0.125,0.125,0.125));
                  const __m256d xi_63 = _mm256_mul_pd(xi_62,xia_17_collide);
                  const __m256d xi_71 = _mm256_mul_pd(rr_1_a_collide,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664));
                  const __m256d xi_307 = _mm256_mul_pd(rr_1_a_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d rr_2_a_collide = _mm256_mul_pd(xia_0_collide,_mm256_add_pd(_mm256_mul_pd(xi_19,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0)),_mm256_set_pd(4.0,4.0,4.0,4.0)));
                  const __m256d xi_28 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_a_collide,xi_26),_mm256_mul_pd(_mm256_mul_pd(xia_17_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_27,xi_27,xi_27,xi_27))),xi_26);
                  const __m256d xi_43 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_a_collide,xi_42),_mm256_mul_pd(_mm256_mul_pd(xia_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_27,xi_27,xi_27,xi_27))),xi_42);
                  const __m256d xi_52 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_a_collide,xi_51),_mm256_mul_pd(_mm256_mul_pd(xia_20_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_27,xi_27,xi_27,xi_27))),xi_51);
                  const __m256d xi_55 = _mm256_mul_pd(rr_2_a_collide,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664));
                  const __m256d xi_56 = _mm256_mul_pd(xi_55,xia_17_collide);
                  const __m256d xi_57 = _mm256_add_pd(_mm256_mul_pd(xi_56,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_54);
                  const __m256d xi_68 = _mm256_mul_pd(xi_55,xia_9_collide);
                  const __m256d xi_69 = _mm256_add_pd(_mm256_mul_pd(xi_67,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_68);
                  const __m256d xi_73 = _mm256_add_pd(_mm256_mul_pd(xi_54,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_56);
                  const __m256d xi_75 = _mm256_add_pd(_mm256_mul_pd(xi_68,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_67);
                  const __m256d xi_78 = _mm256_mul_pd(xi_55,xia_20_collide);
                  const __m256d xi_79 = _mm256_add_pd(_mm256_mul_pd(xi_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_77);
                  const __m256d xi_98 = _mm256_add_pd(_mm256_mul_pd(xi_77,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_78);
                  const __m256d xi_274 = _mm256_mul_pd(rr_2_a_collide,_mm256_set_pd(0.5,0.5,0.5,0.5));
                  const __m256d xi_310 = _mm256_mul_pd(rr_2_a_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d rho_a_collide = xia_18_collide;
                  const __m256d xi_195 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_208 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_210 = _mm256_add_pd(xi_199,xi_208);
                  const __m256d xi_269 = _mm256_mul_pd(rho_a_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_a_collide = xia_11_collide;
                  const __m256d xi_20 = _mm256_mul_pd(u_0_a_collide,xia_9_collide);
                  const __m256d xi_33 = _mm256_mul_pd(xi_20,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_34 = _mm256_mul_pd(xi_20,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_35 = _mm256_mul_pd(rr_1_a_collide,xi_34);
                  const __m256d xi_36 = _mm256_add_pd(_mm256_mul_pd(xi_33,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_35);
                  const __m256d xi_39 = _mm256_add_pd(_mm256_mul_pd(xi_35,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_33);
                  const __m256d xi_41 = _mm256_mul_pd(rr_1_a_collide,xi_33);
                  const __m256d xi_59 = _mm256_mul_pd(_mm256_mul_pd(u_0_a_collide,xi_58),_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_64 = _mm256_mul_pd(u_0_a_collide,xi_63);
                  const __m256d xi_80 = _mm256_mul_pd(xi_34,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_81 = _mm256_mul_pd(xi_20,xi_71);
                  const __m256d xi_90 = _mm256_mul_pd(u_0_a_collide,xia_20_collide);
                  const __m256d xi_91 = _mm256_mul_pd(xi_90,_mm256_set_pd(-0.25,-0.25,-0.25,-0.25));
                  const __m256d xi_94 = _mm256_mul_pd(xi_62,xi_90);
                  const __m256d xi_194 = _mm256_mul_pd(u_0_a_collide,u_0_a_collide);
                  const __m256d xi_201 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_194),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_209 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_203,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_206,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_208,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_25_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_194)));
                  const __m256d xi_287 = _mm256_mul_pd(u_0_a_collide,xi_269);
                  const __m256d xi_288 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_28_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_277),xi_287),xia_6_collide);
                  const __m256d xi_289 = _mm256_mul_pd(xi_274,xi_288);
                  const __m256d xi_292 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_5_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_287),xi_291),xia_26_collide);
                  const __m256d xi_293 = _mm256_mul_pd(xi_274,xi_292);
                  const __m256d xi_311 = _mm256_mul_pd(xi_288,xi_310);
                  const __m256d xi_312 = _mm256_mul_pd(xi_311,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_352 = _mm256_mul_pd(xi_292,xi_310);
                  const __m256d xi_353 = _mm256_mul_pd(xi_352,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_a_collide = xia_21_collide;
                  const __m256d xi_21 = _mm256_mul_pd(u_1_a_collide,xia_17_collide);
                  const __m256d xi_24 = _mm256_mul_pd(xi_21,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_25 = _mm256_mul_pd(rr_1_a_collide,xi_24);
                  const __m256d xi_44 = _mm256_mul_pd(xi_21,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_45 = _mm256_mul_pd(rr_1_a_collide,xi_44);
                  const __m256d xi_46 = _mm256_add_pd(_mm256_mul_pd(xi_45,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_24);
                  const __m256d xi_48 = _mm256_add_pd(_mm256_mul_pd(xi_24,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_45);
                  const __m256d xi_50 = _mm256_add_pd(xi_36,xi_48);
                  const __m256d xi_53 = _mm256_add_pd(xi_39,xi_46);
                  const __m256d xi_60 = _mm256_mul_pd(u_1_a_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_61 = _mm256_mul_pd(_mm256_mul_pd(xi_60,xia_9_collide),_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_65 = _mm256_mul_pd(u_1_a_collide,xi_62);
                  const __m256d xi_66 = _mm256_mul_pd(xi_65,xia_9_collide);
                  const __m256d xi_70 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_59,xi_61),xi_64),xi_66),xi_69);
                  const __m256d xi_76 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_59,xi_61),xi_64),xi_66),xi_75);
                  const __m256d xi_83 = _mm256_mul_pd(xi_60,xia_20_collide);
                  const __m256d xi_85 = _mm256_mul_pd(xi_65,xia_20_collide);
                  const __m256d xi_192 = _mm256_mul_pd(u_1_a_collide,u_1_a_collide);
                  const __m256d xi_193 = _mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_192),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_211 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_206,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_210,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_24_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_192)));
                  const __m256d xi_270 = _mm256_mul_pd(u_1_a_collide,xi_269);
                  const __m256d xi_273 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_7_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_270),xi_272),xia_13_collide);
                  const __m256d xi_275 = _mm256_mul_pd(xi_273,xi_274);
                  const __m256d xi_278 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_6_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_270),xi_277),xia_28_collide);
                  const __m256d xi_279 = _mm256_mul_pd(xi_274,xi_278);
                  const __m256d xi_306 = _mm256_mul_pd(rho_a_collide,u_1_a_collide);
                  const __m256d xi_308 = _mm256_mul_pd(xi_307,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_15_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_a_collide,xi_306)),xi_205),xi_276));
                  const __m256d xi_309 = _mm256_mul_pd(xi_308,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_323 = _mm256_mul_pd(xi_278,xi_310);
                  const __m256d xi_335 = _mm256_mul_pd(xi_273,xi_310);
                  const __m256d xi_345 = _mm256_mul_pd(xi_335,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_a_collide = xia_19_collide;
                  const __m256d xi_22 = _mm256_mul_pd(u_2_a_collide,xia_20_collide);
                  const __m256d xi_29 = _mm256_mul_pd(xi_22,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_30 = _mm256_mul_pd(xi_22,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_31 = _mm256_mul_pd(rr_1_a_collide,xi_30);
                  const __m256d xi_32 = _mm256_add_pd(_mm256_mul_pd(xi_29,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_31);
                  const __m256d xi_37 = _mm256_add_pd(xi_32,xi_36);
                  const __m256d xi_38 = _mm256_add_pd(_mm256_mul_pd(xi_31,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_29);
                  const __m256d xi_40 = _mm256_add_pd(xi_38,xi_39);
                  const __m256d xi_47 = _mm256_add_pd(xi_38,xi_46);
                  const __m256d xi_49 = _mm256_mul_pd(rr_1_a_collide,xi_29);
                  const __m256d xi_72 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_30,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_22,xi_71)),xi_53);
                  const __m256d xi_74 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rr_1_a_collide,u_2_a_collide),xia_20_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),xi_30),xi_50);
                  const __m256d xi_82 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_47,xi_57),xi_80),xi_81);
                  const __m256d xi_84 = _mm256_mul_pd(u_2_a_collide,xi_58);
                  const __m256d xi_86 = _mm256_mul_pd(u_2_a_collide,xi_63);
                  const __m256d xi_87 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_85,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_83),xi_84);
                  const __m256d xi_88 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_83,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_85),xi_86);
                  const __m256d xi_89 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_47,xi_73),xi_80),xi_81);
                  const __m256d xi_92 = _mm256_mul_pd(u_2_a_collide,xia_9_collide);
                  const __m256d xi_93 = _mm256_mul_pd(xi_92,_mm256_set_pd(-0.25,-0.25,-0.25,-0.25));
                  const __m256d xi_95 = _mm256_mul_pd(xi_62,xi_92);
                  const __m256d xi_96 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_69,xi_91),xi_93),xi_94),xi_95);
                  const __m256d xi_97 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_44,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_21,xi_71)),xi_40);
                  const __m256d xi_99 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rr_1_a_collide,u_1_a_collide),xia_17_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),xi_37),xi_44);
                  const __m256d xi_100 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_75,xi_91),xi_93),xi_94),xi_95);
                  const __m256d xi_196 = _mm256_mul_pd(u_2_a_collide,u_2_a_collide);
                  const __m256d xi_197 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_196),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_195);
                  const __m256d xi_200 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_193,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_197,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_199,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_194),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_204 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_197,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_201,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_203,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_192),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_207 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_193,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_195,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_201,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_206,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,xi_196),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_212 = _mm256_mul_pd(rr_1_a_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_203,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_210,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_23_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xia_3_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_a_collide,xi_196)));
                  const __m256d xi_283 = _mm256_mul_pd(xi_200,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_284 = _mm256_mul_pd(xi_207,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_285 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_211,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_283),xi_284);
                  const __m256d xi_295 = _mm256_mul_pd(xi_204,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_296 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_209,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_284),xi_295);
                  const __m256d xi_298 = _mm256_mul_pd(u_2_a_collide,xi_269);
                  const __m256d xi_299 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_13_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_272),xi_298),xia_7_collide);
                  const __m256d xi_300 = _mm256_mul_pd(xi_274,xi_299);
                  const __m256d xi_301 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_26_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_291),xi_298),xia_5_collide);
                  const __m256d xi_302 = _mm256_mul_pd(xi_274,xi_301);
                  const __m256d xi_304 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_212,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_283),xi_295);
                  const __m256d xi_322 = _mm256_mul_pd(xi_207,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_324 = _mm256_add_pd(xi_322,xi_323);
                  const __m256d xi_332 = _mm256_add_pd(_mm256_mul_pd(xi_323,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_322);
                  const __m256d xi_334 = _mm256_mul_pd(xi_307,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_14_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_a_collide,xi_306)),xi_198),xi_271));
                  const __m256d xi_341 = _mm256_mul_pd(xi_200,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_342 = _mm256_mul_pd(xi_299,xi_310);
                  const __m256d xi_343 = _mm256_add_pd(xi_341,xi_342);
                  const __m256d xi_344 = _mm256_mul_pd(xi_334,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_350 = _mm256_mul_pd(xi_307,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xia_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_a_collide,u_0_a_collide),u_2_a_collide)),xi_202),xi_290));
                  const __m256d xi_351 = _mm256_mul_pd(xi_350,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_357 = _mm256_mul_pd(xi_204,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_358 = _mm256_mul_pd(xi_301,xi_310);
                  const __m256d xi_359 = _mm256_add_pd(xi_357,xi_358);
                  const __m256d xi_364 = _mm256_add_pd(_mm256_mul_pd(xi_342,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_341);
                  const __m256d xi_367 = _mm256_add_pd(_mm256_mul_pd(xi_358,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_357);
                  const __m256d forceTerm_0_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_20,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_21,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_22,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_20,xi_23)),_mm256_mul_pd(xi_21,xi_23)),_mm256_mul_pd(xi_22,xi_23));
                  const __m256d forceTerm_1_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_21,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_25,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_28),xi_37);
                  const __m256d forceTerm_2_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_25,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_28,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_40,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_1_a_collide,xia_17_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_3_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_41,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_43,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_47,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_0_a_collide,xia_9_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_4_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_20,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_41,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_32),xi_43),xi_48);
                  const __m256d forceTerm_5_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_22,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_50),xi_52);
                  const __m256d forceTerm_6_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_52,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_53,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_2_a_collide,xia_20_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_7_a_collide = _mm256_add_pd(_mm256_add_pd(xi_57,xi_70),xi_72);
                  const __m256d forceTerm_8_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_70,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_73,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_9_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_74,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_76,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_10_a_collide = _mm256_add_pd(_mm256_add_pd(xi_72,xi_73),xi_76);
                  const __m256d forceTerm_11_a_collide = _mm256_add_pd(_mm256_add_pd(xi_79,xi_82),xi_87);
                  const __m256d forceTerm_12_a_collide = _mm256_add_pd(_mm256_add_pd(xi_79,xi_88),xi_89);
                  const __m256d forceTerm_13_a_collide = _mm256_add_pd(_mm256_add_pd(xi_79,xi_96),xi_97);
                  const __m256d forceTerm_14_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_96,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_98,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_99,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_15_a_collide = _mm256_add_pd(_mm256_add_pd(xi_82,xi_88),xi_98);
                  const __m256d forceTerm_16_a_collide = _mm256_add_pd(_mm256_add_pd(xi_87,xi_89),xi_98);
                  const __m256d forceTerm_17_a_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_100,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_79,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_99,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_18_a_collide = _mm256_add_pd(_mm256_add_pd(xi_100,xi_97),xi_98);
                  const __m256d xib_3_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_4_collide = _mm256_load_pd(& _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const __m256d xi_101 = _mm256_cmp_pd(xib_4_collide,_mm256_set_pd(0.5,0.5,0.5,0.5),_CMP_NLE_UQ);
                  const __m256d xi_102 = _mm256_cmp_pd(xib_4_collide,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5),_CMP_NGE_UQ);
                  const __m256d xi_103 = _mm256_mul_pd(xib_4_collide,xib_4_collide);
                  const __m256d xi_104 = _mm256_cmp_pd(xib_4_collide,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xib_5_collide = _mm256_load_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0]);
                  const __m256d xi_133 = _mm256_mul_pd(xib_5_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_158 = _mm256_mul_pd(xib_5_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xib_6_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_7_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_8_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_9_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_10_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_405 = _mm256_mul_pd(xib_10_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_11_collide = _mm256_load_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0]);
                  const __m256d xib_12_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_13_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_14_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_15_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_392 = _mm256_mul_pd(xib_15_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xib_16_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_17_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0]);
                  const __m256d xib_18_collide = _mm256_loadu_pd(& _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0]);
                  const __m256d xib_19_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_20_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_382 = _mm256_add_pd(xib_20_collide,xib_3_collide);
                  const __m256d xib_21_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0]);
                  const __m256d xi_142 = _mm256_mul_pd(xib_21_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_168 = _mm256_mul_pd(xib_21_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xib_22_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_406 = _mm256_add_pd(xi_405,xib_22_collide);
                  const __m256d xib_23_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_383 = _mm256_add_pd(_mm256_add_pd(xi_382,xib_23_collide),xib_8_collide);
                  const __m256d xi_397 = _mm256_mul_pd(xib_23_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_398 = _mm256_add_pd(xi_397,xib_8_collide);
                  const __m256d xib_24_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_375 = _mm256_add_pd(xib_16_collide,xib_24_collide);
                  const __m256d xib_25_collide = _mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const __m256d xib_26_collide = _mm256_load_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_376 = _mm256_add_pd(_mm256_add_pd(xi_375,xib_15_collide),xib_26_collide);
                  const __m256d xi_393 = _mm256_add_pd(xi_392,xib_26_collide);
                  const __m256d xib_27_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xib_28_collide = _mm256_loadu_pd(& _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0]);
                  const __m256d xi_117 = _mm256_mul_pd(xib_28_collide,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_145 = _mm256_mul_pd(xib_28_collide,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_149 = _mm256_mul_pd(xib_28_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xib_29_collide = _mm256_loadu_pd(& _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0]);
                  const __m256d xi_379 = _mm256_add_pd(xib_19_collide,xib_29_collide);
                  const __m256d xi_380 = _mm256_add_pd(_mm256_add_pd(xi_379,xib_10_collide),xib_22_collide);
                  const __m256d xib_0_collide = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(-0.25,-0.25,-0.25,-0.25),_mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(xi_13,xi_13,xi_13,xi_13)),_mm256_mul_pd(xib_4_collide,_mm256_set_pd(xi_13,xi_13,xi_13,xi_13))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(xi_8,xi_8,xi_8,xi_8)),_mm256_mul_pd(xib_4_collide,_mm256_set_pd(xi_9,xi_9,xi_9,xi_9))),_mm256_set_pd(xi_5,xi_5,xi_5,xi_5)),xi_104),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_102),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_101)),_mm256_set_pd(2.0,2.0,2.0,2.0)));
                  const __m256d xi_110 = _mm256_blendv_pd(_mm256_blendv_pd(_mm256_blendv_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(xi_109,xi_109,xi_109,xi_109)),_mm256_mul_pd(xib_4_collide,_mm256_set_pd(xi_109,xi_109,xi_109,xi_109))),_mm256_set_pd(xi_106,xi_106,xi_106,xi_106)),_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_103,_mm256_set_pd(xi_108,xi_108,xi_108,xi_108)),_mm256_mul_pd(_mm256_mul_pd(xib_4_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_108,xi_108,xi_108,xi_108))),_mm256_set_pd(xi_106,xi_106,xi_106,xi_106)),xi_104),_mm256_set_pd(omega_shear_b,omega_shear_b,omega_shear_b,omega_shear_b),xi_102),_mm256_set_pd(omega_shear_a,omega_shear_a,omega_shear_a,omega_shear_a),xi_101);
                  const __m256d rr_1_b_collide = xi_110;
                  const __m256d xi_114 = _mm256_mul_pd(rr_1_b_collide,_mm256_set_pd(0.5,0.5,0.5,0.5));
                  const __m256d xi_153 = _mm256_mul_pd(rr_1_b_collide,_mm256_set_pd(0.125,0.125,0.125,0.125));
                  const __m256d xi_154 = _mm256_mul_pd(xi_153,xib_28_collide);
                  const __m256d xi_162 = _mm256_mul_pd(rr_1_b_collide,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664));
                  const __m256d xi_420 = _mm256_mul_pd(rr_1_b_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d rr_2_b_collide = _mm256_mul_pd(xib_0_collide,_mm256_add_pd(_mm256_mul_pd(xi_110,_mm256_set_pd(-2.0,-2.0,-2.0,-2.0)),_mm256_set_pd(4.0,4.0,4.0,4.0)));
                  const __m256d xi_119 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_b_collide,xi_117),_mm256_mul_pd(_mm256_mul_pd(xib_28_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_118,xi_118,xi_118,xi_118))),xi_117);
                  const __m256d xi_134 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_b_collide,xi_133),_mm256_mul_pd(_mm256_mul_pd(xib_5_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_118,xi_118,xi_118,xi_118))),xi_133);
                  const __m256d xi_143 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(rr_2_b_collide,xi_142),_mm256_mul_pd(_mm256_mul_pd(xib_21_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_set_pd(xi_118,xi_118,xi_118,xi_118))),xi_142);
                  const __m256d xi_146 = _mm256_mul_pd(rr_2_b_collide,_mm256_set_pd(0.041666666666666664,0.041666666666666664,0.041666666666666664,0.041666666666666664));
                  const __m256d xi_147 = _mm256_mul_pd(xi_146,xib_28_collide);
                  const __m256d xi_148 = _mm256_add_pd(_mm256_mul_pd(xi_147,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_145);
                  const __m256d xi_159 = _mm256_mul_pd(xi_146,xib_5_collide);
                  const __m256d xi_160 = _mm256_add_pd(_mm256_mul_pd(xi_158,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_159);
                  const __m256d xi_164 = _mm256_add_pd(_mm256_mul_pd(xi_145,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_147);
                  const __m256d xi_166 = _mm256_add_pd(_mm256_mul_pd(xi_159,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_158);
                  const __m256d xi_169 = _mm256_mul_pd(xi_146,xib_21_collide);
                  const __m256d xi_170 = _mm256_add_pd(_mm256_mul_pd(xi_169,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_168);
                  const __m256d xi_189 = _mm256_add_pd(_mm256_mul_pd(xi_168,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_169);
                  const __m256d xi_395 = _mm256_mul_pd(rr_2_b_collide,_mm256_set_pd(0.5,0.5,0.5,0.5));
                  const __m256d xi_423 = _mm256_mul_pd(rr_2_b_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d rho_b_collide = xib_25_collide;
                  const __m256d xi_372 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.1111111111111111,-0.1111111111111111,-0.1111111111111111,-0.1111111111111111));
                  const __m256d xi_385 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_386 = _mm256_add_pd(xi_380,xi_385);
                  const __m256d xi_390 = _mm256_mul_pd(rho_b_collide,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331));
                  const __m256d u_0_b_collide = xib_11_collide;
                  const __m256d xi_111 = _mm256_mul_pd(u_0_b_collide,xib_5_collide);
                  const __m256d xi_124 = _mm256_mul_pd(xi_111,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_125 = _mm256_mul_pd(xi_111,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_126 = _mm256_mul_pd(rr_1_b_collide,xi_125);
                  const __m256d xi_127 = _mm256_add_pd(_mm256_mul_pd(xi_124,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_126);
                  const __m256d xi_130 = _mm256_add_pd(_mm256_mul_pd(xi_126,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_124);
                  const __m256d xi_132 = _mm256_mul_pd(rr_1_b_collide,xi_124);
                  const __m256d xi_150 = _mm256_mul_pd(_mm256_mul_pd(u_0_b_collide,xi_149),_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_155 = _mm256_mul_pd(u_0_b_collide,xi_154);
                  const __m256d xi_171 = _mm256_mul_pd(xi_125,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_172 = _mm256_mul_pd(xi_111,xi_162);
                  const __m256d xi_181 = _mm256_mul_pd(u_0_b_collide,xib_21_collide);
                  const __m256d xi_182 = _mm256_mul_pd(xi_181,_mm256_set_pd(-0.25,-0.25,-0.25,-0.25));
                  const __m256d xi_185 = _mm256_mul_pd(xi_153,xi_181);
                  const __m256d xi_371 = _mm256_mul_pd(u_0_b_collide,u_0_b_collide);
                  const __m256d xi_378 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_371),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_387 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_383,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_386,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_12_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_27_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_371)));
                  const __m256d xi_404 = _mm256_mul_pd(u_0_b_collide,xi_390);
                  const __m256d xi_407 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_29_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_404),xi_406),xib_19_collide);
                  const __m256d xi_408 = _mm256_mul_pd(xi_395,xi_407);
                  const __m256d xi_409 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_20_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_398),xi_404),xib_3_collide);
                  const __m256d xi_410 = _mm256_mul_pd(xi_395,xi_409);
                  const __m256d xi_424 = _mm256_mul_pd(xi_409,xi_423);
                  const __m256d xi_425 = _mm256_mul_pd(xi_424,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_439 = _mm256_mul_pd(xi_407,xi_423);
                  const __m256d xi_440 = _mm256_mul_pd(xi_439,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_1_b_collide = xib_18_collide;
                  const __m256d xi_112 = _mm256_mul_pd(u_1_b_collide,xib_28_collide);
                  const __m256d xi_115 = _mm256_mul_pd(xi_112,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_116 = _mm256_mul_pd(rr_1_b_collide,xi_115);
                  const __m256d xi_135 = _mm256_mul_pd(xi_112,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_136 = _mm256_mul_pd(rr_1_b_collide,xi_135);
                  const __m256d xi_137 = _mm256_add_pd(_mm256_mul_pd(xi_136,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_115);
                  const __m256d xi_139 = _mm256_add_pd(_mm256_mul_pd(xi_115,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_136);
                  const __m256d xi_141 = _mm256_add_pd(xi_127,xi_139);
                  const __m256d xi_144 = _mm256_add_pd(xi_130,xi_137);
                  const __m256d xi_151 = _mm256_mul_pd(u_1_b_collide,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_152 = _mm256_mul_pd(_mm256_mul_pd(xi_151,xib_5_collide),_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_156 = _mm256_mul_pd(u_1_b_collide,xi_153);
                  const __m256d xi_157 = _mm256_mul_pd(xi_156,xib_5_collide);
                  const __m256d xi_161 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_150,xi_152),xi_155),xi_157),xi_160);
                  const __m256d xi_167 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_150,xi_152),xi_155),xi_157),xi_166);
                  const __m256d xi_174 = _mm256_mul_pd(xi_151,xib_21_collide);
                  const __m256d xi_176 = _mm256_mul_pd(xi_156,xib_21_collide);
                  const __m256d xi_369 = _mm256_mul_pd(u_1_b_collide,u_1_b_collide);
                  const __m256d xi_370 = _mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_369),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331));
                  const __m256d xi_388 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_376,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_383,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_385,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_14_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_6_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_369)));
                  const __m256d xi_391 = _mm256_mul_pd(u_1_b_collide,xi_390);
                  const __m256d xi_394 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_24_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_391),xi_393),xib_16_collide);
                  const __m256d xi_396 = _mm256_mul_pd(xi_394,xi_395);
                  const __m256d xi_399 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_3_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_391),xi_398),xib_20_collide);
                  const __m256d xi_400 = _mm256_mul_pd(xi_395,xi_399);
                  const __m256d xi_419 = _mm256_mul_pd(rho_b_collide,u_1_b_collide);
                  const __m256d xi_421 = _mm256_mul_pd(xi_420,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_8_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_0_b_collide,xi_419)),xi_382),xi_397));
                  const __m256d xi_422 = _mm256_mul_pd(xi_421,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_427 = _mm256_mul_pd(xi_399,xi_423);
                  const __m256d xi_431 = _mm256_mul_pd(xi_394,xi_423);
                  const __m256d xi_436 = _mm256_mul_pd(xi_431,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d u_2_b_collide = xib_17_collide;
                  const __m256d xi_113 = _mm256_mul_pd(u_2_b_collide,xib_21_collide);
                  const __m256d xi_120 = _mm256_mul_pd(xi_113,_mm256_set_pd(0.16666666666666666,0.16666666666666666,0.16666666666666666,0.16666666666666666));
                  const __m256d xi_121 = _mm256_mul_pd(xi_113,_mm256_set_pd(0.083333333333333329,0.083333333333333329,0.083333333333333329,0.083333333333333329));
                  const __m256d xi_122 = _mm256_mul_pd(rr_1_b_collide,xi_121);
                  const __m256d xi_123 = _mm256_add_pd(_mm256_mul_pd(xi_120,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_122);
                  const __m256d xi_128 = _mm256_add_pd(xi_123,xi_127);
                  const __m256d xi_129 = _mm256_add_pd(_mm256_mul_pd(xi_122,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_120);
                  const __m256d xi_131 = _mm256_add_pd(xi_129,xi_130);
                  const __m256d xi_138 = _mm256_add_pd(xi_129,xi_137);
                  const __m256d xi_140 = _mm256_mul_pd(rr_1_b_collide,xi_120);
                  const __m256d xi_163 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_121,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_113,xi_162)),xi_144);
                  const __m256d xi_165 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rr_1_b_collide,u_2_b_collide),xib_21_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),xi_121),xi_141);
                  const __m256d xi_173 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_138,xi_148),xi_171),xi_172);
                  const __m256d xi_175 = _mm256_mul_pd(u_2_b_collide,xi_149);
                  const __m256d xi_177 = _mm256_mul_pd(u_2_b_collide,xi_154);
                  const __m256d xi_178 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_176,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_177,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_174),xi_175);
                  const __m256d xi_179 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_174,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_175,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_176),xi_177);
                  const __m256d xi_180 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_138,xi_164),xi_171),xi_172);
                  const __m256d xi_183 = _mm256_mul_pd(u_2_b_collide,xib_5_collide);
                  const __m256d xi_184 = _mm256_mul_pd(xi_183,_mm256_set_pd(-0.25,-0.25,-0.25,-0.25));
                  const __m256d xi_186 = _mm256_mul_pd(xi_153,xi_183);
                  const __m256d xi_187 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_160,xi_182),xi_184),xi_185),xi_186);
                  const __m256d xi_188 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_135,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_112,xi_162)),xi_131);
                  const __m256d xi_190 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(rr_1_b_collide,u_1_b_collide),xib_28_collide),_mm256_set_pd(-0.041666666666666664,-0.041666666666666664,-0.041666666666666664,-0.041666666666666664)),xi_128),xi_135);
                  const __m256d xi_191 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_166,xi_182),xi_184),xi_185),xi_186);
                  const __m256d xi_373 = _mm256_mul_pd(u_2_b_collide,u_2_b_collide);
                  const __m256d xi_374 = _mm256_add_pd(_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_373),_mm256_set_pd(-0.33333333333333331,-0.33333333333333331,-0.33333333333333331,-0.33333333333333331)),xi_372);
                  const __m256d xi_377 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_370,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_374,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_376,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_371),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_381 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_374,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_378,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_380,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_369),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_384 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_370,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_372,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_378,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_383,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,xi_373),_mm256_set_pd(-0.16666666666666666,-0.16666666666666666,-0.16666666666666666,-0.16666666666666666))));
                  const __m256d xi_389 = _mm256_mul_pd(rr_1_b_collide,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_376,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_386,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_7_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xib_9_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(rho_b_collide,xi_373)));
                  const __m256d xi_401 = _mm256_mul_pd(xi_377,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_402 = _mm256_mul_pd(xi_384,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_403 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_388,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_401),xi_402);
                  const __m256d xi_411 = _mm256_mul_pd(xi_381,_mm256_set_pd(-0.5,-0.5,-0.5,-0.5));
                  const __m256d xi_412 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_387,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_402),xi_411);
                  const __m256d xi_413 = _mm256_mul_pd(u_2_b_collide,xi_390);
                  const __m256d xi_414 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_19_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_406),xi_413),xib_29_collide);
                  const __m256d xi_415 = _mm256_mul_pd(xi_395,xi_414);
                  const __m256d xi_416 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_16_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_393),xi_413),xib_24_collide);
                  const __m256d xi_417 = _mm256_mul_pd(xi_395,xi_416);
                  const __m256d xi_418 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_389,_mm256_set_pd(0.5,0.5,0.5,0.5)),xi_401),xi_411);
                  const __m256d xi_426 = _mm256_mul_pd(xi_384,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_428 = _mm256_add_pd(xi_426,xi_427);
                  const __m256d xi_429 = _mm256_add_pd(_mm256_mul_pd(xi_427,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_426);
                  const __m256d xi_430 = _mm256_mul_pd(xi_420,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_26_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(u_2_b_collide,xi_419)),xi_375),xi_392));
                  const __m256d xi_432 = _mm256_mul_pd(xi_377,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_433 = _mm256_mul_pd(xi_416,xi_423);
                  const __m256d xi_434 = _mm256_add_pd(xi_432,xi_433);
                  const __m256d xi_435 = _mm256_mul_pd(xi_430,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_437 = _mm256_mul_pd(xi_420,_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xib_22_collide,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(_mm256_mul_pd(rho_b_collide,u_0_b_collide),u_2_b_collide)),xi_379),xi_405));
                  const __m256d xi_438 = _mm256_mul_pd(xi_437,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xi_441 = _mm256_mul_pd(xi_381,_mm256_set_pd(0.25,0.25,0.25,0.25));
                  const __m256d xi_442 = _mm256_mul_pd(xi_414,xi_423);
                  const __m256d xi_443 = _mm256_add_pd(xi_441,xi_442);
                  const __m256d xi_444 = _mm256_add_pd(_mm256_mul_pd(xi_433,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_432);
                  const __m256d xi_445 = _mm256_add_pd(_mm256_mul_pd(xi_442,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_441);
                  const __m256d forceTerm_0_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_111,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_112,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_113,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_111,xi_114)),_mm256_mul_pd(xi_112,xi_114)),_mm256_mul_pd(xi_113,xi_114));
                  const __m256d forceTerm_1_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_112,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_119),xi_128);
                  const __m256d forceTerm_2_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_116,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_119,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_131,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_1_b_collide,xib_28_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_3_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_134,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_138,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_0_b_collide,xib_5_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_4_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_111,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_132,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_123),xi_134),xi_139);
                  const __m256d forceTerm_5_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_113,_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)),_mm256_mul_pd(xi_140,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_141),xi_143);
                  const __m256d forceTerm_6_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_140,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_143,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_144,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(_mm256_mul_pd(u_2_b_collide,xib_21_collide),_mm256_set_pd(0.33333333333333331,0.33333333333333331,0.33333333333333331,0.33333333333333331)));
                  const __m256d forceTerm_7_b_collide = _mm256_add_pd(_mm256_add_pd(xi_148,xi_161),xi_163);
                  const __m256d forceTerm_8_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_161,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_164,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_9_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_148,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_165,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_167,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_10_b_collide = _mm256_add_pd(_mm256_add_pd(xi_163,xi_164),xi_167);
                  const __m256d forceTerm_11_b_collide = _mm256_add_pd(_mm256_add_pd(xi_170,xi_173),xi_178);
                  const __m256d forceTerm_12_b_collide = _mm256_add_pd(_mm256_add_pd(xi_170,xi_179),xi_180);
                  const __m256d forceTerm_13_b_collide = _mm256_add_pd(_mm256_add_pd(xi_170,xi_187),xi_188);
                  const __m256d forceTerm_14_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_187,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_189,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_190,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_15_b_collide = _mm256_add_pd(_mm256_add_pd(xi_173,xi_179),xi_189);
                  const __m256d forceTerm_16_b_collide = _mm256_add_pd(_mm256_add_pd(xi_178,xi_180),xi_189);
                  const __m256d forceTerm_17_b_collide = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_170,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_190,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_191,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d forceTerm_18_b_collide = _mm256_add_pd(_mm256_add_pd(xi_188,xi_189),xi_191);
                  const __m256d tmp_a0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_209,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_211,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_212,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_a_collide),xi_200),xi_204),xi_207),xi_268),xia_27_collide);
                  const __m256d tmp_a1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_275,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_279,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_a_collide),xi_282),xi_285),xia_8_collide);
                  const __m256d tmp_a2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_a_collide,xi_275),xi_279),xi_285),xi_286),xia_24_collide);
                  const __m256d tmp_a3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_a_collide,xi_289),xi_293),xi_294),xi_296),xia_25_collide);
                  const __m256d tmp_a4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_289,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_293,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_a_collide),xi_296),xi_297),xia_22_collide);
                  const __m256d tmp_a5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_300,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_302,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_a_collide),xi_303),xi_304),xia_23_collide);
                  const __m256d tmp_a6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_a_collide,xi_300),xi_302),xi_304),xi_305),xia_3_collide);
                  const __m256d tmp_a7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_a_collide,xi_309),xi_312),xi_321),xi_324),xia_6_collide);
                  const __m256d tmp_a8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_a_collide,xi_308),xi_311),xi_324),xi_330),xia_29_collide);
                  const __m256d tmp_a9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_a_collide,xi_308),xi_312),xi_331),xi_332),xia_15_collide);
                  const __m256d tmp_a10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_a_collide,xi_309),xi_311),xi_332),xi_333),xia_28_collide);
                  const __m256d tmp_a11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_a_collide,xi_334),xi_335),xi_340),xi_343),xia_10_collide);
                  const __m256d tmp_a12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_a_collide,xi_343),xi_344),xi_345),xi_349),xia_13_collide);
                  const __m256d tmp_a13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_a_collide,xi_351),xi_353),xi_356),xi_359),xia_26_collide);
                  const __m256d tmp_a14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_a_collide,xi_350),xi_352),xi_359),xi_362),xia_12_collide);
                  const __m256d tmp_a15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_a_collide,xi_335),xi_344),xi_363),xi_364),xia_7_collide);
                  const __m256d tmp_a16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_a_collide,xi_334),xi_345),xi_364),xi_365),xia_14_collide);
                  const __m256d tmp_a17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_a_collide,xi_350),xi_353),xi_366),xi_367),xia_16_collide);
                  const __m256d tmp_a18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_a_collide,xi_351),xi_352),xi_367),xi_368),xia_5_collide);
                  const __m256d tmp_b0 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_387,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_388,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xi_389,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_0_b_collide),xi_268),xi_377),xi_381),xi_384),xib_13_collide);
                  const __m256d tmp_b1 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_396,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_400,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_1_b_collide),xi_282),xi_403),xib_6_collide);
                  const __m256d tmp_b2 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_2_b_collide,xi_286),xi_396),xi_400),xi_403),xib_14_collide);
                  const __m256d tmp_b3 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_3_b_collide,xi_294),xi_408),xi_410),xi_412),xib_27_collide);
                  const __m256d tmp_b4 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_408,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_410,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_4_b_collide),xi_297),xi_412),xib_12_collide);
                  const __m256d tmp_b5 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_415,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_417,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),forceTerm_5_b_collide),xi_303),xi_418),xib_7_collide);
                  const __m256d tmp_b6 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_6_b_collide,xi_305),xi_415),xi_417),xi_418),xib_9_collide);
                  const __m256d tmp_b7 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_7_b_collide,xi_321),xi_422),xi_425),xi_428),xib_3_collide);
                  const __m256d tmp_b8 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_8_b_collide,xi_330),xi_421),xi_424),xi_428),xib_23_collide);
                  const __m256d tmp_b9 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_9_b_collide,xi_331),xi_421),xi_425),xi_429),xib_8_collide);
                  const __m256d tmp_b10 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_10_b_collide,xi_333),xi_422),xi_424),xi_429),xib_20_collide);
                  const __m256d tmp_b11 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_11_b_collide,xi_340),xi_430),xi_431),xi_434),xib_15_collide);
                  const __m256d tmp_b12 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_12_b_collide,xi_349),xi_434),xi_435),xi_436),xib_16_collide);
                  const __m256d tmp_b13 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_13_b_collide,xi_356),xi_438),xi_440),xi_443),xib_19_collide);
                  const __m256d tmp_b14 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_14_b_collide,xi_362),xi_437),xi_439),xi_443),xib_10_collide);
                  const __m256d tmp_b15 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_15_b_collide,xi_363),xi_431),xi_435),xi_444),xib_24_collide);
                  const __m256d tmp_b16 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_16_b_collide,xi_365),xi_430),xi_436),xi_444),xib_26_collide);
                  const __m256d tmp_b17 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_17_b_collide,xi_366),xi_437),xi_440),xi_445),xib_22_collide);
                  const __m256d tmp_b18 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(forceTerm_18_b_collide,xi_368),xi_438),xi_439),xi_445),xib_29_collide);
                  const __m256d xirecolor_0 = _mm256_add_pd(tmp_a0,tmp_b0);
                  const __m256d xirecolor_1 = _mm256_add_pd(_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_2 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_1);
                  const __m256d xi_467 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_3 = _mm256_mul_pd(xirecolor_2,_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_4 = _mm256_add_pd(tmp_a1,tmp_b1);
                  const __m256d xirecolor_5 = xi_223;
                  const __m256d xirecolor_6 = xi_221;
                  const __m256d xirecolor_7 = xi_241;
                  const __m256d xi_446 = _mm256_mul_pd(xirecolor_7,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_8 = xi_446;
                  const __m256d xirecolor_9 = xi_228;
                  const __m256d xirecolor_10 = _mm256_mul_pd(xirecolor_9,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_11 = xi_230;
                  const __m256d xi_447 = _mm256_add_pd(xirecolor_10,xirecolor_11);
                  const __m256d xirecolor_12 = xi_243;
                  const __m256d xi_448 = _mm256_add_pd(xirecolor_12,xirecolor_8);
                  const __m256d xirecolor_13 = _mm256_add_pd(xi_447,xi_448);
                  const __m256d xirecolor_14 = xi_245;
                  const __m256d xirecolor_15 = xi_246;
                  const __m256d xi_449 = _mm256_add_pd(_mm256_mul_pd(xirecolor_15,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_14);
                  const __m256d xirecolor_16 = _mm256_mul_pd(xi_449,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_17 = xi_215;
                  const __m256d xirecolor_18 = xi_216;
                  const __m256d xirecolor_19 = xi_217;
                  const __m256d xi_450 = _mm256_add_pd(xirecolor_18,xirecolor_19);
                  const __m256d xirecolor_20 = xi_218;
                  const __m256d xirecolor_21 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_17,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_450),xirecolor_20);
                  const __m256d xirecolor_22 = xi_214;
                  const __m256d xirecolor_23 = xi_213;
                  const __m256d xi_451 = _mm256_add_pd(_mm256_mul_pd(xirecolor_23,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_22);
                  const __m256d xirecolor_24 = _mm256_mul_pd(xi_451,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_25 = xi_249;
                  const __m256d xirecolor_26 = _mm256_mul_pd(xirecolor_25,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_27 = xi_235;
                  const __m256d xirecolor_28 = xi_234;
                  const __m256d xirecolor_29 = xi_251;
                  const __m256d xi_452 = _mm256_add_pd(xirecolor_26,xirecolor_29);
                  const __m256d xirecolor_30 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_27,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_452),xirecolor_28);
                  const __m256d xirecolor_31 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_5,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_13),xirecolor_16),xirecolor_21),xirecolor_24),xirecolor_30),xirecolor_6);
                  const __m256d xirecolor_32 = xi_220;
                  const __m256d xirecolor_33 = xi_222;
                  const __m256d xi_453 = _mm256_add_pd(xirecolor_32,xirecolor_33);
                  const __m256d xirecolor_34 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_6,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_453),xirecolor_5);
                  const __m256d xirecolor_35 = xi_225;
                  const __m256d xirecolor_36 = xi_226;
                  const __m256d xi_454 = _mm256_add_pd(_mm256_mul_pd(xirecolor_36,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_35);
                  const __m256d xirecolor_37 = _mm256_mul_pd(xi_454,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_38 = xi_447;
                  const __m256d xirecolor_39 = xi_232;
                  const __m256d xirecolor_40 = xi_233;
                  const __m256d xi_455 = _mm256_add_pd(xirecolor_39,xirecolor_40);
                  const __m256d xirecolor_41 = _mm256_add_pd(_mm256_mul_pd(xirecolor_28,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_27);
                  const __m256d xirecolor_42 = _mm256_add_pd(xi_455,xirecolor_41);
                  const __m256d xirecolor_43 = _mm256_add_pd(xirecolor_38,xirecolor_42);
                  const __m256d xirecolor_44 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_451,xirecolor_21),xirecolor_34),xirecolor_37),xirecolor_43);
                  const __m256d xirecolor_45 = xi_449;
                  const __m256d xirecolor_46 = xi_448;
                  const __m256d xirecolor_47 = xi_257;
                  const __m256d xirecolor_48 = _mm256_mul_pd(xirecolor_47,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0));
                  const __m256d xirecolor_49 = xi_259;
                  const __m256d xi_456 = _mm256_add_pd(xirecolor_48,xirecolor_49);
                  const __m256d xirecolor_50 = _mm256_add_pd(xi_454,xi_456);
                  const __m256d xirecolor_51 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_20,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_450),xirecolor_17),xirecolor_24),xirecolor_34),xirecolor_45),xirecolor_46),xirecolor_50);
                  const __m256d xirecolor_52 = _mm256_sqrt_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_31,xirecolor_31),_mm256_mul_pd(xirecolor_44,xirecolor_44)),_mm256_mul_pd(xirecolor_51,xirecolor_51)));
                  const __m256d xirecolor_53 = _mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),xirecolor_52);
                  const __m256d xi_457 = _mm256_mul_pd(xirecolor_31,xirecolor_53);
                  const __m256d xi_458 = _mm256_mul_pd(xirecolor_51,xirecolor_53);
                  const __m256d xi_459 = _mm256_mul_pd(xirecolor_44,xirecolor_53);
                  const __m256d xirecolor_54 = _mm256_cmp_pd(xirecolor_52,_mm256_set_pd(0.0,0.0,0.0,0.0),_CMP_NLE_UQ);
                  const __m256d xirecolor_55 = _mm256_mul_pd(_mm256_mul_pd(_mm256_mul_pd(_mm256_set_pd(beta,beta,beta,beta),_mm256_div_pd(_mm256_set_pd(1.0,1.0,1.0,1.0),_mm256_mul_pd(xirecolor_1,xirecolor_1))),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]));
                  const __m256d xirecolor_56 = _mm256_mul_pd(xirecolor_55,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.055555555555555552,0.055555555555555552,0.055555555555555552,0.055555555555555552),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_57 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_457,xirecolor_54));
                  const __m256d xirecolor_58 = _mm256_add_pd(tmp_a2,tmp_b2);
                  const __m256d xirecolor_59 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_457,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_60 = _mm256_add_pd(tmp_a3,tmp_b3);
                  const __m256d xirecolor_61 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_458,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_62 = _mm256_add_pd(tmp_a4,tmp_b4);
                  const __m256d xirecolor_63 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_458,xirecolor_54));
                  const __m256d xirecolor_64 = _mm256_add_pd(tmp_a5,tmp_b5);
                  const __m256d xirecolor_65 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_459,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_66 = _mm256_add_pd(tmp_a6,tmp_b6);
                  const __m256d xirecolor_67 = _mm256_mul_pd(xirecolor_56,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_459,xirecolor_54));
                  const __m256d xirecolor_68 = _mm256_add_pd(tmp_a7,tmp_b7);
                  const __m256d xirecolor_69 = xi_313;
                  const __m256d xirecolor_70 = xi_314;
                  const __m256d xirecolor_71 = xi_453;
                  const __m256d xirecolor_72 = _mm256_add_pd(xirecolor_50,xirecolor_71);
                  const __m256d xirecolor_73 = xi_318;
                  const __m256d xirecolor_74 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_29,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_25),xirecolor_73);
                  const __m256d xirecolor_75 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_11,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xirecolor_69,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),xi_315),xirecolor_41),xirecolor_70),xirecolor_72),xirecolor_74),xirecolor_9);
                  const __m256d xirecolor_76 = _mm256_mul_pd(xirecolor_53,_mm256_set_pd(0.70710678118654757,0.70710678118654757,0.70710678118654757,0.70710678118654757));
                  const __m256d xi_460 = _mm256_mul_pd(xirecolor_75,xirecolor_76);
                  const __m256d xirecolor_77 = _mm256_mul_pd(xirecolor_55,_mm256_add_pd(_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0])),_mm256_mul_pd(_mm256_set_pd(0.027777777777777776,0.027777777777777776,0.027777777777777776,0.027777777777777776),_mm256_load_pd(& _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]))));
                  const __m256d xirecolor_78 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_460,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_79 = _mm256_add_pd(tmp_a8,tmp_b8);
                  const __m256d xirecolor_80 = xi_325;
                  const __m256d xirecolor_81 = xi_326;
                  const __m256d xi_461 = _mm256_add_pd(_mm256_mul_pd(xirecolor_80,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_81);
                  const __m256d xirecolor_82 = xi_327;
                  const __m256d xirecolor_83 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_328,xi_461),xirecolor_30),xirecolor_38),xirecolor_72),xirecolor_82);
                  const __m256d xi_462 = _mm256_mul_pd(xirecolor_76,xirecolor_83);
                  const __m256d xirecolor_84 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_462,xirecolor_54));
                  const __m256d xirecolor_85 = _mm256_add_pd(tmp_a9,tmp_b9);
                  const __m256d xirecolor_86 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_462,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_87 = _mm256_add_pd(tmp_a10,tmp_b10);
                  const __m256d xirecolor_88 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_460,xirecolor_54));
                  const __m256d xirecolor_89 = _mm256_add_pd(tmp_a11,tmp_b11);
                  const __m256d xirecolor_90 = _mm256_add_pd(_mm256_add_pd(xi_455,xirecolor_37),xirecolor_71);
                  const __m256d xirecolor_91 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xi_446,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_461,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0))),_mm256_mul_pd(xirecolor_12,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)));
                  const __m256d xirecolor_92 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_336,xirecolor_45),xirecolor_74),xirecolor_90),xirecolor_91);
                  const __m256d xi_463 = _mm256_mul_pd(xirecolor_76,xirecolor_92);
                  const __m256d xirecolor_93 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_463,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_94 = _mm256_add_pd(tmp_a12,tmp_b12);
                  const __m256d xirecolor_95 = _mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_70,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_16),xirecolor_69);
                  const __m256d xirecolor_96 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_346,xi_452),xirecolor_46),xirecolor_82),xirecolor_90),xirecolor_95);
                  const __m256d xi_464 = _mm256_mul_pd(xirecolor_76,xirecolor_96);
                  const __m256d xirecolor_97 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_464,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_98 = _mm256_add_pd(tmp_a13,tmp_b13);
                  const __m256d xirecolor_99 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(xi_354,xi_456),xirecolor_13),xirecolor_42),xirecolor_45),xirecolor_73),xirecolor_82);
                  const __m256d xi_465 = _mm256_mul_pd(xirecolor_76,xirecolor_99);
                  const __m256d xirecolor_100 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_465,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_101 = _mm256_add_pd(tmp_a14,tmp_b14);
                  const __m256d xirecolor_102 = _mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_add_pd(_mm256_mul_pd(xirecolor_49,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xi_360),xirecolor_43),xirecolor_47),xirecolor_91),xirecolor_95);
                  const __m256d xi_466 = _mm256_mul_pd(xirecolor_102,xirecolor_76);
                  const __m256d xirecolor_103 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),_mm256_mul_pd(xi_466,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),xirecolor_54));
                  const __m256d xirecolor_104 = _mm256_add_pd(tmp_a15,tmp_b15);
                  const __m256d xirecolor_105 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_464,xirecolor_54));
                  const __m256d xirecolor_106 = _mm256_add_pd(tmp_a16,tmp_b16);
                  const __m256d xirecolor_107 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_463,xirecolor_54));
                  const __m256d xirecolor_108 = _mm256_add_pd(tmp_a17,tmp_b17);
                  const __m256d xirecolor_109 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_466,xirecolor_54));
                  const __m256d xirecolor_110 = _mm256_add_pd(tmp_a18,tmp_b18);
                  const __m256d xirecolor_111 = _mm256_mul_pd(xirecolor_77,_mm256_blendv_pd(_mm256_set_pd(0.0,0.0,0.0,0.0),xi_465,xirecolor_54));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0],_mm256_mul_pd(xirecolor_0,xirecolor_3));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_4),xirecolor_57));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_58),xirecolor_59));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_60),xirecolor_61));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_62),xirecolor_63));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_64),xirecolor_65));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_66),xirecolor_67));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_68),xirecolor_78));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_79),xirecolor_84));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_85),xirecolor_86));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_87),xirecolor_88));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_89),xirecolor_93));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_94),xirecolor_97));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_3,xirecolor_98),xirecolor_100));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_101,xirecolor_3),xirecolor_103));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_104,xirecolor_3),xirecolor_105));
                  _mm256_store_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_106,xirecolor_3),xirecolor_107));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_108,xirecolor_3),xirecolor_109));
                  _mm256_storeu_pd(&_data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_110,xirecolor_3),xirecolor_111));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0],_mm256_mul_pd(xi_467,xirecolor_0));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_57,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_4)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_59,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_58)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_61,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_60)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_63,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_62)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_65,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_64)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_67,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_66)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_78,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_68)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_84,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_79)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_86,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_85)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_88,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_87)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_93,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_89)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_97,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_94)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_100,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_98)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_103,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_101)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_105,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_104)));
                  _mm256_store_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_107,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_106)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_109,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_108)));
                  _mm256_storeu_pd(&_data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0],_mm256_add_pd(_mm256_mul_pd(xirecolor_111,_mm256_set_pd(-1.0,-1.0,-1.0,-1.0)),_mm256_mul_pd(xi_467,xirecolor_110)));
               }
               for (int64_t ctr_0 = (int64_t)((_size_force_a_0 - 2) / (4)) * (4) + 1; ctr_0 < _size_force_a_0 - 1; ctr_0 += 1)
               {
                  const double xi_213 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_214 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_215 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_216 = -0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_217 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_218 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_219 = -xi_215 + xi_216 + xi_217 + xi_218;
                  const double xi_220 = -0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_221 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_222 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_223 = 0.013888888888888888*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_224 = xi_220 - xi_221 + xi_222 + xi_223;
                  const double xi_225 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_226 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_227 = -xi_225 + xi_226;
                  const double xi_228 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_229 = -xi_228;
                  const double xi_230 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_231 = xi_229 + xi_230;
                  const double xi_232 = -0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_233 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_234 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0];
                  const double xi_235 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_236 = -xi_234 + xi_235;
                  const double xi_237 = xi_232 + xi_233 + xi_236;
                  const double xi_238 = xi_231 + xi_237;
                  const double xi_239 = -xi_213 + xi_214 + xi_219 + xi_224 + xi_227 + xi_238;
                  const double xi_240 = (xi_239*xi_239);
                  const double xi_241 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_242 = -xi_241;
                  const double xi_243 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_244 = xi_229 + xi_230 + xi_242 + xi_243;
                  const double xi_245 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_246 = 0.055555555555555552*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_247 = -xi_245 + xi_246;
                  const double xi_248 = xi_213 - xi_214;
                  const double xi_249 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const double xi_250 = -xi_249;
                  const double xi_251 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const double xi_252 = xi_234 - xi_235 + xi_250 + xi_251;
                  const double xi_253 = xi_219 + xi_221 - xi_223 + xi_244 + xi_247 + xi_248 + xi_252;
                  const double xi_254 = (xi_253*xi_253);
                  const double xi_255 = xi_245 - xi_246;
                  const double xi_256 = xi_242 + xi_243;
                  const double xi_257 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1];
                  const double xi_258 = -xi_257;
                  const double xi_259 = 0.22222222222222221*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_260 = xi_225 - xi_226 + xi_258 + xi_259;
                  const double xi_261 = xi_215 + xi_216 + xi_217 - xi_218 + xi_224 + xi_248 + xi_255 + xi_256 + xi_260;
                  const double xi_262 = (xi_261*xi_261);
                  const double xi_263 = xi_240 + xi_254 + xi_262;
                  const double xi_264 = pow(xi_263, 0.5);
                  const double xi_265 = (_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]);
                  const double xi_266 = sigma*xi_264*((0.5 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_a): ((-0.5 > _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (omega_shear_b): ((0.0 < _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]) ? (xi_265*xi_8 + xi_5 + xi_9*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0]): (xi_13*xi_265 + xi_13*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0] + xi_5))));
                  const bool xi_267 = xi_264 > 0.0;
                  const double xi_268 = ((xi_267) ? (xi_266*0.25): (0.0));
                  const double xi_280 = ((1.0) / (xi_263));
                  const double xi_281 = xi_266*1.125;
                  const double xi_282 = ((xi_267) ? (xi_281*(xi_254*xi_280*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_286 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552*(xi_253*xi_253) - 0.018518518518518517)): (0.0));
                  const double xi_294 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552*(xi_261*xi_261) - 0.018518518518518517)): (0.0));
                  const double xi_297 = ((xi_267) ? (xi_281*(xi_262*xi_280*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_303 = ((xi_267) ? (xi_281*(xi_280*0.055555555555555552*(xi_239*xi_239) - 0.018518518518518517)): (0.0));
                  const double xi_305 = ((xi_267) ? (xi_281*(xi_240*xi_280*0.055555555555555552 - 0.018518518518518517)): (0.0));
                  const double xi_313 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_314 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_315 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_316 = xi_220 + xi_222;
                  const double xi_317 = xi_260 + xi_316;
                  const double xi_318 = -0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_319 = xi_249 - xi_251 + xi_318;
                  const double xi_320 = xi_228 - xi_230 + xi_236 - xi_313 + xi_314 + xi_315 + xi_317 + xi_319;
                  const double xi_321 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_320*xi_320) - 0.037037037037037035)): (0.0));
                  const double xi_325 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_326 = 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_327 = -0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.027777777777777776*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_328 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + ctr_0 + 1];
                  const double xi_329 = xi_231 + xi_252 + xi_317 - xi_325 + xi_326 + xi_327 + xi_328;
                  const double xi_330 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_329*xi_329) - 0.037037037037037035)): (0.0));
                  const double xi_331 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_329*xi_329) - 0.037037037037037035)): (0.0));
                  const double xi_333 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_320*xi_320) - 0.037037037037037035)): (0.0));
                  const double xi_336 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_337 = xi_227 + xi_232 + xi_233 + xi_316;
                  const double xi_338 = xi_241 - xi_243 + xi_325 - xi_326;
                  const double xi_339 = xi_255 + xi_319 + xi_336 + xi_337 + xi_338;
                  const double xi_340 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_339*xi_339) - 0.037037037037037035)): (0.0));
                  const double xi_346 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 - _stride_phasefield_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0];
                  const double xi_347 = xi_247 + xi_313 - xi_314;
                  const double xi_348 = xi_250 + xi_251 + xi_256 + xi_327 + xi_337 + xi_346 + xi_347;
                  const double xi_349 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_348*xi_348) - 0.037037037037037035)): (0.0));
                  const double xi_354 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 - 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 + 1];
                  const double xi_355 = xi_237 + xi_244 + xi_255 + xi_258 + xi_259 + xi_318 + xi_327 + xi_354;
                  const double xi_356 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_355*xi_355) - 0.037037037037037035)): (0.0));
                  const double xi_360 = -0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + _stride_phasefield_2 + ctr_0 + 1] + 0.1111111111111111*_data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 - _stride_phasefield_2 + ctr_0 - 1];
                  const double xi_361 = xi_238 + xi_257 - xi_259 + xi_338 + xi_347 + xi_360;
                  const double xi_362 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_361*xi_361) - 0.037037037037037035)): (0.0));
                  const double xi_363 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_348*xi_348) - 0.037037037037037035)): (0.0));
                  const double xi_365 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_339*xi_339) - 0.037037037037037035)): (0.0));
                  const double xi_366 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_361*xi_361) - 0.037037037037037035)): (0.0));
                  const double xi_368 = ((xi_267) ? (xi_281*(xi_280*0.027777777777777776*(xi_355*xi_355) - 0.037037037037037035)): (0.0));
                  const double xia_3_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0];
                  const double xia_4_collide = _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const bool xi_0 = xia_4_collide > 0.5;
                  const bool xi_1 = xia_4_collide < -0.5;
                  const double xi_10 = (xia_4_collide*xia_4_collide);
                  const bool xi_11 = xia_4_collide > 0.0;
                  const double xia_5_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0];
                  const double xia_6_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0];
                  const double xia_7_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0];
                  const double xia_8_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0];
                  const double xia_9_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + ctr_0];
                  const double xi_42 = xia_9_collide*0.16666666666666666;
                  const double xi_67 = xia_9_collide*0.083333333333333329;
                  const double xia_10_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0];
                  const double xi_271 = -xia_10_collide;
                  const double xia_11_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xia_12_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0];
                  const double xi_290 = -xia_12_collide;
                  const double xia_13_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0];
                  const double xi_198 = xia_13_collide + xia_7_collide;
                  const double xia_14_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0];
                  const double xi_199 = xi_198 + xia_10_collide + xia_14_collide;
                  const double xi_272 = xi_271 + xia_14_collide;
                  const double xia_15_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0];
                  const double xia_16_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0];
                  const double xi_291 = xi_290 + xia_16_collide;
                  const double xia_17_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + _stride_force_a_3 + ctr_0];
                  const double xi_26 = xia_17_collide*0.16666666666666666;
                  const double xi_54 = xia_17_collide*0.083333333333333329;
                  const double xi_58 = xia_17_collide*0.25;
                  const double xia_18_collide = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xia_19_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xia_20_collide = _data_force_a[_stride_force_a_1*ctr_1 + _stride_force_a_2*ctr_2 + 2*_stride_force_a_3 + ctr_0];
                  const double xi_51 = xia_20_collide*0.16666666666666666;
                  const double xi_77 = xia_20_collide*0.083333333333333329;
                  const double xia_21_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xia_22_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0];
                  const double xia_23_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0];
                  const double xia_24_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0];
                  const double xia_25_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0];
                  const double xia_26_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0];
                  const double xi_202 = xia_26_collide + xia_5_collide;
                  const double xi_203 = xi_202 + xia_12_collide + xia_16_collide;
                  const double xia_27_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0];
                  const double xia_28_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0];
                  const double xi_205 = xia_28_collide + xia_6_collide;
                  const double xia_29_collide = _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0];
                  const double xi_206 = xi_205 + xia_15_collide + xia_29_collide;
                  const double xi_276 = -xia_29_collide;
                  const double xi_277 = xi_276 + xia_15_collide;
                  const double xia_0_collide = ((1.0) / (-0.25*((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_8 + xi_5 + xi_9*xia_4_collide): (xi_10*xi_13 + xi_13*xia_4_collide + xi_5)))) + 2.0));
                  const double xi_19 = ((xi_0) ? (omega_shear_a): ((xi_1) ? (omega_shear_b): ((xi_11) ? (xi_10*xi_17 + xi_15 - xi_17*xia_4_collide): (xi_10*xi_18 + xi_15 + xi_18*xia_4_collide))));
                  const double rr_1_a_collide = xi_19;
                  const double xi_23 = rr_1_a_collide*0.5;
                  const double xi_62 = rr_1_a_collide*0.125;
                  const double xi_63 = xi_62*xia_17_collide;
                  const double xi_71 = rr_1_a_collide*0.041666666666666664;
                  const double xi_307 = rr_1_a_collide*0.25;
                  const double rr_2_a_collide = xia_0_collide*(xi_19*-2.0 + 4.0);
                  const double xi_28 = rr_2_a_collide*xi_26 + xi_26 - xi_27*xia_17_collide;
                  const double xi_43 = rr_2_a_collide*xi_42 - xi_27*xia_9_collide + xi_42;
                  const double xi_52 = rr_2_a_collide*xi_51 - xi_27*xia_20_collide + xi_51;
                  const double xi_55 = rr_2_a_collide*0.041666666666666664;
                  const double xi_56 = xi_55*xia_17_collide;
                  const double xi_57 = xi_54 - xi_56;
                  const double xi_68 = xi_55*xia_9_collide;
                  const double xi_69 = -xi_67 + xi_68;
                  const double xi_73 = -xi_54 + xi_56;
                  const double xi_75 = xi_67 - xi_68;
                  const double xi_78 = xi_55*xia_20_collide;
                  const double xi_79 = xi_77 - xi_78;
                  const double xi_98 = -xi_77 + xi_78;
                  const double xi_274 = rr_2_a_collide*0.5;
                  const double xi_310 = rr_2_a_collide*0.25;
                  const double rho_a_collide = xia_18_collide;
                  const double xi_195 = rho_a_collide*-0.1111111111111111;
                  const double xi_208 = rho_a_collide*-0.33333333333333331;
                  const double xi_210 = xi_199 + xi_208;
                  const double xi_269 = rho_a_collide*0.33333333333333331;
                  const double u_0_a_collide = xia_11_collide;
                  const double xi_20 = u_0_a_collide*xia_9_collide;
                  const double xi_33 = xi_20*0.16666666666666666;
                  const double xi_34 = xi_20*0.083333333333333329;
                  const double xi_35 = rr_1_a_collide*xi_34;
                  const double xi_36 = -xi_33 + xi_35;
                  const double xi_39 = xi_33 - xi_35;
                  const double xi_41 = rr_1_a_collide*xi_33;
                  const double xi_59 = -u_0_a_collide*xi_58;
                  const double xi_64 = u_0_a_collide*xi_63;
                  const double xi_80 = -xi_34;
                  const double xi_81 = xi_20*xi_71;
                  const double xi_90 = u_0_a_collide*xia_20_collide;
                  const double xi_91 = xi_90*-0.25;
                  const double xi_94 = xi_62*xi_90;
                  const double xi_194 = (u_0_a_collide*u_0_a_collide);
                  const double xi_201 = rho_a_collide*xi_194*-0.33333333333333331;
                  const double xi_209 = rr_1_a_collide*(rho_a_collide*xi_194 - xi_203 - xi_206 - xi_208 - xia_22_collide - xia_25_collide);
                  const double xi_287 = u_0_a_collide*xi_269;
                  const double xi_288 = xi_277 + xi_287 - xia_28_collide + xia_6_collide;
                  const double xi_289 = xi_274*xi_288;
                  const double xi_292 = xi_287 + xi_291 + xia_26_collide - xia_5_collide;
                  const double xi_293 = xi_274*xi_292;
                  const double xi_311 = xi_288*xi_310;
                  const double xi_312 = -xi_311;
                  const double xi_352 = xi_292*xi_310;
                  const double xi_353 = -xi_352;
                  const double u_1_a_collide = xia_21_collide;
                  const double xi_21 = u_1_a_collide*xia_17_collide;
                  const double xi_24 = xi_21*0.16666666666666666;
                  const double xi_25 = rr_1_a_collide*xi_24;
                  const double xi_44 = xi_21*0.083333333333333329;
                  const double xi_45 = rr_1_a_collide*xi_44;
                  const double xi_46 = xi_24 - xi_45;
                  const double xi_48 = -xi_24 + xi_45;
                  const double xi_50 = xi_36 + xi_48;
                  const double xi_53 = xi_39 + xi_46;
                  const double xi_60 = u_1_a_collide*0.25;
                  const double xi_61 = -xi_60*xia_9_collide;
                  const double xi_65 = u_1_a_collide*xi_62;
                  const double xi_66 = xi_65*xia_9_collide;
                  const double xi_70 = xi_59 + xi_61 + xi_64 + xi_66 + xi_69;
                  const double xi_76 = xi_59 + xi_61 + xi_64 + xi_66 + xi_75;
                  const double xi_83 = xi_60*xia_20_collide;
                  const double xi_85 = xi_65*xia_20_collide;
                  const double xi_192 = (u_1_a_collide*u_1_a_collide);
                  const double xi_193 = rho_a_collide*xi_192*-0.33333333333333331;
                  const double xi_211 = rr_1_a_collide*(rho_a_collide*xi_192 - xi_206 - xi_210 - xia_24_collide - xia_8_collide);
                  const double xi_270 = u_1_a_collide*xi_269;
                  const double xi_273 = xi_270 + xi_272 + xia_13_collide - xia_7_collide;
                  const double xi_275 = xi_273*xi_274;
                  const double xi_278 = xi_270 + xi_277 + xia_28_collide - xia_6_collide;
                  const double xi_279 = xi_274*xi_278;
                  const double xi_306 = rho_a_collide*u_1_a_collide;
                  const double xi_308 = xi_307*(u_0_a_collide*xi_306 + xi_205 + xi_276 - xia_15_collide);
                  const double xi_309 = -xi_308;
                  const double xi_323 = xi_278*xi_310;
                  const double xi_335 = xi_273*xi_310;
                  const double xi_345 = -xi_335;
                  const double u_2_a_collide = xia_19_collide;
                  const double xi_22 = u_2_a_collide*xia_20_collide;
                  const double xi_29 = xi_22*0.16666666666666666;
                  const double xi_30 = xi_22*0.083333333333333329;
                  const double xi_31 = rr_1_a_collide*xi_30;
                  const double xi_32 = -xi_29 + xi_31;
                  const double xi_37 = xi_32 + xi_36;
                  const double xi_38 = xi_29 - xi_31;
                  const double xi_40 = xi_38 + xi_39;
                  const double xi_47 = xi_38 + xi_46;
                  const double xi_49 = rr_1_a_collide*xi_29;
                  const double xi_72 = xi_22*xi_71 - xi_30 + xi_53;
                  const double xi_74 = rr_1_a_collide*u_2_a_collide*xia_20_collide*-0.041666666666666664 + xi_30 + xi_50;
                  const double xi_82 = xi_47 + xi_57 + xi_80 + xi_81;
                  const double xi_84 = u_2_a_collide*xi_58;
                  const double xi_86 = u_2_a_collide*xi_63;
                  const double xi_87 = xi_83 + xi_84 - xi_85 - xi_86;
                  const double xi_88 = -xi_83 - xi_84 + xi_85 + xi_86;
                  const double xi_89 = xi_47 + xi_73 + xi_80 + xi_81;
                  const double xi_92 = u_2_a_collide*xia_9_collide;
                  const double xi_93 = xi_92*-0.25;
                  const double xi_95 = xi_62*xi_92;
                  const double xi_96 = xi_69 + xi_91 + xi_93 + xi_94 + xi_95;
                  const double xi_97 = xi_21*xi_71 + xi_40 - xi_44;
                  const double xi_99 = rr_1_a_collide*u_1_a_collide*xia_17_collide*-0.041666666666666664 + xi_37 + xi_44;
                  const double xi_100 = xi_75 + xi_91 + xi_93 + xi_94 + xi_95;
                  const double xi_196 = (u_2_a_collide*u_2_a_collide);
                  const double xi_197 = rho_a_collide*xi_196*-0.33333333333333331 + xi_195;
                  const double xi_200 = rr_1_a_collide*(rho_a_collide*xi_194*-0.16666666666666666 - xi_193 - xi_197 - xi_199);
                  const double xi_204 = rr_1_a_collide*(rho_a_collide*xi_192*-0.16666666666666666 - xi_197 - xi_201 - xi_203);
                  const double xi_207 = rr_1_a_collide*(rho_a_collide*xi_196*-0.16666666666666666 - xi_193 - xi_195 - xi_201 - xi_206);
                  const double xi_212 = rr_1_a_collide*(rho_a_collide*xi_196 - xi_203 - xi_210 - xia_23_collide - xia_3_collide);
                  const double xi_283 = xi_200*-0.5;
                  const double xi_284 = xi_207*-0.5;
                  const double xi_285 = xi_211*0.5 + xi_283 + xi_284;
                  const double xi_295 = xi_204*-0.5;
                  const double xi_296 = xi_209*0.5 + xi_284 + xi_295;
                  const double xi_298 = u_2_a_collide*xi_269;
                  const double xi_299 = xi_272 + xi_298 - xia_13_collide + xia_7_collide;
                  const double xi_300 = xi_274*xi_299;
                  const double xi_301 = xi_291 + xi_298 - xia_26_collide + xia_5_collide;
                  const double xi_302 = xi_274*xi_301;
                  const double xi_304 = xi_212*0.5 + xi_283 + xi_295;
                  const double xi_322 = xi_207*0.25;
                  const double xi_324 = xi_322 + xi_323;
                  const double xi_332 = xi_322 - xi_323;
                  const double xi_334 = xi_307*(u_2_a_collide*xi_306 + xi_198 + xi_271 - xia_14_collide);
                  const double xi_341 = xi_200*0.25;
                  const double xi_342 = xi_299*xi_310;
                  const double xi_343 = xi_341 + xi_342;
                  const double xi_344 = -xi_334;
                  const double xi_350 = xi_307*(rho_a_collide*u_0_a_collide*u_2_a_collide + xi_202 + xi_290 - xia_16_collide);
                  const double xi_351 = -xi_350;
                  const double xi_357 = xi_204*0.25;
                  const double xi_358 = xi_301*xi_310;
                  const double xi_359 = xi_357 + xi_358;
                  const double xi_364 = xi_341 - xi_342;
                  const double xi_367 = xi_357 - xi_358;
                  const double forceTerm_0_a_collide = xi_20*xi_23 - xi_20 + xi_21*xi_23 - xi_21 + xi_22*xi_23 - xi_22;
                  const double forceTerm_1_a_collide = xi_21*0.33333333333333331 - xi_25 + xi_28 + xi_37;
                  const double forceTerm_2_a_collide = u_1_a_collide*xia_17_collide*0.33333333333333331 - xi_25 - xi_28 - xi_40;
                  const double forceTerm_3_a_collide = u_0_a_collide*xia_9_collide*0.33333333333333331 - xi_41 - xi_43 - xi_47;
                  const double forceTerm_4_a_collide = xi_20*0.33333333333333331 + xi_32 - xi_41 + xi_43 + xi_48;
                  const double forceTerm_5_a_collide = xi_22*0.33333333333333331 - xi_49 + xi_50 + xi_52;
                  const double forceTerm_6_a_collide = u_2_a_collide*xia_20_collide*0.33333333333333331 - xi_49 - xi_52 - xi_53;
                  const double forceTerm_7_a_collide = xi_57 + xi_70 + xi_72;
                  const double forceTerm_8_a_collide = -xi_70 - xi_73 - xi_74;
                  const double forceTerm_9_a_collide = -xi_57 - xi_74 - xi_76;
                  const double forceTerm_10_a_collide = xi_72 + xi_73 + xi_76;
                  const double forceTerm_11_a_collide = xi_79 + xi_82 + xi_87;
                  const double forceTerm_12_a_collide = xi_79 + xi_88 + xi_89;
                  const double forceTerm_13_a_collide = xi_79 + xi_96 + xi_97;
                  const double forceTerm_14_a_collide = -xi_96 - xi_98 - xi_99;
                  const double forceTerm_15_a_collide = xi_82 + xi_88 + xi_98;
                  const double forceTerm_16_a_collide = xi_87 + xi_89 + xi_98;
                  const double forceTerm_17_a_collide = -xi_100 - xi_79 - xi_99;
                  const double forceTerm_18_a_collide = xi_100 + xi_97 + xi_98;
                  const double xib_3_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0];
                  const double xib_4_collide = _data_phasefield[_stride_phasefield_1*ctr_1 + _stride_phasefield_2*ctr_2 + ctr_0];
                  const bool xi_101 = xib_4_collide > 0.5;
                  const bool xi_102 = xib_4_collide < -0.5;
                  const double xi_103 = (xib_4_collide*xib_4_collide);
                  const bool xi_104 = xib_4_collide > 0.0;
                  const double xib_5_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + ctr_0];
                  const double xi_133 = xib_5_collide*0.16666666666666666;
                  const double xi_158 = xib_5_collide*0.083333333333333329;
                  const double xib_6_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0];
                  const double xib_7_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0];
                  const double xib_8_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0];
                  const double xib_9_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0];
                  const double xib_10_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0];
                  const double xi_405 = -xib_10_collide;
                  const double xib_11_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + ctr_0];
                  const double xib_12_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0];
                  const double xib_13_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0];
                  const double xib_14_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0];
                  const double xib_15_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0];
                  const double xi_392 = -xib_15_collide;
                  const double xib_16_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0];
                  const double xib_17_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + 2*_stride_velocity_3 + ctr_0];
                  const double xib_18_collide = _data_velocity[_stride_velocity_1*ctr_1 + _stride_velocity_2*ctr_2 + _stride_velocity_3 + ctr_0];
                  const double xib_19_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0];
                  const double xib_20_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0];
                  const double xi_382 = xib_20_collide + xib_3_collide;
                  const double xib_21_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + 2*_stride_force_b_3 + ctr_0];
                  const double xi_142 = xib_21_collide*0.16666666666666666;
                  const double xi_168 = xib_21_collide*0.083333333333333329;
                  const double xib_22_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0];
                  const double xi_406 = xi_405 + xib_22_collide;
                  const double xib_23_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0];
                  const double xi_383 = xi_382 + xib_23_collide + xib_8_collide;
                  const double xi_397 = -xib_23_collide;
                  const double xi_398 = xi_397 + xib_8_collide;
                  const double xib_24_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0];
                  const double xi_375 = xib_16_collide + xib_24_collide;
                  const double xib_25_collide = _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xib_26_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0];
                  const double xi_376 = xi_375 + xib_15_collide + xib_26_collide;
                  const double xi_393 = xi_392 + xib_26_collide;
                  const double xib_27_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0];
                  const double xib_28_collide = _data_force_b[_stride_force_b_1*ctr_1 + _stride_force_b_2*ctr_2 + _stride_force_b_3 + ctr_0];
                  const double xi_117 = xib_28_collide*0.16666666666666666;
                  const double xi_145 = xib_28_collide*0.083333333333333329;
                  const double xi_149 = xib_28_collide*0.25;
                  const double xib_29_collide = _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0];
                  const double xi_379 = xib_19_collide + xib_29_collide;
                  const double xi_380 = xi_379 + xib_10_collide + xib_22_collide;
                  const double xib_0_collide = ((1.0) / (-0.25*((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_8 + xi_5 + xi_9*xib_4_collide): (xi_103*xi_13 + xi_13*xib_4_collide + xi_5)))) + 2.0));
                  const double xi_110 = ((xi_101) ? (omega_shear_a): ((xi_102) ? (omega_shear_b): ((xi_104) ? (xi_103*xi_108 + xi_106 - xi_108*xib_4_collide): (xi_103*xi_109 + xi_106 + xi_109*xib_4_collide))));
                  const double rr_1_b_collide = xi_110;
                  const double xi_114 = rr_1_b_collide*0.5;
                  const double xi_153 = rr_1_b_collide*0.125;
                  const double xi_154 = xi_153*xib_28_collide;
                  const double xi_162 = rr_1_b_collide*0.041666666666666664;
                  const double xi_420 = rr_1_b_collide*0.25;
                  const double rr_2_b_collide = xib_0_collide*(xi_110*-2.0 + 4.0);
                  const double xi_119 = rr_2_b_collide*xi_117 + xi_117 - xi_118*xib_28_collide;
                  const double xi_134 = rr_2_b_collide*xi_133 - xi_118*xib_5_collide + xi_133;
                  const double xi_143 = rr_2_b_collide*xi_142 - xi_118*xib_21_collide + xi_142;
                  const double xi_146 = rr_2_b_collide*0.041666666666666664;
                  const double xi_147 = xi_146*xib_28_collide;
                  const double xi_148 = xi_145 - xi_147;
                  const double xi_159 = xi_146*xib_5_collide;
                  const double xi_160 = -xi_158 + xi_159;
                  const double xi_164 = -xi_145 + xi_147;
                  const double xi_166 = xi_158 - xi_159;
                  const double xi_169 = xi_146*xib_21_collide;
                  const double xi_170 = xi_168 - xi_169;
                  const double xi_189 = -xi_168 + xi_169;
                  const double xi_395 = rr_2_b_collide*0.5;
                  const double xi_423 = rr_2_b_collide*0.25;
                  const double rho_b_collide = xib_25_collide;
                  const double xi_372 = rho_b_collide*-0.1111111111111111;
                  const double xi_385 = rho_b_collide*-0.33333333333333331;
                  const double xi_386 = xi_380 + xi_385;
                  const double xi_390 = rho_b_collide*0.33333333333333331;
                  const double u_0_b_collide = xib_11_collide;
                  const double xi_111 = u_0_b_collide*xib_5_collide;
                  const double xi_124 = xi_111*0.16666666666666666;
                  const double xi_125 = xi_111*0.083333333333333329;
                  const double xi_126 = rr_1_b_collide*xi_125;
                  const double xi_127 = -xi_124 + xi_126;
                  const double xi_130 = xi_124 - xi_126;
                  const double xi_132 = rr_1_b_collide*xi_124;
                  const double xi_150 = -u_0_b_collide*xi_149;
                  const double xi_155 = u_0_b_collide*xi_154;
                  const double xi_171 = -xi_125;
                  const double xi_172 = xi_111*xi_162;
                  const double xi_181 = u_0_b_collide*xib_21_collide;
                  const double xi_182 = xi_181*-0.25;
                  const double xi_185 = xi_153*xi_181;
                  const double xi_371 = (u_0_b_collide*u_0_b_collide);
                  const double xi_378 = rho_b_collide*xi_371*-0.33333333333333331;
                  const double xi_387 = rr_1_b_collide*(rho_b_collide*xi_371 - xi_383 - xi_386 - xib_12_collide - xib_27_collide);
                  const double xi_404 = u_0_b_collide*xi_390;
                  const double xi_407 = xi_404 + xi_406 + xib_19_collide - xib_29_collide;
                  const double xi_408 = xi_395*xi_407;
                  const double xi_409 = xi_398 + xi_404 - xib_20_collide + xib_3_collide;
                  const double xi_410 = xi_395*xi_409;
                  const double xi_424 = xi_409*xi_423;
                  const double xi_425 = -xi_424;
                  const double xi_439 = xi_407*xi_423;
                  const double xi_440 = -xi_439;
                  const double u_1_b_collide = xib_18_collide;
                  const double xi_112 = u_1_b_collide*xib_28_collide;
                  const double xi_115 = xi_112*0.16666666666666666;
                  const double xi_116 = rr_1_b_collide*xi_115;
                  const double xi_135 = xi_112*0.083333333333333329;
                  const double xi_136 = rr_1_b_collide*xi_135;
                  const double xi_137 = xi_115 - xi_136;
                  const double xi_139 = -xi_115 + xi_136;
                  const double xi_141 = xi_127 + xi_139;
                  const double xi_144 = xi_130 + xi_137;
                  const double xi_151 = u_1_b_collide*0.25;
                  const double xi_152 = -xi_151*xib_5_collide;
                  const double xi_156 = u_1_b_collide*xi_153;
                  const double xi_157 = xi_156*xib_5_collide;
                  const double xi_161 = xi_150 + xi_152 + xi_155 + xi_157 + xi_160;
                  const double xi_167 = xi_150 + xi_152 + xi_155 + xi_157 + xi_166;
                  const double xi_174 = xi_151*xib_21_collide;
                  const double xi_176 = xi_156*xib_21_collide;
                  const double xi_369 = (u_1_b_collide*u_1_b_collide);
                  const double xi_370 = rho_b_collide*xi_369*-0.33333333333333331;
                  const double xi_388 = rr_1_b_collide*(rho_b_collide*xi_369 - xi_376 - xi_383 - xi_385 - xib_14_collide - xib_6_collide);
                  const double xi_391 = u_1_b_collide*xi_390;
                  const double xi_394 = xi_391 + xi_393 + xib_16_collide - xib_24_collide;
                  const double xi_396 = xi_394*xi_395;
                  const double xi_399 = xi_391 + xi_398 + xib_20_collide - xib_3_collide;
                  const double xi_400 = xi_395*xi_399;
                  const double xi_419 = rho_b_collide*u_1_b_collide;
                  const double xi_421 = xi_420*(u_0_b_collide*xi_419 + xi_382 + xi_397 - xib_8_collide);
                  const double xi_422 = -xi_421;
                  const double xi_427 = xi_399*xi_423;
                  const double xi_431 = xi_394*xi_423;
                  const double xi_436 = -xi_431;
                  const double u_2_b_collide = xib_17_collide;
                  const double xi_113 = u_2_b_collide*xib_21_collide;
                  const double xi_120 = xi_113*0.16666666666666666;
                  const double xi_121 = xi_113*0.083333333333333329;
                  const double xi_122 = rr_1_b_collide*xi_121;
                  const double xi_123 = -xi_120 + xi_122;
                  const double xi_128 = xi_123 + xi_127;
                  const double xi_129 = xi_120 - xi_122;
                  const double xi_131 = xi_129 + xi_130;
                  const double xi_138 = xi_129 + xi_137;
                  const double xi_140 = rr_1_b_collide*xi_120;
                  const double xi_163 = xi_113*xi_162 - xi_121 + xi_144;
                  const double xi_165 = rr_1_b_collide*u_2_b_collide*xib_21_collide*-0.041666666666666664 + xi_121 + xi_141;
                  const double xi_173 = xi_138 + xi_148 + xi_171 + xi_172;
                  const double xi_175 = u_2_b_collide*xi_149;
                  const double xi_177 = u_2_b_collide*xi_154;
                  const double xi_178 = xi_174 + xi_175 - xi_176 - xi_177;
                  const double xi_179 = -xi_174 - xi_175 + xi_176 + xi_177;
                  const double xi_180 = xi_138 + xi_164 + xi_171 + xi_172;
                  const double xi_183 = u_2_b_collide*xib_5_collide;
                  const double xi_184 = xi_183*-0.25;
                  const double xi_186 = xi_153*xi_183;
                  const double xi_187 = xi_160 + xi_182 + xi_184 + xi_185 + xi_186;
                  const double xi_188 = xi_112*xi_162 + xi_131 - xi_135;
                  const double xi_190 = rr_1_b_collide*u_1_b_collide*xib_28_collide*-0.041666666666666664 + xi_128 + xi_135;
                  const double xi_191 = xi_166 + xi_182 + xi_184 + xi_185 + xi_186;
                  const double xi_373 = (u_2_b_collide*u_2_b_collide);
                  const double xi_374 = rho_b_collide*xi_373*-0.33333333333333331 + xi_372;
                  const double xi_377 = rr_1_b_collide*(rho_b_collide*xi_371*-0.16666666666666666 - xi_370 - xi_374 - xi_376);
                  const double xi_381 = rr_1_b_collide*(rho_b_collide*xi_369*-0.16666666666666666 - xi_374 - xi_378 - xi_380);
                  const double xi_384 = rr_1_b_collide*(rho_b_collide*xi_373*-0.16666666666666666 - xi_370 - xi_372 - xi_378 - xi_383);
                  const double xi_389 = rr_1_b_collide*(rho_b_collide*xi_373 - xi_376 - xi_386 - xib_7_collide - xib_9_collide);
                  const double xi_401 = xi_377*-0.5;
                  const double xi_402 = xi_384*-0.5;
                  const double xi_403 = xi_388*0.5 + xi_401 + xi_402;
                  const double xi_411 = xi_381*-0.5;
                  const double xi_412 = xi_387*0.5 + xi_402 + xi_411;
                  const double xi_413 = u_2_b_collide*xi_390;
                  const double xi_414 = xi_406 + xi_413 - xib_19_collide + xib_29_collide;
                  const double xi_415 = xi_395*xi_414;
                  const double xi_416 = xi_393 + xi_413 - xib_16_collide + xib_24_collide;
                  const double xi_417 = xi_395*xi_416;
                  const double xi_418 = xi_389*0.5 + xi_401 + xi_411;
                  const double xi_426 = xi_384*0.25;
                  const double xi_428 = xi_426 + xi_427;
                  const double xi_429 = xi_426 - xi_427;
                  const double xi_430 = xi_420*(u_2_b_collide*xi_419 + xi_375 + xi_392 - xib_26_collide);
                  const double xi_432 = xi_377*0.25;
                  const double xi_433 = xi_416*xi_423;
                  const double xi_434 = xi_432 + xi_433;
                  const double xi_435 = -xi_430;
                  const double xi_437 = xi_420*(rho_b_collide*u_0_b_collide*u_2_b_collide + xi_379 + xi_405 - xib_22_collide);
                  const double xi_438 = -xi_437;
                  const double xi_441 = xi_381*0.25;
                  const double xi_442 = xi_414*xi_423;
                  const double xi_443 = xi_441 + xi_442;
                  const double xi_444 = xi_432 - xi_433;
                  const double xi_445 = xi_441 - xi_442;
                  const double forceTerm_0_b_collide = xi_111*xi_114 - xi_111 + xi_112*xi_114 - xi_112 + xi_113*xi_114 - xi_113;
                  const double forceTerm_1_b_collide = xi_112*0.33333333333333331 - xi_116 + xi_119 + xi_128;
                  const double forceTerm_2_b_collide = u_1_b_collide*xib_28_collide*0.33333333333333331 - xi_116 - xi_119 - xi_131;
                  const double forceTerm_3_b_collide = u_0_b_collide*xib_5_collide*0.33333333333333331 - xi_132 - xi_134 - xi_138;
                  const double forceTerm_4_b_collide = xi_111*0.33333333333333331 + xi_123 - xi_132 + xi_134 + xi_139;
                  const double forceTerm_5_b_collide = xi_113*0.33333333333333331 - xi_140 + xi_141 + xi_143;
                  const double forceTerm_6_b_collide = u_2_b_collide*xib_21_collide*0.33333333333333331 - xi_140 - xi_143 - xi_144;
                  const double forceTerm_7_b_collide = xi_148 + xi_161 + xi_163;
                  const double forceTerm_8_b_collide = -xi_161 - xi_164 - xi_165;
                  const double forceTerm_9_b_collide = -xi_148 - xi_165 - xi_167;
                  const double forceTerm_10_b_collide = xi_163 + xi_164 + xi_167;
                  const double forceTerm_11_b_collide = xi_170 + xi_173 + xi_178;
                  const double forceTerm_12_b_collide = xi_170 + xi_179 + xi_180;
                  const double forceTerm_13_b_collide = xi_170 + xi_187 + xi_188;
                  const double forceTerm_14_b_collide = -xi_187 - xi_189 - xi_190;
                  const double forceTerm_15_b_collide = xi_173 + xi_179 + xi_189;
                  const double forceTerm_16_b_collide = xi_178 + xi_180 + xi_189;
                  const double forceTerm_17_b_collide = -xi_170 - xi_190 - xi_191;
                  const double forceTerm_18_b_collide = xi_188 + xi_189 + xi_191;
                  const double tmp_a0 = forceTerm_0_a_collide + xi_200 + xi_204 + xi_207 - xi_209 - xi_211 - xi_212 + xi_268 + xia_27_collide;
                  const double tmp_a1 = forceTerm_1_a_collide - xi_275 - xi_279 + xi_282 + xi_285 + xia_8_collide;
                  const double tmp_a2 = forceTerm_2_a_collide + xi_275 + xi_279 + xi_285 + xi_286 + xia_24_collide;
                  const double tmp_a3 = forceTerm_3_a_collide + xi_289 + xi_293 + xi_294 + xi_296 + xia_25_collide;
                  const double tmp_a4 = forceTerm_4_a_collide - xi_289 - xi_293 + xi_296 + xi_297 + xia_22_collide;
                  const double tmp_a5 = forceTerm_5_a_collide - xi_300 - xi_302 + xi_303 + xi_304 + xia_23_collide;
                  const double tmp_a6 = forceTerm_6_a_collide + xi_300 + xi_302 + xi_304 + xi_305 + xia_3_collide;
                  const double tmp_a7 = forceTerm_7_a_collide + xi_309 + xi_312 + xi_321 + xi_324 + xia_6_collide;
                  const double tmp_a8 = forceTerm_8_a_collide + xi_308 + xi_311 + xi_324 + xi_330 + xia_29_collide;
                  const double tmp_a9 = forceTerm_9_a_collide + xi_308 + xi_312 + xi_331 + xi_332 + xia_15_collide;
                  const double tmp_a10 = forceTerm_10_a_collide + xi_309 + xi_311 + xi_332 + xi_333 + xia_28_collide;
                  const double tmp_a11 = forceTerm_11_a_collide + xi_334 + xi_335 + xi_340 + xi_343 + xia_10_collide;
                  const double tmp_a12 = forceTerm_12_a_collide + xi_343 + xi_344 + xi_345 + xi_349 + xia_13_collide;
                  const double tmp_a13 = forceTerm_13_a_collide + xi_351 + xi_353 + xi_356 + xi_359 + xia_26_collide;
                  const double tmp_a14 = forceTerm_14_a_collide + xi_350 + xi_352 + xi_359 + xi_362 + xia_12_collide;
                  const double tmp_a15 = forceTerm_15_a_collide + xi_335 + xi_344 + xi_363 + xi_364 + xia_7_collide;
                  const double tmp_a16 = forceTerm_16_a_collide + xi_334 + xi_345 + xi_364 + xi_365 + xia_14_collide;
                  const double tmp_a17 = forceTerm_17_a_collide + xi_350 + xi_353 + xi_366 + xi_367 + xia_16_collide;
                  const double tmp_a18 = forceTerm_18_a_collide + xi_351 + xi_352 + xi_367 + xi_368 + xia_5_collide;
                  const double tmp_b0 = forceTerm_0_b_collide + xi_268 + xi_377 + xi_381 + xi_384 - xi_387 - xi_388 - xi_389 + xib_13_collide;
                  const double tmp_b1 = forceTerm_1_b_collide + xi_282 - xi_396 - xi_400 + xi_403 + xib_6_collide;
                  const double tmp_b2 = forceTerm_2_b_collide + xi_286 + xi_396 + xi_400 + xi_403 + xib_14_collide;
                  const double tmp_b3 = forceTerm_3_b_collide + xi_294 + xi_408 + xi_410 + xi_412 + xib_27_collide;
                  const double tmp_b4 = forceTerm_4_b_collide + xi_297 - xi_408 - xi_410 + xi_412 + xib_12_collide;
                  const double tmp_b5 = forceTerm_5_b_collide + xi_303 - xi_415 - xi_417 + xi_418 + xib_7_collide;
                  const double tmp_b6 = forceTerm_6_b_collide + xi_305 + xi_415 + xi_417 + xi_418 + xib_9_collide;
                  const double tmp_b7 = forceTerm_7_b_collide + xi_321 + xi_422 + xi_425 + xi_428 + xib_3_collide;
                  const double tmp_b8 = forceTerm_8_b_collide + xi_330 + xi_421 + xi_424 + xi_428 + xib_23_collide;
                  const double tmp_b9 = forceTerm_9_b_collide + xi_331 + xi_421 + xi_425 + xi_429 + xib_8_collide;
                  const double tmp_b10 = forceTerm_10_b_collide + xi_333 + xi_422 + xi_424 + xi_429 + xib_20_collide;
                  const double tmp_b11 = forceTerm_11_b_collide + xi_340 + xi_430 + xi_431 + xi_434 + xib_15_collide;
                  const double tmp_b12 = forceTerm_12_b_collide + xi_349 + xi_434 + xi_435 + xi_436 + xib_16_collide;
                  const double tmp_b13 = forceTerm_13_b_collide + xi_356 + xi_438 + xi_440 + xi_443 + xib_19_collide;
                  const double tmp_b14 = forceTerm_14_b_collide + xi_362 + xi_437 + xi_439 + xi_443 + xib_10_collide;
                  const double tmp_b15 = forceTerm_15_b_collide + xi_363 + xi_431 + xi_435 + xi_444 + xib_24_collide;
                  const double tmp_b16 = forceTerm_16_b_collide + xi_365 + xi_430 + xi_436 + xi_444 + xib_26_collide;
                  const double tmp_b17 = forceTerm_17_b_collide + xi_366 + xi_437 + xi_440 + xi_445 + xib_22_collide;
                  const double tmp_b18 = forceTerm_18_b_collide + xi_368 + xi_438 + xi_439 + xi_445 + xib_29_collide;
                  const double xirecolor_0 = tmp_a0 + tmp_b0;
                  const double xirecolor_1 = _data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + _data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_2 = ((1.0) / (xirecolor_1));
                  const double xi_467 = xirecolor_2*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_3 = xirecolor_2*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0];
                  const double xirecolor_4 = tmp_a1 + tmp_b1;
                  const double xirecolor_5 = xi_223;
                  const double xirecolor_6 = xi_221;
                  const double xirecolor_7 = xi_241;
                  const double xi_446 = -xirecolor_7;
                  const double xirecolor_8 = xi_446;
                  const double xirecolor_9 = xi_228;
                  const double xirecolor_10 = -xirecolor_9;
                  const double xirecolor_11 = xi_230;
                  const double xi_447 = xirecolor_10 + xirecolor_11;
                  const double xirecolor_12 = xi_243;
                  const double xi_448 = xirecolor_12 + xirecolor_8;
                  const double xirecolor_13 = xi_447 + xi_448;
                  const double xirecolor_14 = xi_245;
                  const double xirecolor_15 = xi_246;
                  const double xi_449 = xirecolor_14 - xirecolor_15;
                  const double xirecolor_16 = -xi_449;
                  const double xirecolor_17 = xi_215;
                  const double xirecolor_18 = xi_216;
                  const double xirecolor_19 = xi_217;
                  const double xi_450 = xirecolor_18 + xirecolor_19;
                  const double xirecolor_20 = xi_218;
                  const double xirecolor_21 = xi_450 - xirecolor_17 + xirecolor_20;
                  const double xirecolor_22 = xi_214;
                  const double xirecolor_23 = xi_213;
                  const double xi_451 = xirecolor_22 - xirecolor_23;
                  const double xirecolor_24 = -xi_451;
                  const double xirecolor_25 = xi_249;
                  const double xirecolor_26 = -xirecolor_25;
                  const double xirecolor_27 = xi_235;
                  const double xirecolor_28 = xi_234;
                  const double xirecolor_29 = xi_251;
                  const double xi_452 = xirecolor_26 + xirecolor_29;
                  const double xirecolor_30 = xi_452 - xirecolor_27 + xirecolor_28;
                  const double xirecolor_31 = xirecolor_13 + xirecolor_16 + xirecolor_21 + xirecolor_24 + xirecolor_30 - xirecolor_5 + xirecolor_6;
                  const double xirecolor_32 = xi_220;
                  const double xirecolor_33 = xi_222;
                  const double xi_453 = xirecolor_32 + xirecolor_33;
                  const double xirecolor_34 = xi_453 + xirecolor_5 - xirecolor_6;
                  const double xirecolor_35 = xi_225;
                  const double xirecolor_36 = xi_226;
                  const double xi_454 = xirecolor_35 - xirecolor_36;
                  const double xirecolor_37 = -xi_454;
                  const double xirecolor_38 = xi_447;
                  const double xirecolor_39 = xi_232;
                  const double xirecolor_40 = xi_233;
                  const double xi_455 = xirecolor_39 + xirecolor_40;
                  const double xirecolor_41 = xirecolor_27 - xirecolor_28;
                  const double xirecolor_42 = xi_455 + xirecolor_41;
                  const double xirecolor_43 = xirecolor_38 + xirecolor_42;
                  const double xirecolor_44 = xi_451 + xirecolor_21 + xirecolor_34 + xirecolor_37 + xirecolor_43;
                  const double xirecolor_45 = xi_449;
                  const double xirecolor_46 = xi_448;
                  const double xirecolor_47 = xi_257;
                  const double xirecolor_48 = -xirecolor_47;
                  const double xirecolor_49 = xi_259;
                  const double xi_456 = xirecolor_48 + xirecolor_49;
                  const double xirecolor_50 = xi_454 + xi_456;
                  const double xirecolor_51 = xi_450 + xirecolor_17 - xirecolor_20 + xirecolor_24 + xirecolor_34 + xirecolor_45 + xirecolor_46 + xirecolor_50;
                  const double xirecolor_52 = pow((xirecolor_31*xirecolor_31) + (xirecolor_44*xirecolor_44) + (xirecolor_51*xirecolor_51), 0.5);
                  const double xirecolor_53 = ((1.0) / (xirecolor_52));
                  const double xi_457 = xirecolor_31*xirecolor_53;
                  const double xi_458 = xirecolor_51*xirecolor_53;
                  const double xi_459 = xirecolor_44*xirecolor_53;
                  const bool xirecolor_54 = xirecolor_52 > 0.0;
                  const double xirecolor_55 = beta*((1.0) / ((xirecolor_1*xirecolor_1)))*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0]*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0];
                  const double xirecolor_56 = xirecolor_55*(0.055555555555555552*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.055555555555555552*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_57 = xirecolor_56*((xirecolor_54) ? (xi_457): (0.0));
                  const double xirecolor_58 = tmp_a2 + tmp_b2;
                  const double xirecolor_59 = xirecolor_56*((xirecolor_54) ? (-xi_457): (0.0));
                  const double xirecolor_60 = tmp_a3 + tmp_b3;
                  const double xirecolor_61 = xirecolor_56*((xirecolor_54) ? (-xi_458): (0.0));
                  const double xirecolor_62 = tmp_a4 + tmp_b4;
                  const double xirecolor_63 = xirecolor_56*((xirecolor_54) ? (xi_458): (0.0));
                  const double xirecolor_64 = tmp_a5 + tmp_b5;
                  const double xirecolor_65 = xirecolor_56*((xirecolor_54) ? (-xi_459): (0.0));
                  const double xirecolor_66 = tmp_a6 + tmp_b6;
                  const double xirecolor_67 = xirecolor_56*((xirecolor_54) ? (xi_459): (0.0));
                  const double xirecolor_68 = tmp_a7 + tmp_b7;
                  const double xirecolor_69 = xi_313;
                  const double xirecolor_70 = xi_314;
                  const double xirecolor_71 = xi_453;
                  const double xirecolor_72 = xirecolor_50 + xirecolor_71;
                  const double xirecolor_73 = xi_318;
                  const double xirecolor_74 = xirecolor_25 - xirecolor_29 + xirecolor_73;
                  const double xirecolor_75 = xi_315 - xirecolor_11 + xirecolor_41 - xirecolor_69 + xirecolor_70 + xirecolor_72 + xirecolor_74 + xirecolor_9;
                  const double xirecolor_76 = xirecolor_53*0.70710678118654757;
                  const double xi_460 = xirecolor_75*xirecolor_76;
                  const double xirecolor_77 = xirecolor_55*(0.027777777777777776*_data_rho_a[_stride_rho_a_1*ctr_1 + _stride_rho_a_2*ctr_2 + ctr_0] + 0.027777777777777776*_data_rho_b[_stride_rho_b_1*ctr_1 + _stride_rho_b_2*ctr_2 + ctr_0]);
                  const double xirecolor_78 = xirecolor_77*((xirecolor_54) ? (-xi_460): (0.0));
                  const double xirecolor_79 = tmp_a8 + tmp_b8;
                  const double xirecolor_80 = xi_325;
                  const double xirecolor_81 = xi_326;
                  const double xi_461 = -xirecolor_80 + xirecolor_81;
                  const double xirecolor_82 = xi_327;
                  const double xirecolor_83 = xi_328 + xi_461 + xirecolor_30 + xirecolor_38 + xirecolor_72 + xirecolor_82;
                  const double xi_462 = xirecolor_76*xirecolor_83;
                  const double xirecolor_84 = xirecolor_77*((xirecolor_54) ? (xi_462): (0.0));
                  const double xirecolor_85 = tmp_a9 + tmp_b9;
                  const double xirecolor_86 = xirecolor_77*((xirecolor_54) ? (-xi_462): (0.0));
                  const double xirecolor_87 = tmp_a10 + tmp_b10;
                  const double xirecolor_88 = xirecolor_77*((xirecolor_54) ? (xi_460): (0.0));
                  const double xirecolor_89 = tmp_a11 + tmp_b11;
                  const double xirecolor_90 = xi_455 + xirecolor_37 + xirecolor_71;
                  const double xirecolor_91 = -xi_446 - xi_461 - xirecolor_12;
                  const double xirecolor_92 = xi_336 + xirecolor_45 + xirecolor_74 + xirecolor_90 + xirecolor_91;
                  const double xi_463 = xirecolor_76*xirecolor_92;
                  const double xirecolor_93 = xirecolor_77*((xirecolor_54) ? (-xi_463): (0.0));
                  const double xirecolor_94 = tmp_a12 + tmp_b12;
                  const double xirecolor_95 = xirecolor_16 + xirecolor_69 - xirecolor_70;
                  const double xirecolor_96 = xi_346 + xi_452 + xirecolor_46 + xirecolor_82 + xirecolor_90 + xirecolor_95;
                  const double xi_464 = xirecolor_76*xirecolor_96;
                  const double xirecolor_97 = xirecolor_77*((xirecolor_54) ? (-xi_464): (0.0));
                  const double xirecolor_98 = tmp_a13 + tmp_b13;
                  const double xirecolor_99 = xi_354 + xi_456 + xirecolor_13 + xirecolor_42 + xirecolor_45 + xirecolor_73 + xirecolor_82;
                  const double xi_465 = xirecolor_76*xirecolor_99;
                  const double xirecolor_100 = xirecolor_77*((xirecolor_54) ? (-xi_465): (0.0));
                  const double xirecolor_101 = tmp_a14 + tmp_b14;
                  const double xirecolor_102 = xi_360 + xirecolor_43 + xirecolor_47 - xirecolor_49 + xirecolor_91 + xirecolor_95;
                  const double xi_466 = xirecolor_102*xirecolor_76;
                  const double xirecolor_103 = xirecolor_77*((xirecolor_54) ? (-xi_466): (0.0));
                  const double xirecolor_104 = tmp_a15 + tmp_b15;
                  const double xirecolor_105 = xirecolor_77*((xirecolor_54) ? (xi_464): (0.0));
                  const double xirecolor_106 = tmp_a16 + tmp_b16;
                  const double xirecolor_107 = xirecolor_77*((xirecolor_54) ? (xi_463): (0.0));
                  const double xirecolor_108 = tmp_a17 + tmp_b17;
                  const double xirecolor_109 = xirecolor_77*((xirecolor_54) ? (xi_466): (0.0));
                  const double xirecolor_110 = tmp_a18 + tmp_b18;
                  const double xirecolor_111 = xirecolor_77*((xirecolor_54) ? (xi_465): (0.0));
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + ctr_0] = xirecolor_0*xirecolor_3;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + _stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_4 + xirecolor_57;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 2*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_58 + xirecolor_59;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 3*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_60 + xirecolor_61;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 4*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_62 + xirecolor_63;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 5*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_64 + xirecolor_65;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 6*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_66 + xirecolor_67;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 7*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_68 + xirecolor_78;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 8*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_79 + xirecolor_84;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 9*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_85 + xirecolor_86;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 10*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_87 + xirecolor_88;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 11*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_89 + xirecolor_93;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 12*_stride_pdfs_a_3 + ctr_0] = xirecolor_3*xirecolor_94 + xirecolor_97;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 13*_stride_pdfs_a_3 + ctr_0] = xirecolor_100 + xirecolor_3*xirecolor_98;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 14*_stride_pdfs_a_3 + ctr_0] = xirecolor_101*xirecolor_3 + xirecolor_103;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 15*_stride_pdfs_a_3 + ctr_0] = xirecolor_104*xirecolor_3 + xirecolor_105;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 16*_stride_pdfs_a_3 + ctr_0] = xirecolor_106*xirecolor_3 + xirecolor_107;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 17*_stride_pdfs_a_3 + ctr_0] = xirecolor_108*xirecolor_3 + xirecolor_109;
                  _data_pdfs_a[_stride_pdfs_a_1*ctr_1 + _stride_pdfs_a_2*ctr_2 + 18*_stride_pdfs_a_3 + ctr_0] = xirecolor_110*xirecolor_3 + xirecolor_111;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + ctr_0] = xi_467*xirecolor_0;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + _stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_4 - xirecolor_57;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 2*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_58 - xirecolor_59;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 3*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_60 - xirecolor_61;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 4*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_62 - xirecolor_63;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 5*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_64 - xirecolor_65;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 6*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_66 - xirecolor_67;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 7*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_68 - xirecolor_78;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 8*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_79 - xirecolor_84;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 9*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_85 - xirecolor_86;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 10*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_87 - xirecolor_88;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 11*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_89 - xirecolor_93;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 12*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_94 - xirecolor_97;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 13*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_98 - xirecolor_100;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 14*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_101 - xirecolor_103;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 15*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_104 - xirecolor_105;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 16*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_106 - xirecolor_107;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 17*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_108 - xirecolor_109;
                  _data_pdfs_b[_stride_pdfs_b_1*ctr_1 + _stride_pdfs_b_2*ctr_2 + 18*_stride_pdfs_b_3 + ctr_0] = xi_467*xirecolor_110 - xirecolor_111;
               }
            }
         }
      }
   }
}
}


void ColorGradientCollideSweepDoublePrecisionAVX::run(IBlock * block)
{
   
    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(-1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(-1, -1, -1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(force_a->xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(force_a->xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(force_a->ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(force_a->ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(force_a->zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(force_a->zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}


void ColorGradientCollideSweepDoublePrecisionAVX::runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block)
{
   
    CellInterval ci = globalCellInterval;
    CellInterval blockBB = blocks->getBlockCellBB( *block);
    blockBB.expand( ghostLayers );
    ci.intersect( blockBB );
    blocks->transformGlobalToBlockLocalCellInterval( ci, *block );
    if( ci.empty() )
        return;

    auto pdfs_a = block->getData< field::GhostLayerField<double, 19> >(pdfs_aID);
    auto pdfs_b = block->getData< field::GhostLayerField<double, 19> >(pdfs_bID);
    auto phasefield = block->getData< field::GhostLayerField<double, 1> >(phasefieldID);
    auto rho_a = block->getData< field::GhostLayerField<double, 1> >(rho_aID);
    auto force_a = block->getData< field::GhostLayerField<double, 3> >(force_aID);
    auto velocity = block->getData< field::GhostLayerField<double, 3> >(velocityID);
    auto force_b = block->getData< field::GhostLayerField<double, 3> >(force_bID);
    auto rho_b = block->getData< field::GhostLayerField<double, 1> >(rho_bID);

    auto & omega_shear_a = this->omega_shear_a_;
    auto & beta = this->beta_;
    auto & sigma = this->sigma_;
    auto & omega_shear_b = this->omega_shear_b_;
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_a->nrOfGhostLayers()))
    double * RESTRICT const _data_force_a = force_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(force_b->nrOfGhostLayers()))
    double * RESTRICT const _data_force_b = force_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(force_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_a->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_a = pdfs_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(pdfs_b->nrOfGhostLayers()))
    double * RESTRICT  _data_pdfs_b = pdfs_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(pdfs_b->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) pdfs_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(phasefield->nrOfGhostLayers()))
    double * RESTRICT const _data_phasefield = phasefield->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) phasefield->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_a->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_a = rho_a->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(rho_b->nrOfGhostLayers()))
    double * RESTRICT const _data_rho_b = rho_b->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL((uintptr_t) rho_b->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(ci.xMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.yMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    WALBERLA_ASSERT_GREATER_EQUAL(ci.zMin() - 1, -int_c(velocity->nrOfGhostLayers()))
    double * RESTRICT const _data_velocity = velocity->dataAt(ci.xMin() - 1, ci.yMin() - 1, ci.zMin() - 1, 0);
    WALBERLA_ASSERT_EQUAL(velocity->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) velocity->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->xSizeWithGhostLayer(), int64_t(int64_c(ci.xSize()) + 2))
    const int64_t _size_force_a_0 = int64_t(int64_c(ci.xSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->ySizeWithGhostLayer(), int64_t(int64_c(ci.ySize()) + 2))
    const int64_t _size_force_a_1 = int64_t(int64_c(ci.ySize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    WALBERLA_ASSERT_GREATER_EQUAL(force_a->zSizeWithGhostLayer(), int64_t(int64_c(ci.zSize()) + 2))
    const int64_t _size_force_a_2 = int64_t(int64_c(ci.zSize()) + 2);
    WALBERLA_ASSERT_EQUAL(force_a->layout(), field::fzyx)
    WALBERLA_ASSERT_EQUAL((uintptr_t) force_a->dataAt(0, 0, 0, 0) %32, 0)
    const int64_t _stride_force_a_1 = int64_t(force_a->yStride());
    const int64_t _stride_force_a_2 = int64_t(force_a->zStride());
    const int64_t _stride_force_a_3 = int64_t(1 * int64_t(force_a->fStride()));
    const int64_t _stride_force_b_1 = int64_t(force_b->yStride());
    const int64_t _stride_force_b_2 = int64_t(force_b->zStride());
    const int64_t _stride_force_b_3 = int64_t(1 * int64_t(force_b->fStride()));
    const int64_t _stride_pdfs_a_1 = int64_t(pdfs_a->yStride());
    const int64_t _stride_pdfs_a_2 = int64_t(pdfs_a->zStride());
    const int64_t _stride_pdfs_a_3 = int64_t(1 * int64_t(pdfs_a->fStride()));
    const int64_t _stride_pdfs_b_1 = int64_t(pdfs_b->yStride());
    const int64_t _stride_pdfs_b_2 = int64_t(pdfs_b->zStride());
    const int64_t _stride_pdfs_b_3 = int64_t(1 * int64_t(pdfs_b->fStride()));
    const int64_t _stride_phasefield_1 = int64_t(phasefield->yStride());
    const int64_t _stride_phasefield_2 = int64_t(phasefield->zStride());
    const int64_t _stride_rho_a_1 = int64_t(rho_a->yStride());
    const int64_t _stride_rho_a_2 = int64_t(rho_a->zStride());
    const int64_t _stride_rho_b_1 = int64_t(rho_b->yStride());
    const int64_t _stride_rho_b_2 = int64_t(rho_b->zStride());
    const int64_t _stride_velocity_1 = int64_t(velocity->yStride());
    const int64_t _stride_velocity_2 = int64_t(velocity->zStride());
    const int64_t _stride_velocity_3 = int64_t(1 * int64_t(velocity->fStride()));
    internal_6ed8e43fbb9656b49ba1d2a7ec0e653c::colorgradientcollidesweepdoubleprecisionavx_colorgradientcollidesweepdoubleprecisionavx(_data_force_a, _data_force_b, _data_pdfs_a, _data_pdfs_b, _data_phasefield, _data_rho_a, _data_rho_b, _data_velocity, _size_force_a_0, _size_force_a_1, _size_force_a_2, _stride_force_a_1, _stride_force_a_2, _stride_force_a_3, _stride_force_b_1, _stride_force_b_2, _stride_force_b_3, _stride_pdfs_a_1, _stride_pdfs_a_2, _stride_pdfs_a_3, _stride_pdfs_b_1, _stride_pdfs_b_2, _stride_pdfs_b_3, _stride_phasefield_1, _stride_phasefield_2, _stride_rho_a_1, _stride_rho_a_2, _stride_rho_b_1, _stride_rho_b_2, _stride_velocity_1, _stride_velocity_2, _stride_velocity_3, beta, omega_shear_a, omega_shear_b, sigma);
    
}



} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_INTEL )
#pragma warning pop
#endif
