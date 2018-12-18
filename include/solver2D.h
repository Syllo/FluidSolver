/*
 * Copyright (c) 2018 Maxime Schmitt <max.schmitt@unistra.fr>
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors
 * may be used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef SOLVER2D_H_
#define SOLVER2D_H_

#include <stddef.h>

void dens_step2D(size_t M, size_t N, float d[M][N], float d0[M][N],
                 float u[M][N], float v[M][N], float diff, float dt);

void vel_step2D(size_t M, size_t N, float u[M][N], float v[M][N],
                float u0[M][N], float v0[M][N], float visc, float dt);

void add_force_source2D(
    size_t sim_step, size_t M, size_t N, float d[restrict M][N],
    float u[restrict M][N], float v[restrict M][N], size_t num_dens,
    float dens_val[restrict num_dens], size_t dens_sourceX[restrict num_dens],
    size_t dens_sourceY[restrict num_dens],
    size_t dens_frequency[restrict num_dens], size_t num_vel,
    float vel_force[restrict num_vel], size_t vel_sourceX[restrict num_vel],
    size_t vel_sourceY[restrict num_vel], float vel_angleXY[restrict num_vel],
    size_t vel_frequency[restrict num_vel]);

#endif // SOLVER2D_H_
