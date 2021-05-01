// cpd - Coherent Point Drift
// Copyright (C) 2017 Pete Gadomski <pete.gadomski@gmail.com>
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

/// \file
///
/// Generic coherent point drift transform.
///
/// Downstreams shouldn't need to include this file directly, use a realization
/// of a transform (e.g. `Rigid`) instead.

#pragma once

#include <chrono>
#include "cpd_gauss_transform.h"
#include "cpd_matrix.h"
#include "cpd_normalization.h"
#include "cpd_utils.h"
#include <memory>
#include <iostream>
#include <algorithm>
#include <fstream>

namespace cpd {

/// The default number of iterations allowed.
const size_t DEFAULT_MAX_ITERATIONS = 150;
/// Whether points should be normalized by default.
const bool DEFAULT_NORMALIZE = false;
/// The default outlier weight.
const double DEFAULT_OUTLIERS = 0.1;
/// The default tolerance.
const double DEFAULT_TOLERANCE = 1e-8;
/// The default sigma2.
const double DEFAULT_SIGMA2 = 0.0;
/// Whether correspondence vector should be computed by default.
const bool DEFAULT_CORRESPONDENCE = false;
/// Are the scalings of the two datasets linked by default?
const bool DEFAULT_LINKED = true;

/// The result of a generic transform run.
struct Result {
    /// The final moved points.
    Matrix points;
    /// The final sigma2 value.
    double sigma2;
    /// The correspondence vector.
    IndexVector correspondence;
    /// The runtime.
    std::chrono::microseconds runtime;
    /// The number of iterations.
    size_t iterations;

    // The transformation used to produce result
    Matrix tform;

    /// De-normalize this result.
    ///
    /// Generally, this scales the points back, and sometimes adjust transforms
    /// or shifts or the like.
    virtual void denormalize(const Normalization& normalization);
};

/// Generic coherent point drift transform.
///
/// An abstract base class for real transforms, e.g. `Rigid` or `Nonrigid`.
template <typename Result>
class Transform {
public:
    Transform()
      : m_correspondence(DEFAULT_CORRESPONDENCE)
      , m_gauss_transform(GaussTransform::make_default())
      , m_max_iterations(DEFAULT_MAX_ITERATIONS)
      , m_normalize(DEFAULT_NORMALIZE)
      , m_outliers(DEFAULT_OUTLIERS)
      , m_sigma2(DEFAULT_SIGMA2)
      , m_tolerance(DEFAULT_TOLERANCE) {}

    inline void set_mesh_num(const int& mesh_num){
      m_mesh_num=mesh_num;
    }

    /// Sets whether the correspondence vector will be computed.
    Transform& correspondence(bool correspondence) {
        m_correspondence = correspondence;
        return *this;
    }

    /// Sets the gauss transform.
    Transform& gauss_transform(
        std::unique_ptr<GaussTransform> gauss_transform) {
        m_gauss_transform = std::move(gauss_transform);
        return *this;
    }

    /// Sets the max iterations for this transform.
    Transform& max_iterations(double max_iterations) {
        m_max_iterations = max_iterations;
        return *this;
    }

    /// Sets whether to normalize the points before running cpd.
    Transform& normalize(bool normalize) {
        m_normalize = normalize;
        return *this;
    }

    /// Sets the outlier tolerance.
    Transform& outliers(double outliers) {
        m_outliers = outliers;
        return *this;
    }

    /// Sets the sigma2 value for this transform.
    Transform& sigma2(double sigma2) {
        m_sigma2 = sigma2;
        return *this;
    }

    /// Sets the final tolerance.
    Transform& tolerance(double tolerance) {
        m_tolerance = tolerance;
        return *this;
    }

    /// Sets the tracked mesh number.
    Transform& mesh_num(int mesh_num) {
        m_mesh_num = mesh_num;
        return *this;
    }

    /// Runs this transform for the provided matrices.
    Result run(Matrix fixed, Matrix moving) {
        auto tic = std::chrono::high_resolution_clock::now();
        Normalization normalization(fixed, moving, linked());
        if (m_normalize) {
            fixed = normalization.fixed;
            moving = normalization.moving;
        }
        this->init(fixed, moving);
        Result result;
        result.points = moving;

        // Sigma2 : Isotropic Covariance Matrix
        if (m_sigma2 == 0.0) {
            result.sigma2 = cpd::default_sigma2(fixed, moving);
        } else if (m_normalize) {
            result.sigma2 = m_sigma2 / normalization.fixed_scale;
        } else {
            result.sigma2 = m_sigma2;
        }

        int iter = 0;
        double ntol = m_tolerance + 10.0;
        double l = 0.;

         while (iter < m_max_iterations && ntol > m_tolerance &&
        // while (iter < 1 && ntol > m_tolerance &&
               result.sigma2 > 10 * std::numeric_limits<double>::epsilon()) {
            Probabilities probabilities = m_gauss_transform->compute(
                fixed, result.points, result.sigma2, m_outliers);

            visualize_probabilities(probabilities,moving,m_mesh_num,iter);
            this->modify_probabilities(probabilities);

            ntol = std::abs((probabilities.l - l) / probabilities.l);
            l = probabilities.l;
            result =
                this->compute_one(fixed, moving, probabilities, result.sigma2);
            ++iter;
        }

        if (m_normalize) {
            result.denormalize(normalization);
        }

        if (m_correspondence) {
            GaussTransformDirect direct;
            Probabilities probabilities = direct.compute(
                fixed, result.points, result.sigma2, m_outliers);
            result.correspondence = probabilities.correspondence;
            assert(result.correspondence.rows() > 0);
        }

        auto toc = std::chrono::high_resolution_clock::now();
        result.runtime =
            std::chrono::duration_cast<std::chrono::microseconds>(toc - tic);
        result.iterations = iter;

        return result;
    }

    /// Initialize this transform for the provided matrices.
    ///
    /// This is called before beginning each run, but after normalization. In
    /// general, transforms shouldn't need to be initialized, but the nonrigid
    /// does.
    virtual void init(const Matrix& fixed, const Matrix& moving) {
    }

    /// Modifies `Probabilities` in some way.
    ///
    /// Some types of transform need to tweak the probabilities before moving on
    /// with an interation. The default behavior is to do nothing.
    virtual void modify_probabilities(Probabilities& probabilities) const {}

    /// Computes one iteration of the transform.
    virtual Result compute_one(const Matrix& fixed, const Matrix& moving,
                               const Probabilities& probabilities,
                               double sigma2) const = 0;

    /// Returns true if the normalization should be linked.
    ///
    /// No effect if no normalization is applied.
    virtual bool linked() const = 0;

    // used in probability visualization
    static void val2rgb(float min, float max, float val, int (&rgb)[3]){
      float ratio = 2* (val-min)/(max-min);
      rgb[0] = int(std::max(0.0f,255*(1-ratio)));
      rgb[2] = int(std::max(0.0f,255*(ratio-1)));
      rgb[1] = 255 - rgb[2] - rgb[0];
    }

    // TODO: Modify to accept only prob matrix that we're interested in and
    //  perform a check to make sure fixed.size == prob_mat.size
    static void visualize_probabilities(const Probabilities& probabilities,
      const Matrix& fixed, const int& mesh_num, const int& iteration){

      float p1_max =probabilities.p1(0);
      float p1_min =probabilities.p1(0);
      for (size_t index = 0; index < probabilities.p1.rows(); index++) {
        if(probabilities.p1(index)>p1_max) p1_max = probabilities.p1(index);
        if(probabilities.p1(index)<p1_min) p1_min = probabilities.p1(index);
      }

    //   //Output target with indicated probabilities
    //   std::ofstream target_p_out;
    //   std::string outname = "/MeshTracker/resource/tests/t_"+
    //                         std::to_string(mesh_num)+"_"+
    //                         std::to_string(iteration)+"_p.ply";
    //   target_p_out.open(outname);
    //   target_p_out<<"ply\n";
    //   target_p_out<<"format ascii 1.0\n";
    //   target_p_out<<"element vertex "<<fixed.rows()<<"\n";
    //   target_p_out<<"property float x\n";
    //   target_p_out<<"property float y\n";
    //   target_p_out<<"property float z\n";
    //   target_p_out<<"property uchar red\n";
    //   target_p_out<<"property uchar green\n";
    //   target_p_out<<"property uchar blue\n";
    //   target_p_out<<"end_header\n";

    //   for (size_t row = 0; row < fixed.rows(); row++) {
    //     int rgb[3] = {0,0,0};
    //     val2rgb(p1_min,p1_max,probabilities.p1(row),rgb);
    //     target_p_out<<" "<<fixed(row,0)<<" "
    //                 <<fixed(row,1)<<" "
    //                 <<fixed(row,2)<<" "
    //                 <<rgb[0]<<" "
    //                 <<rgb[1]<<" "
    //                 <<rgb[2]<<" "<<std::endl;
    //   }
    //   target_p_out.close();
    }



private:
    bool m_correspondence;
    std::unique_ptr<GaussTransform> m_gauss_transform;
    size_t m_max_iterations;
    bool m_normalize;
    double m_outliers;
    double m_sigma2;
    double m_tolerance;
    int m_mesh_num=0;
};
} // namespace cpd
