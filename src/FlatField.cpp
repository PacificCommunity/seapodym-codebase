// FlatField.cpp
// Implementation of RawSlice2D, FlatField3D, and FlatField4D.
// See FlatField.h for design notes.

#include "FlatField.h"
#include <stdexcept>

// ---------------------------------------------------------------------------
// RawSlice2D — dmatrix interop
// ---------------------------------------------------------------------------

RawSlice2D& RawSlice2D::operator=(const dmatrix& m)
{
    for (int i = imin; i <= imax; ++i) {
        const int jlo = jinf[i - imin];
        const int jhi = jsup[i - imin];
        double* row = ptr + irow_base[i - imin];
        for (int j = jlo; j <= jhi; ++j) {
            row[j - jlo] = m[i][j];
        }
    }
    return *this;
}

RawSlice2D::operator dmatrix() const
{
    // Build index vectors with explicit allocate() — ADMB ivector has
    // no (lb,ub) constructor in all versions.
    ivector jlo_vec, jhi_vec;
    jlo_vec.allocate(imin, imax);
    jhi_vec.allocate(imin, imax);
    for (int i = imin; i <= imax; ++i) {
        jlo_vec[i] = jinf[i - imin];
        jhi_vec[i] = jsup[i - imin];
    }
    dmatrix m;
    m.allocate(imin, imax, jlo_vec, jhi_vec);
    for (int i = imin; i <= imax; ++i) {
        const int jlo = jinf[i - imin];
        const int jhi = jsup[i - imin];
        const double* row = ptr + irow_base[i - imin];
        for (int j = jlo; j <= jhi; ++j) {
            m[i][j] = row[j - jlo];
        }
    }
    return m;
}

dmatrix RawSlice2D::operator*(double scale) const
{
    dmatrix m = static_cast<dmatrix>(*this);
    for (int i = imin; i <= imax; ++i) {
        const int jlo = jinf[i - imin];
        const int jhi = jsup[i - imin];
        for (int j = jlo; j <= jhi; ++j) {
            m[i][j] *= scale;
        }
    }
    return m;
}


// ---------------------------------------------------------------------------
// FlatField3D
// ---------------------------------------------------------------------------

void FlatField3D::buildIndex(const PMap& map, int t0, int nbt)
{
    t0_   = t0;
    T_    = static_cast<std::size_t>(nbt - t0 + 1);
    imin_ = map.imin;
    imax_ = map.imax;

    int ni = imax_ - imin_ + 1;
    jinf_.resize(ni);
    jsup_.resize(ni);
    irow_base_.resize(ni);

    std::size_t cell = 0;
    for (int i = imin_; i <= imax_; ++i) {
        int idx = i - imin_;
        jinf_[idx]      = map.jinf[i];
        jsup_[idx]      = map.jsup[i];
        irow_base_[idx] = static_cast<int>(cell);
        int ncol = map.jsup[i] - map.jinf[i] + 1;
        if (ncol > 0) cell += static_cast<std::size_t>(ncol);
    }
    mapCells_ = cell;

    // Build staging dmatrix using the map's existing ivectors directly —
    // they are already indexed [imin..imax], matching our range.
    stage_.allocate(imin_, imax_, map.jinf, map.jsup);
}

void FlatField3D::allocate(const PMap& map, int t0, int nbt,
                           double* externalBuf)
{
    buildIndex(map, t0, nbt);

    if (externalBuf) {
        data_ = externalBuf;
        owns_ = false;
    } else {
        std::size_t n = T_ * mapCells_;
        data_ = new double[n];
        owns_ = true;
    }
}

void FlatField3D::commit(int t)
{
    double* dst = data_ + static_cast<std::size_t>(t - t0_) * mapCells_;
    for (int i = imin_; i <= imax_; ++i) {
        int idx = i - imin_;
        const int jlo = jinf_[idx];
        const int jhi = jsup_[idx];
        double* row = dst + irow_base_[idx];
        for (int j = jlo; j <= jhi; ++j) {
            row[j - jlo] = stage_[i][j];
        }
    }
}


// ---------------------------------------------------------------------------
// FlatField4D
// ---------------------------------------------------------------------------

void FlatField4D::buildIndex(const PMap& map, int t0, int nbt, int K)
{
    t0_ = t0;
    T_  = static_cast<std::size_t>(nbt - t0 + 1);
    K_  = static_cast<std::size_t>(K);
    imin_ = map.imin;
    imax_ = map.imax;

    int ni = imax_ - imin_ + 1;
    jinf_.resize(ni);
    jsup_.resize(ni);
    irow_base_.resize(ni);

    std::size_t cell = 0;
    for (int i = imin_; i <= imax_; ++i) {
        int idx = i - imin_;
        jinf_[idx]      = map.jinf[i];
        jsup_[idx]      = map.jsup[i];
        irow_base_[idx] = static_cast<int>(cell);
        int ncol = map.jsup[i] - map.jinf[i] + 1;
        if (ncol > 0) cell += static_cast<std::size_t>(ncol);
    }
    mapCells_ = cell;

    stage_.allocate(imin_, imax_, map.jinf, map.jsup);
}

void FlatField4D::allocate(const PMap& map, int t0, int nbt, int K,
                           double* externalBuf)
{
    buildIndex(map, t0, nbt, K);

    if (externalBuf) {
        data_ = externalBuf;
        owns_ = false;
    } else {
        std::size_t n = T_ * K_ * mapCells_;
        data_ = new double[n];
        owns_ = true;
    }
}

void FlatField4D::commit(int t, int k)
{
    double* dst = data_
                + (static_cast<std::size_t>(t - t0_) * K_ + k) * mapCells_;
    for (int i = imin_; i <= imax_; ++i) {
        int idx = i - imin_;
        const int jlo = jinf_[idx];
        const int jhi = jsup_[idx];
        double* row = dst + irow_base_[idx];
        for (int j = jlo; j <= jhi; ++j) {
            row[j - jlo] = stage_[i][j];
        }
    }
}
