// FlatField.h
// Contiguous-memory replacements for ADMB d3_array / d4_array for ocean
// forcing data (np1, sst, vld, un, vn, tempn, oxygen, forage).
//
// Motivation
// ----------
// ADMB's d3_array/d4_array are ragged: each row i has its own j range
// (jinf[i]..jsup[i]), so the data is a chain of pointer-of-pointers.
// FlatField3D and FlatField4D pack valid ocean cells into a single
// contiguous double[] buffer, which:
//   1. Removes three levels of pointer indirection on every access.
//   2. Allows the buffer to be aliased directly to a DataProvider
//      shared-memory window (Step 3 of the MPI memory-reduction plan),
//      eliminating all serialisation/deserialisation copies.
//
// Layout
// ------
//   FlatField3D  buf[(t - t0) * mapCells + cell]
//   FlatField4D  buf[((t - t0) * K + k) * mapCells + cell]
//
//   cell = irow_base[i - imin] + (j - jinf[i])
//   irow_base[i - imin] = number of valid cells in rows imin .. i-1
//
// Operator interface
// ------------------
// Both classes expose the same [t][i][j] / [t][k][i][j] bracket syntax
// as the ADMB originals through lightweight proxy structs (RawSlice2D,
// Slice3D).  No call-site changes are needed for scalar reads/writes.
// Whole-slice operations (assignment, scaling, dmatrix conversion) are
// provided on the proxy types to cover the handful of non-scalar uses in
// ReadAll and the physics functions.
//
// DataProvider aliasing
// ---------------------
// Call  allocate(map, t0, nbt, externalBuf)  with the DataProvider
// window pointer so that the field aliases that buffer directly (owns_
// stays false, no allocation, no copy ever needed).

#ifndef __FlatField_h__
#define __FlatField_h__

#include <vector>
#include <cstring>
#include <cstddef>
#include <fvar.hpp>   // dmatrix, ivector
#include "Map.h"      // PMap

// ---------------------------------------------------------------------------
// RawSlice2D
// A lightweight non-owning view of one 2D (i,j) slice stored in the flat
// buffer.  Returned by FlatField3D::operator[](t) and by
// FlatField4D::Slice3D::operator[](k).
// ---------------------------------------------------------------------------
struct RawSlice2D {
    double*     ptr;        // pointer to mapCells doubles in cell-major order
    std::size_t mapCells;
    int         imin, imax;
    const int*  jinf;       // jinf[i - imin]  (lower j bound for row i)
    const int*  jsup;       // jsup[i - imin]  (upper j bound for row i)
    const int*  irow_base;  // irow_base[i - imin] = cell offset for row i

    // --- element access --------------------------------------------------

    // Inner proxy for the final [j] subscript
    struct RowProxy {
        double* row;   // ptr + irow_base[i-imin]
        int     j0;    // jinf[i-imin]
        double& operator[](int j) const { return row[j - j0]; }
    };

    RowProxy operator[](int i) const {
        return {ptr + irow_base[i - imin], jinf[i - imin]};
    }
    double& operator()(int i, int j) const {
        return ptr[irow_base[i - imin] + (j - jinf[i - imin])];
    }

    // --- slice-level assignment / arithmetic ----------------------------

    // Fill all cells with a scalar  (mat.vld[t] = 1.0)
    RawSlice2D& operator=(double val) {
        for (std::size_t c = 0; c < mapCells; ++c) ptr[c] = val;
        return *this;
    }
    // Copy from another same-layout slice  (mat.sst[t] = mat.tempn[t][0])
    RawSlice2D& operator=(const RawSlice2D& o) {
        std::memcpy(ptr, o.ptr, mapCells * sizeof(double));
        return *this;
    }
    // Copy from a dmatrix  (used in commit() and direct assignments)
    RawSlice2D& operator=(const dmatrix& m);

    // Scale in place  (mat.vld[t] /= 1000)
    void scaleInPlace(double factor) {
        for (std::size_t c = 0; c < mapCells; ++c) ptr[c] *= factor;
    }
    // Add a scalar to every cell  (forage += 0.000001)
    void addInPlace(double val) {
        for (std::size_t c = 0; c < mapCells; ++c) ptr[c] += val;
    }
    // Add a ragged dmatrix element-wise into this slice  (forage[t][n] += mats[n])
    RawSlice2D& operator+=(const dmatrix& m) {
        for (int i = imin; i <= imax; ++i) {
            const int jlo = jinf[i - imin];
            const int jhi = jsup[i - imin];
            double* row = ptr + irow_base[i - imin];
            for (int j = jlo; j <= jhi; ++j) row[j - jlo] += m[i][j];
        }
        return *this;
    }

    // Materialise a copy as a ragged dmatrix
    // (enables implicit conversion for  const dmatrix&  parameters and
    //  for  mat.u = mat.un[t][k]  assignments)
    operator dmatrix() const;

    // Return a scaled copy as dmatrix  (elarvae_dt * mat.un[t][k])
    dmatrix operator*(double scale) const;
};

// Free operator so  scalar * RawSlice2D  also works
inline dmatrix operator*(double scale, const RawSlice2D& s) { return s * scale; }


// ---------------------------------------------------------------------------
// FlatField3D  —  contiguous [t][i][j] forcing field
// Replaces d3_array for np1, sst, ph1, vld.
// ---------------------------------------------------------------------------
class FlatField3D {
public:
    FlatField3D()  = default;
    ~FlatField3D() { if (owns_ && data_) delete[] data_; }

    FlatField3D(const FlatField3D&)            = delete;
    FlatField3D& operator=(const FlatField3D&) = delete;

    // Move constructor: transfers ownership of the flat buffer.
    FlatField3D(FlatField3D&& o) noexcept
        : data_(o.data_), owns_(o.owns_), t0_(o.t0_), T_(o.T_),
          mapCells_(o.mapCells_), imin_(o.imin_), imax_(o.imax_),
          jinf_(std::move(o.jinf_)), jsup_(std::move(o.jsup_)),
          irow_base_(std::move(o.irow_base_)), stage_(o.stage_)
    { o.data_ = nullptr; o.owns_ = false; }
    FlatField3D& operator=(FlatField3D&&) = delete;

    // Allocate own flat buffer (externalBuf == nullptr → new double[])
    // or alias an external buffer (DataProvider shared-memory window).
    void allocate(const PMap& map, int t0, int nbt,
                  double* externalBuf = nullptr);

    void initialize() {
        if (data_) std::memset(data_, 0, T_ * mapCells_ * sizeof(double));
    }

    // --- scalar access ---------------------------------------------------
    double& operator()(int t, int i, int j) {
        return data_[off3(t, i, j)];
    }
    double operator()(int t, int i, int j) const {
        return data_[off3(t, i, j)];
    }

    // --- 2D slice proxy --------------------------------------------------
    RawSlice2D operator[](int t) {
        return makeSlice(data_ + (std::size_t)(t - t0_) * mapCells_);
    }
    RawSlice2D operator()(int t) { return (*this)[t]; }
    // const version (ptr cast away const — forcing fields are write-once)
    RawSlice2D operator[](int t) const {
        return makeSlice(const_cast<double*>(data_)
                         + (std::size_t)(t - t0_) * mapCells_);
    }

    // --- staging for rbin_input2d ----------------------------------------
    // Returns a pre-allocated ragged dmatrix for rbin_input2d to fill.
    dmatrix& stage() { return stage_; }
    // Pack stage_ into the flat buffer at time step t.
    void commit(int t);

    // --- bulk operations -------------------------------------------------
    // Scale every value in time step t by factor.
    void scaleTimeStep(int t, double factor) {
        double* p = data_ + (std::size_t)(t - t0_) * mapCells_;
        for (std::size_t c = 0; c < mapCells_; ++c) p[c] *= factor;
    }

private:
    double*          data_     = nullptr;
    bool             owns_     = false;
    int              t0_       = 0;
    std::size_t      T_        = 0;
    std::size_t      mapCells_ = 0;
    int              imin_     = 0, imax_ = 0;
    std::vector<int> jinf_;       // [i - imin]
    std::vector<int> jsup_;       // [i - imin]
    std::vector<int> irow_base_;  // [i - imin]
    dmatrix          stage_;      // staging dmatrix for rbin_input2d

    std::size_t off3(int t, int i, int j) const {
        return (std::size_t)(t - t0_) * mapCells_
             + irow_base_[i - imin_]
             + (j - jinf_[i - imin_]);
    }
    RawSlice2D makeSlice(double* p) const {
        return {p, mapCells_, imin_, imax_,
                jinf_.data(), jsup_.data(), irow_base_.data()};
    }
    void buildIndex(const PMap& map, int t0, int nbt);
};


// ---------------------------------------------------------------------------
// FlatField4D  —  contiguous [t][k][i][j] forcing field
// Replaces d4_array for un, vn, tempn, oxygen, forage.
// K = nb_layer (currents / temperature / oxygen) or nb_forage.
// ---------------------------------------------------------------------------
class FlatField4D {
public:
    FlatField4D()  = default;
    ~FlatField4D() { if (owns_ && data_) delete[] data_; }

    FlatField4D(const FlatField4D&)            = delete;
    FlatField4D& operator=(const FlatField4D&) = delete;

    // Move constructor: transfers ownership of the flat buffer.
    FlatField4D(FlatField4D&& o) noexcept
        : data_(o.data_), owns_(o.owns_), t0_(o.t0_), T_(o.T_), K_(o.K_),
          mapCells_(o.mapCells_), imin_(o.imin_), imax_(o.imax_),
          jinf_(std::move(o.jinf_)), jsup_(std::move(o.jsup_)),
          irow_base_(std::move(o.irow_base_)), stage_(o.stage_)
    { o.data_ = nullptr; o.owns_ = false; }
    FlatField4D& operator=(FlatField4D&&) = delete;

    // allocate(map, t0, nbt, K [, externalBuf])
    void allocate(const PMap& map, int t0, int nbt, int K,
                  double* externalBuf = nullptr);

    void initialize() {
        if (data_) std::memset(data_, 0, T_ * K_ * mapCells_ * sizeof(double));
    }

    // --- scalar access ---------------------------------------------------
    double& operator()(int t, int k, int i, int j) {
        return data_[off4(t, k, i, j)];
    }
    double operator()(int t, int k, int i, int j) const {
        return data_[off4(t, k, i, j)];
    }

    // --- 3D slice proxy (one time step, all K layers) --------------------
    struct Slice3D {
        double*     ptr;        // data_ + (t - t0) * K * mapCells
        std::size_t K, mapCells;
        int         imin, imax;
        const int*  jinf, *jsup, *irow_base;

        RawSlice2D operator[](int k) const {
            return {ptr + (std::size_t)k * mapCells,
                    mapCells, imin, imax, jinf, jsup, irow_base};
        }
        // Scale all K layers in this time step in place.
        void scaleInPlace(double factor) {
            double* p = ptr;
            for (std::size_t c = 0; c < K * mapCells; ++c) p[c] *= factor;
        }
    };

    Slice3D operator[](int t) { return makeSlice3(t); }
    Slice3D operator()(int t) { return makeSlice3(t); }
    Slice3D operator[](int t) const { return makeSlice3(t); }

    // 2-arg form: mat.oxygen(t, k) → RawSlice2D  (replaces ADMB d4_array(t,k) → dmatrix&)
    RawSlice2D operator()(int t, int k) { return makeSlice3(t)[k]; }
    RawSlice2D operator()(int t, int k) const { return makeSlice3(t)[k]; }

    // --- staging for rbin_input2d ----------------------------------------
    dmatrix& stage() { return stage_; }
    // Pack stage_ into the flat buffer at (t, k).
    void commit(int t, int k);

    // --- bulk operations -------------------------------------------------
    void scaleTimeStep(int t, double factor) {
        double* p = data_ + (std::size_t)(t - t0_) * K_ * mapCells_;
        for (std::size_t c = 0; c < K_ * mapCells_; ++c) p[c] *= factor;
    }

private:
    double*          data_     = nullptr;
    bool             owns_     = false;
    int              t0_       = 0;
    std::size_t      T_        = 0;
    std::size_t      K_        = 0;
    std::size_t      mapCells_ = 0;
    int              imin_     = 0, imax_ = 0;
    std::vector<int> jinf_;
    std::vector<int> jsup_;
    std::vector<int> irow_base_;
    dmatrix          stage_;

    std::size_t off4(int t, int k, int i, int j) const {
        return ((std::size_t)(t - t0_) * K_ + k) * mapCells_
             + irow_base_[i - imin_]
             + (j - jinf_[i - imin_]);
    }
    Slice3D makeSlice3(int t) const {
        return {const_cast<double*>(data_)
                    + (std::size_t)(t - t0_) * K_ * mapCells_,
                K_, mapCells_, imin_, imax_,
                jinf_.data(), jsup_.data(), irow_base_.data()};
    }
    void buildIndex(const PMap& map, int t0, int nbt, int K);
};

#endif // __FlatField_h__
