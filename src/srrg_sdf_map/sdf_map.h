#pragma once
#include <map>
#include <unordered_map>
#include <queue>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>

#include "srrg_txt_io/pinhole_image_message.h"
#include "srrg_core_map/cloud.h"

namespace srrg_sdf_map {


class Cell {
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Cell(const Eigen::Vector3i& idx=Eigen::Vector3i::Zero()):
        _idx(idx){

    }

    inline bool operator < (const Cell& c) const {
        for (int i=0; i<3; i++){
            if (_idx[i]<c._idx[i])
                return true;
            if (_idx[i]>c._idx[i])
                return false;
        }
        return false;
    }

    inline bool operator == (const Cell& c) const {
        for (int i=0; i<3; i++)
            if(_idx[i] != c._idx[i])
                return false;
        return true;
    }

    inline void setCenter(Eigen::Vector3f origin, float resolution) {
        _center = origin + _idx.cast<float>()*resolution + Eigen::Vector3f(resolution/2,resolution/2,resolution/2);
    }


    Eigen::Vector3i _idx;
    Eigen::Vector3f _center;
    srrg_core::Vector3fVector _points;
    Cell* _parent;
    int _closest;
    float _distance;
    int _sign;
};


struct QEntry{
    QEntry(Cell* c=0, float d=std::numeric_limits<float>::max()) {
        _cell = c;
        _distance = d;
    }

    inline bool operator < (const QEntry& e) const {
        return e._distance < _distance ;
    }

    float _distance;
    Cell* _cell;
};

struct CellQueue : public std::priority_queue<QEntry> {
    typedef typename std::priority_queue<QEntry>::size_type size_type;
    CellQueue(size_type capacity = 0) { reserve(capacity); }
    inline void reserve(size_type capacity) { this->c.reserve(capacity); }
    inline size_type capacity() const { return this->c.capacity(); }
    inline Cell* top() { return std::priority_queue<QEntry>::top()._cell;}
    inline void push(Cell* c) { return std::priority_queue<QEntry>::push(QEntry(c, c->_distance));}
};


template<typename T>
struct matrix_hash : std::unary_function<T, size_t> {
    std::size_t operator()(T const& matrix) const {
        size_t seed = 0;
        for (size_t i = 0; i < matrix.size(); ++i) {
            auto elem = *(matrix.data() + i);
            seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};
typedef std::unordered_map<Eigen::Vector3i,Cell*,matrix_hash<Eigen::Vector3i> > Vector3iCellPtrMap;

class SdfMap : public Vector3iCellPtrMap {
public:
    SdfMap (float resolution_ = 0.05,
            Eigen::Vector3f origin_ = Eigen::Vector3f::Zero(),
            Eigen::Vector3i dimensions_ = Eigen::Vector3i::Zero());

    inline float resolution(){ return _resolution;}
    inline const Eigen::Vector3i dimensions(){ return _dimensions;}
    inline const Eigen::Vector3f origin(){ return _origin;}
    inline int numCells(){ return _num_cells;}

    inline const Eigen::Vector3i toGrid(const Eigen::Vector3f& point) const {
        return ((point - _origin)*_inverse_resolution).cast<int>();
    }
    inline const Eigen::Vector3f toWorld(const Eigen::Vector3i& cell) const{
        return (_origin + cell.cast<float>() *_resolution);
    }

    inline const float euclideanDistance(const Eigen::Vector3f& a, const Eigen::Vector3f& b){
        return (pow(a.x()-b.x(),2)+pow(a.y()-b.y(),2)+pow(a.z()-b.z(),2));
    }

    inline const int computeSign(const Eigen::Vector3f& a, const Eigen::Vector3f& b){
        return (a.norm()-b.norm() < 0) ? -1 : 1;
    }

    void integrateScan(srrg_core::PinholeImageMessage* image_msg, float max_distance, float min_distance);

protected:
    float _resolution;
    float _inverse_resolution;
    Eigen::Vector3f _origin;
    Eigen::Vector3i _dimensions;
    int _num_cells;

    inline const bool hasCell(const Eigen::Vector3i& idx){
        Vector3iCellPtrMap::iterator it = find(idx);
        return (it != end()) ? true : false;
    }

    int findNeighbors(Cell **neighbors, Cell *c);
};

}
