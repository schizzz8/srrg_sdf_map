#include "sdf_map.h"

namespace std {
template<>
bool std::less<Eigen::Vector3f>::operator ()(const Eigen::Vector3f& a,const Eigen::Vector3f& b) const {
    for(size_t i=0;i<3;++i) {
        if(a[i]<b[i]) return true;
        if(a[i]>b[i]) return false;
    }
    return false;
}
}

using namespace srrg_sdf_map;
using namespace srrg_core;
using namespace srrg_core_map;
using namespace std;

SdfMap::SdfMap(float resolution_, Eigen::Vector3f origin_, Eigen::Vector3i dimensions_)
    :_resolution(resolution_),
      _origin(origin_),
      _dimensions(dimensions_){
    _inverse_resolution = 1./_resolution;
    _num_cells = _dimensions.x()*_dimensions.y()*_dimensions.z();
}

void SdfMap::integrateScan(PinholeImageMessage *image_msg, float max_distance, float min_distance){

    Eigen::Matrix3f invK = image_msg->cameraMatrix().inverse();
    const cv::Mat& image = image_msg->image();
    const float& depth_scale = image_msg->depthScale();
    Eigen::Isometry3f global_transform = image_msg->odometry()*image_msg->offset();
    Eigen::Isometry3f inverse_transform = image_msg->offset().inverse()*image_msg->odometry().inverse();

    for(int r=0; r<image.rows;r++) {
        const unsigned short* id_ptr  = image.ptr<unsigned short>(r);
        for(int c=0; c<image.cols;c++){
            unsigned short id = *id_ptr;
            float d = id * depth_scale;
            Eigen::Vector3f camera_point = invK * Eigen::Vector3f(c*d,r*d,d);

            if (id>0 && d<max_distance && d>min_distance){
                Eigen::Vector3f world_point = global_transform*camera_point;
                Eigen::Vector3i idx = toGrid(world_point);

                if(hasCell(idx)){
                    at(idx)->_points.push_back(world_point);
                    float dist = euclideanDistance(at(idx)->_center,world_point);
                    if(dist < at(idx)->_distance){
                        at(idx)->_distance = dist;
                        at(idx)->_sign = computeSign(inverse_transform*at(idx)->_center,world_point);
                        at(idx)->_closest = at(idx)->_points.size()-1;
                    }
                } else {
                    Cell* cell = new Cell(idx);
                    cell->setCenter(_origin,_resolution);
                    cell->_points.push_back(world_point);
                    float dist = euclideanDistance(cell->_center,world_point);
                    if(dist < cell->_distance) {
                        cell->_distance = dist;
                        cell->_sign = computeSign(inverse_transform*at(idx)->_center,camera_point);
                        cell->_closest = 0;
                    }
                    Vector3iCellPtrMap::iterator it = begin();
                    insert(it,std::pair<Eigen::Vector3i,Cell*>(idx,cell));
                }
            }

            id_ptr++;
        }
    }

    CellQueue q;
    for(Vector3iCellPtrMap::iterator it = begin(); it != end(); ++it) {
        Cell* cell = it->second;
        cell->_parent = cell;
        q.push(cell);
    }

    bool stop = false;
    int loop = 1;
    Cell* last = 0;
    CellQueue dummy(q);
    while (!dummy.empty()) {
        last = dummy.top();
        dummy.pop();
    }

    Cell* neighbors[26];

    while(stop == false && loop > 0) {
        Cell* current = q.top();
        Cell* parent = current->_parent;
        q.pop();
        int k = findNeighbors(neighbors,current);
        for(int ii=0; ii<k; ii++) {
            Cell* child = neighbors[ii];
            Eigen::Vector3f world_point = parent->_points.at(parent->_closest);
            float dist = euclideanDistance(child->_center,world_point);
            if(dist<child->_distance) {
                child->_parent = parent;
                child->_distance = dist;
                child->_closest = parent->_closest;
                child->_sign = computeSign(inverse_transform*parent->_center,inverse_transform*world_point);
                q.push(child);
            }
        }
        if(*last == *current) {
            dummy = q;
            while (!dummy.empty()) {
                last = dummy.top();
                dummy.pop();
            }
            loop--;
        }
    }

}

int SdfMap::findNeighbors(Cell **neighbors, Cell *c) {
    int x = c->_idx.x();
    int y = c->_idx.y();
    int z = c->_idx.z();
    int xmin = (x-1<0) ? 0 : x-1;
    int xmax = (x+1>_dimensions.x()-1) ? _dimensions.x()-1 : x+1;
    int ymin = (y-1<0) ? 0 : y-1;
    int ymax = (y+1>_dimensions.y()-1) ? _dimensions.y()-1 : y+1;
    int zmin = (z-1<0) ? 0 : z-1;
    int zmax = (z+1>_dimensions.z()-1) ? _dimensions.z()-1 : z+1;
    int k=0;
    for(size_t xx=xmin; xx <= xmax; xx++)
        for(size_t yy=ymin; yy <= ymax; yy++)
            for(size_t zz=zmin; zz <= zmax; zz++)
                if(xx != x || yy != y || zz != z) {
                    Eigen::Vector3i idx(xx,yy,zz);
                    if(hasCell(idx)) {
                        neighbors[k] = at(idx);
                        k++;
                    }
                    else {
                        Cell* cell = new Cell(idx);
                        cell->setCenter(_origin,_resolution);
                        Vector3iCellPtrMap::iterator it = begin();
                        insert(it,std::pair<Eigen::Vector3i,Cell*>(idx,cell));
                        neighbors[k] = at(idx);
                        k++;
                    }
                }
    return k;
}
