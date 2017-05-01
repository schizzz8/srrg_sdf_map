#include <cmath>
#include <fstream>
#include "srrg_system_utils/system_utils.h"
#include "srrg_txt_io/message_reader.h"
#include "srrg_txt_io/pinhole_image_message.h"
#include "srrg_sdf_map/sdf_map.h"

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <vtkSmartPointer.h>
#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkXMLImageDataWriter.h>


using namespace std;
using namespace srrg_core;
using namespace srrg_sdf_map;



// Help objects to force linking
PinholeImageMessage p;

const char* banner[] = {
    "srrg_sdf_map_example: example on how to integrate depth images into a signed distance function",
    "",
    "usage: srrg_sdf_map_example [options] <dump_file>",
    0
};

float resolution = 0.05f, range=5.0f, max_distance = 7.0f, min_distance = 0.4f;


void readParameters(string filename){
    cerr << "Reading configuration file!" << endl << endl;
    cv::FileStorage fs;
    fs.open(filename, cv::FileStorage::READ);

    resolution = (float) fs["Resolution"];
    cerr << "resolution: " << resolution << endl;

    max_distance = (float) fs["MaxDistance"];
    cerr << "max distance: " << max_distance << endl;

    min_distance = (float) fs["MinDistance"];
    cerr << "min distance: " << min_distance << endl;

    range = (float) fs["Range"];
    cerr << "range: " << range << endl;

    cerr << endl;
}

void writeToVtkFile(const SdfMap& map) {
    vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    imageData->SetDimensions(map.dimensions().x(),map.dimensions().y(),map.dimensions().z());
    imageData->SetOrigin(map.origin().x(),map.origin().y(),map.origin().z());
    imageData->SetSpacing(map.resolution(),map.resolution(),map.resolution());
    imageData->SetExtent(0,map.dimensions().x()-1,0,map.dimensions().y()-1,0,map.dimensions().z()-1);
#if VTK_MAJOR_VERSION <= 5
    imageData->SetNumberOfScalarComponents(1);
    imageData->SetScalarTypeToFloat();
#else
    imageData->AllocateScalars(VTK_FLOAT, 1);
#endif

    vtkSmartPointer<vtkFloatArray> signed_distance = vtkSmartPointer<vtkFloatArray>::New();
    signed_distance->SetName("signed_distance");
    signed_distance->SetNumberOfComponents(1);

    for(int z = 0; z < map.dimensions().z(); z++)
        for(int y = 0; y < map.dimensions().y(); y++)
            for(int x = 0; x < map.dimensions().x(); x++){
                Eigen::Vector3i idx = Eigen::Vector3i (x,y,z);
                if (map.hasCell(idx)) {
                    signed_distance->InsertNextValue(map.at(idx)->_sign*map.at(idx)->_distance);
                }
                else {
                    signed_distance->InsertNextValue(numeric_limits<float>::quiet_NaN());
                }
            }
    imageData->GetPointData()->AddArray(signed_distance);
    imageData->Update();

    vtkSmartPointer<vtkXMLImageDataWriter> writer = vtkSmartPointer<vtkXMLImageDataWriter>::New();
    writer->SetFileName("solution.vti");
#if VTK_MAJOR_VERSION <= 5
    writer->SetInputConnection(imageData->GetProducerPort());
#else
    writer->SetInputData(imageData);
#endif
    writer->Write();

}


int main(int argc, char ** argv) {
    if (argc < 2 || !strcmp(argv[1], "-h")) {
        printBanner(banner);
        return 0;
    }

    string settings_filename = argv[1];
    string input_filename = argv[2];

    readParameters(settings_filename.c_str());

    Eigen::Vector3f lower,higher;
    float xmin=std::numeric_limits<float>::max();
    float xmax=std::numeric_limits<float>::min();
    float ymin=std::numeric_limits<float>::max();
    float ymax=std::numeric_limits<float>::min();
    float zmin=std::numeric_limits<float>::max();
    float zmax=std::numeric_limits<float>::min();

    vector<PinholeImageMessage*> image_msgs;

    MessageReader reader;
    reader.open(input_filename.c_str());

    BaseMessage* msg = 0;
    while ((msg = reader.readMessage())) {
        msg->untaint();
        PinholeImageMessage* image_msg = dynamic_cast<PinholeImageMessage*>(msg);
        if (image_msg) {
            float x = (image_msg->odometry()*image_msg->offset()).translation().x();
            float y = (image_msg->odometry()*image_msg->offset()).translation().y();
            float z = (image_msg->odometry()*image_msg->offset()).translation().z();

            xmax = xmax > x+range ? xmax : x+range;
            ymax = ymax > y+range ? ymax : y+range;
            zmax = zmax > z+range ? zmax : z+range;
            xmin = xmin < x-range ? xmin : x-range;
            ymin = ymin < y-range ? ymin : y-range;
            zmin = zmin < z-range ? zmin : z-range;

            image_msgs.push_back(image_msg);
        }
    }
    if(image_msgs.size() == 0)  {
        cerr << "No depth images found ... quitting!" << endl;
        return 0;
    }
    cerr << "Read " << image_msgs.size() << " depth images"<< endl;

    lower = Eigen::Vector3f (xmin,ymin,zmin);
    higher = Eigen::Vector3f (xmax,ymax,zmax);
    Eigen::Vector3f origin = lower - Eigen::Vector3f(5*resolution,5*resolution,5*resolution);
    Eigen::Vector3i dimensions = (((higher + Eigen::Vector3f(5*resolution,5*resolution,5*resolution))-lower)/resolution).cast<int> ();

    cerr << endl << "Bounding Box: " << endl;
    cerr << "\t>>Lower: " << lower.transpose() << endl;
    cerr << "\t>>Higher: " << higher.transpose() << endl;

    cerr << "SdfMap size: " << dimensions.x() << "x" << dimensions.y() << "x" << dimensions.z() << endl;
    if(dimensions.x() == 0 || dimensions.y() == 0 || dimensions.z() == 0) {
        cerr << "Zero map size ... quitting!" << endl;
        return 0;
    }

    cerr << "Building the SdfMap" << endl;
    SdfMap sdf_map (resolution,origin,dimensions);
    for(size_t i = 0; i < image_msgs.size(); ++i) {
        PinholeImageMessage* msg = image_msgs[i];
        sdf_map.integrateScan(msg);
        cerr << ".";
    }
    cerr << endl << "done" << endl;

    cerr << "Saving to disk..." << endl;
    writeToVtkFile(sdf_map);

    return 0;
}
