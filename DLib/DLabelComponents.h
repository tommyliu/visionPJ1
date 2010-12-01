//
//  DLib: A simple image processing library.
//
//  David Crandall, 2003-2005
//  crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
//
//
#include <DPlane.h>
#include <vector>

class DRegion
{
  public:
    DRegion(DRect &bb, double val) : bounding_box(bb), pix_value(val) {}
    DRegion() 
    { 
      count = 0; 
      bounding_box = DRect(1000000,1000000,-1,-1);
    }

    DRect getBoundingBox() { return bounding_box; } 
    double getPixelValue() { return pix_value; }
    bool isValid() { return true; }

    DRect bounding_box;
    double pix_value;
    int count;
};


class DLabelComponents
{
  public:
    DLabelComponents(const DPlane &in_image)
    {
      label_connected_components(in_image);
    }

    void label_connected_components(const DPlane &in);
    int size() { return region_count; }
    
    DRegion &operator[](int index) { return(regions[index]); }

    DPlane GetLabelPlane() { return label_plane; }

  protected:
    int region_count;
    std::vector<DRegion> regions;    
    DPlane label_plane; 
};
