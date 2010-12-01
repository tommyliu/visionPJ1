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
#include <DImage.h>
#include <DTriple.h>
#include <vector>

using namespace std;

class DRGBCluster 
{
    public:

        DRGBCluster() {}

        DRGBCluster(DTriple _mean)
        {
            mean = _mean;
        }
        
        DTriple get_mean() 
        {
            return mean;
        }

        void update_mean()
        {
            DTriple sum(0,0,0);
            int count=0;

            vector<DTriple>::iterator iter=members.begin();
            for(; iter != members.end(); ++iter, count++)
                sum += *iter;
            sum /= count;
            mean=sum;
            return;
        }

        void add_sample(const DTriple &sample)
        {
            members.push_back(sample);
        }

        double distance_to(const DTriple &sample)
        {
            DTriple m = get_mean();
            return(sqrt((m.r-sample.r)*(m.r-sample.r) +
                        (m.g-sample.g)*(m.g-sample.g) +
                        (m.b-sample.b)*(m.b-sample.b)));
        }

        void rewind()
        {
            member_iter = members.begin();
        }

        DTriple get_current_sample()
        {
            return *member_iter;
        }

        void next_sample()
        {
            ++member_iter;
        }

        void remove_all()
        {
            members.erase(members.begin(), members.end());
        }


        void set_mean(const DTriple &new_mean)
        {
            mean = new_mean;
        }

        char label[10];

    protected:
        vector<DTriple> members;
        vector<DTriple>::iterator member_iter;
        DTriple mean;
};


class DColorCluster
{
 public:
  DColorCluster(DImage &input, int k, bool ignore_black)
  {
    cluster_plane = do_clustering(input, k, ignore_black);
  }

  DPlane do_clustering(DImage &input, int k, bool ignore_black);
  DPlane apply_clustering(const DImage &input);

  DPlane get_clustered_plane()
  {
    return cluster_plane;
  }

 protected:
  DPlane cluster_plane;
  vector<DRGBCluster> clusters; 
};
