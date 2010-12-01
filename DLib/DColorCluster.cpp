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
#include <DColorCluster.h>

// k clusters per channel = k^3 clusters
DPlane DColorCluster::do_clustering(DImage &input, int k, bool ignore_black) 
{
    int cluster_count = k*k*k;
    for(int i=0; i<cluster_count; i++)
      clusters.push_back(DRGBCluster());

    DPlane result(input.rows(), input.cols());
    
    int delta = 256/(k*2);

    for(int a=0, cl=0; a<k; a++)
        for(int b=0; b<k; b++)
            for(int c=0; c<k; c++, cl++)
                clusters[cl].set_mean(
                    DTriple((a+1)*delta, (b+1)*delta, (c+1)*delta));
               

    int changes=1000000000;

    bool done=false;

    while(!done) {
        changes=0;

    for(int i=0; i<input.rows(); i++)
        for(int j=0; j<input.cols(); j++)
        {
            if(ignore_black && input[0][i][j] == 0 && input[1][i][j] == 0 &&
               input[2][i][j] == 0)
            {
                result[i][j] = 0;
                continue;
            }
            int closest_cluster=0;
            double min_dist=1000000000;
            DTriple sample(input[0][i][j], input[1][i][j], input[2][i][j]);
            for(int c=0; c<cluster_count; c++)
                if(clusters[c].distance_to(sample) < min_dist)
                {
                    min_dist = clusters[c].distance_to(sample);
                    closest_cluster = c;
                }
            clusters[closest_cluster].add_sample(sample);
            if(closest_cluster+1 != result[i][j]) 
            {
                result[i][j] = closest_cluster+1;
                changes++;
            }
        }

    if( changes < input.rows() * input.cols() * 0.001) 
        done=true;
    else
        for(int i=0; i<cluster_count; i++) 
        {
            clusters[i].update_mean();
            clusters[i].remove_all();
        }
    
    printf("there were %d changes\n",changes);
    }

    return result;
}

DPlane DColorCluster::apply_clustering(const DImage &input)
{
    DPlane result(input.rows(), input.cols());
    
    for(int i=0; i<input.rows(); i++)
        for(int j=0; j<input.cols(); j++)
        {
            int closest_cluster=0;
            double min_dist=1000000000;
            DTriple sample(input[0][i][j], input[1][i][j], input[2][i][j]);
            for(int c=0; c<clusters.size(); c++)
                if(clusters[c].distance_to(sample) < min_dist)
                {
                    min_dist = clusters[c].distance_to(sample);
                    closest_cluster = c;
                }
            result[i][j] = closest_cluster+1;
        }    
    
    return result;
}
