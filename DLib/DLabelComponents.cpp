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
#include <DLabelComponents.h>

int equiv_find(int *equiv,int i)
{
  int j;
  for(j=i;equiv[j];j=equiv[j]) ;
  int result=j;
  int oldequivj;
  for(j=i;equiv[j];j=oldequivj)
    {
      oldequivj=equiv[j];
      equiv[j]=result;
    }
  return(j);
}

void equiv_union(int *equiv,int a, int b)
{
  if(equiv_find(equiv,a)!=equiv_find(equiv,b))
    {
      //      equiv[a]=b;
      equiv[equiv_find(equiv,a)]=equiv_find(equiv,b);

      printf( "marking %d and %d \n",a,b);
    }
}


//
// Label 8-connected plateaus
//
void DLabelComponents::label_connected_components(const DPlane &in)
{
    int pixels = in.rows() * in.cols();
    int equiv_table[pixels];
    memset(equiv_table,0,pixels*sizeof(int));

    DPlane label(in.rows(),in.cols());

    int nextlabel=1;
    
    for(int i=0; i<in.rows(); i++)
        for(int j=0; j<in.cols(); j++)
            if(in[i][j] == 0)
                label[i][j]=0;
            else 
            {
                bool lB=true, ulB=true, uB=true, urB=true;

		if(i==39 && j == 157)
		  printf("here \n");

                if(j == 0) 
                    lB=ulB=false;
                else if(j == in.cols()-1) 
                    urB=false;
                if(i == 0) 
                    uB=ulB=urB=false;

                  // check left, upper-left, upper, and upper-right neighbors 
                  // mark current pixel with the first label found
                  // and update equivalence table 
                double this_altitude=in[i][j];
                
		//		std::cout << i << " " << j << std::endl;

                if(lB && in[i][j-1] == this_altitude) {
                    label[i][j]=label[i][j-1];
                    if(ulB && in[i-1][j-1] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j],(int)label[i-1][j-1]);
                    if(uB && in[i-1][j] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j], (int)label[i-1][j]);
                    if(urB && in[i-1][j+1] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j], (int)label[i-1][j+1]);
                }
                else if(ulB && in[i-1][j-1] == this_altitude) {
                    label[i][j]=label[i-1][j-1];
                    if(uB && in[i-1][j] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j], (int)label[i-1][j]);
                    if(urB && in[i-1][j+1] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j], (int)label[i-1][j+1]);
                }
                else if(uB && in[i-1][j] == this_altitude) {
                    label[i][j]=label[i-1][j];
                    if(urB && in[i-1][j+1] == this_altitude)
                        equiv_union(equiv_table, (int)label[i][j],(int)label[i-1][j+1]);
                }
                else if(urB && in[i-1][j+1] == this_altitude)
                    label[i][j]=label[i-1][j+1];
                else
                    label[i][j]=nextlabel++;
            }
                
    int *labelmap = new int[nextlabel+1];
    int *lowestmap = new int[nextlabel];
    memset(lowestmap, 0, sizeof(int)*nextlabel);
    int nextlow=1;

    regions.push_back(DRegion());
    regions[0].pix_value = 0;
    for(int i=1; i<nextlabel; i++) {
        int root=equiv_find(equiv_table,i);
        if(lowestmap[root])
            labelmap[i]=lowestmap[root];
        else
        {
            lowestmap[root]=nextlow;
            labelmap[i]=nextlow;
            nextlow++;
        }
        regions.push_back(DRegion());

    }
            
    for(int i=0; i<in.rows(); i++)
        for(int j=0; j<in.cols();j++)
            if(label[i][j] > 0)
            {
                label[i][j]=labelmap[(int)label[i][j]];
                regions[(int)label[i][j]].pix_value = in[i][j];
            }

    int lastused=nextlow-1;
    // +1 for background
    region_count = lastused+1;

    for(int i=0; i<in.rows(); i++)
        for(int j=0; j<in.cols();j++)
        {
            int this_label = (int)label[i][j];

            regions[this_label].count++; 
            
            if(regions[this_label].bounding_box.left() > j)
              regions[this_label].bounding_box.left(j);
            if(regions[this_label].bounding_box.top() > i)
              regions[this_label].bounding_box.top(i);
            if(regions[this_label].bounding_box.bottom() < i)
              regions[this_label].bounding_box.bottom(i);
            if(regions[this_label].bounding_box.right() < j)
              regions[this_label].bounding_box.right(j);
        }

    delete[] labelmap;
    delete[] lowestmap;

    label_plane = label;

}
