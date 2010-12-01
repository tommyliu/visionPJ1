//
//  k-Fan matching and training code, as described in the paper:
//
//    Crandall, Felzenszwalb, Huttenlocher, "Spatial 
//      Priors for Part-based Recognition using Statistical Models," CVPR 2005.
//     
//
//  Author: David Crandall, 2003-2005
//          crandall@cs.cornell.edu
//
//  Please do not redistribute this code.
//
//
// 
// 
#define TYPE float
#include <Appear.h>
#include <KFan.h>
#include <Train_KFan.h>
#include <Train_Appear.h>
#include <unistd.h>
#include <stdlib.h>

#ifdef DEBUG
std::ostream dout(std::cout.rdbuf()); 
#else 
std::ostream dout(static_cast<std::streambuf*>(0)); 
#endif 

using namespace std;

bool analyze_g = false;
char *appear_file_g = 0, *spatial_file_g = 0, *outfile_g = 0;
int img_count_g = 0;
char **img_list_g = 0;
TYPE cm_ratio_g=0.01;
char *train_file_g = 0;
int K_g=1;
int tile_rows_g=50, tile_cols_g = 50;
TYPE quantile_g = 1.0;

pair<TYPE, _DMatrix<TYPE> > match_process_image(char *fname, AppearModel<TYPE> *appear, KFan<TYPE> &kfan);

int this_proc=0;


unsigned char *load_image(char *fname, int &rows, int &cols)
{
  FILE *fp = fopen(fname, "rb");
  if(!fp)
    throw std::string(string("can't open file ") + string(fname) + string("\n"));

  fread(&rows, 1, sizeof(int), fp);
  fread(&cols, 1, sizeof(int), fp);
  
  unsigned char *in_im = new unsigned char[rows*cols];
  fread(in_im, rows*cols, sizeof(char), fp);
  fclose(fp);

  return in_im;
}

pair<TYPE, _DMatrix<TYPE> > match_process_image(char *fname, AppearModel<TYPE> *appear, KFan<TYPE> &kfan)
{
  pair<TYPE, _DMatrix<TYPE> > result;
  
  int rows, cols;
  
  dout <<  "  loading image " << endl;
  unsigned char *in_im = load_image(fname, rows, cols);
  
  
  int f_rows = appear->rows();
  int f_cols = appear->cols();
  dout <<  "  calling appear " << endl;
  unsigned char *in_im2 = appear->embed_img(in_im, rows, cols, rows+f_rows, 
					    cols+f_cols);
  
  _DMultiDMatrix<TYPE> _cost_maps;
  _cost_maps = appear->run_model(in_im2, rows+f_rows, cols+f_cols);
  
  _DMultiDMatrix<TYPE> cost_maps(3, _cost_maps.planes(),0,0);
  for(int i=0; i<_cost_maps.planes(); i++)
    {
      dout << i << endl;
      DRect rect(f_rows/2, f_cols/2, rows + f_rows - f_rows/2-1, cols +f_cols - f_cols/2-1);
      cost_maps.get(i) = _cost_maps.get(i).extract(rect) * cm_ratio_g;
    }

  //  cout << cost_maps[0] << endl;
  dout <<  "  calling kfan " << endl;
  result = kfan.do_matching(cost_maps);
  dout <<  "  calling kfan done " << endl;


  delete[] in_im2;
  delete[] in_im;

  return result;
}


class TrainInput
{
public:
  int image_no;
  string image_file;
  TYPE likelihood;
  DMatrix localization;
};

int compar_ti(const void *a, const void *b)
{
  TrainInput *A = (TrainInput *)a;
  TrainInput *B = (TrainInput *)b;

  if(A->likelihood > B->likelihood) return -1;
  else if(A->likelihood < B->likelihood) return 1;
  else return 0;
}

void do_training(int tile_rows, int tile_cols)
{
  dout << "in do_training" << endl;
  ifstream ifs;
  ifs.open(train_file_g, ios::in);
  assert(ifs.is_open());

  _DMatrix<TYPE> part_locations;
  vector<string> filenames;

  int part_count;

  int i_count=0;
  TrainInput inputs[2000];
  while(!ifs.eof())
    {  
      int image_no;
      string image_file;
      TYPE likelihood;
      DMatrix localization;
      
      ifs >> image_no >> image_file >> likelihood;
      if(ifs.eof())
	break;

      ifs >> localization;

      if(image_no / 100000 == 7 || likelihood < -1e90)
	continue; // bg image

      inputs[i_count].image_no = image_no;
      inputs[i_count].image_file = image_file;
      inputs[i_count].likelihood = likelihood;
      inputs[i_count].localization = localization;

      i_count++;
    }

  qsort(inputs, i_count, sizeof(TrainInput), compar_ti);

  part_count = inputs[0].localization.rows();
  part_locations = _DMatrix<TYPE>(int(i_count*quantile_g), part_count*2);

  for(int j=0; j<int(i_count*quantile_g); j++)
    {
      
      for(int i=0; i<part_count; i++)
	{
	  part_locations[j][i*2] = inputs[j].localization[i][1];
	  part_locations[j][i*2+1] = inputs[j].localization[i][0];
	}

      filenames.push_back(inputs[j].image_file);

    }

  //  part_locations = part_locations.extract(DRect(0, 0, i_count-1,  
  //						part_count*2-1));



  // first train appearance model
  Train_Appear<TYPE> ta;
  const int dilate=5;

  cout << "train appear" << endl;
  Appear<TYPE> appear = ta.train(part_locations, filenames, dilate, 
			   tile_rows, tile_cols);

  Train_KFan<TYPE> tk;
  cout << "train kfan" << endl;
  KFan<TYPE> kfan = tk.train(part_locations, K_g);

  char temp[1024];
  sprintf(temp, "%s.appear", outfile_g);
  ofstream ofs;
  ofs.open(temp, ios::out);
  ofs << appear << endl;
  ofs.close();

  sprintf(temp, "%s.%dfan", outfile_g, K_g);
  ofs.open(temp, ios::out);
  ofs << kfan << endl;
  ofs.close();

}

int double_compar(const void *a, const void *b)
{
  TYPE A = *((TYPE *)(a));
  TYPE B = *((TYPE *)(b));

  if(A<B) return -1;
  else if(A>B) return 1;
  else return 0;
}

TYPE analyze(char *description)
{
  ifstream ifs;
  ifs.open(outfile_g, ios::in);
  
  int img;
  TYPE likelihood;
  DMatrix result;
  TYPE tp_like[img_count_g], fp_like[img_count_g];
  int tp_count=0, fp_count=0;
  
  while(!ifs.eof())
    {

      char temp[1024];
      ifs >> img;
      
      ifs.getline(temp,1024);
      ifs.getline(temp,1024);
      
      ifs >> likelihood >> result;
      
      if(img / 100000 == 7)
	fp_like[fp_count++] = likelihood;
      else
	tp_like[tp_count++] = likelihood;
    }
  ifs.close();
  if(fp_count != tp_count)
    cout << "ANALYSIS NOT RELIABLE SINCE fp_count != tp_count" << endl;
  
  qsort(fp_like, fp_count, sizeof(TYPE), double_compar);
  qsort(tp_like, tp_count, sizeof(TYPE), double_compar);
  
  int i;
  for(i=0; i<tp_count; i++)
    {
      if(tp_like[i] > fp_like[fp_count-i-1])
	break;
    }
  
  cout << "equal ROC point: " << endl;
  cout << "  threshold = " << tp_like[i] << endl;
  cout << "  correct detect = " << 1.0 - TYPE(i) / tp_count << endl;
  cout << "  false alarm = " << TYPE(i) / fp_count << endl;

  char temp[1024];
  sprintf(temp, "%s.analysis", outfile_g);
  FILE *fp = fopen(temp, "a");
  if(fp)
    {
      fprintf(fp, "--- Results of run on %s\n",description);
      fprintf(fp, "models: %s %s\n", appear_file_g, spatial_file_g);
      fprintf(fp, "equal ROC point:\n");
      fprintf(fp, "  threshold = %f\n",tp_like[i]);
      fprintf(fp, "  correct detect = %f\n",1.0 - TYPE(i) / tp_count);
      fprintf(fp, "  false alarm = %f\n", TYPE(i) / fp_count);
      fprintf(fp, "\n\n");
      fclose(fp);
    }

  return (1.0 - TYPE(i)/ tp_count);
}

void help()
{
  cout << "* " << endl;
  cout << "* This is the simple (\"lite\") version of the k-fan matching and training code. :)" << endl;
  cout << "* See Crandall, Felzenszwalb, Huttenlocher, 'Spatial Priors for Part-based Recognition " << endl;
  cout << "*      using Statistical Models', CVPR 2005 for more details." << endl;
  cout << "* " << endl;
  cout << "* It only supports single-threaded matching and supervised training. Use the full version" << endl;
  cout << "* (match) for multi-threaded runs, unsupervised learning, generating ROC curves, affine invariance, etc." << endl;
  cout << "*" << endl;
  cout << "* usage: match_lite [options] input_file1 input_file2 ... " << endl;
  cout << "*" << endl;
  cout << endl;
  cout << "* K-FAN MATCHING: (localization) " << endl;
  cout << "   usage: match_lite [options] -a appear_model -s kfan_model  edge_map1 edge_map2 ... " << endl;
  cout <<  endl;
  cout << "   -c n        set cm ratio (default 0.01)" << endl;
  cout << "   -A          output some simple analysis of run" << endl;
  cout <<  endl;
  cout <<  endl;
  cout << "* SUPERVISED LEARNING: " << endl;
  cout << "   usage: match [options] -t train_partlist_pt edge_map1 edge_map2 ... " << endl;
  cout <<  endl;
  cout << "   -K n        set number of distinguished parts (k) (default 1)" << endl;
  cout << "   -q f        only use f% highest-likelihood training data (default 1.0)" << endl;
  cout << "   -z n        use n x n patches for appearance models (default 50)" << endl;
  cout <<  endl;
  cout <<  endl;
}

int parse_opts(int argc, char *argv[])
{
  int ii=1;
  for(ii=1; ii<argc; ii++)
    {
      if(!strcmp(argv[ii], "-a"))
	appear_file_g = argv[++ii];
      else if(!strcmp(argv[ii], "-A"))
	analyze_g = true;
      else if(!strcmp(argv[ii], "-s"))
	spatial_file_g = argv[++ii];
      else if(!strcmp(argv[ii], "-K"))
	K_g = atoi(argv[++ii]);
      else if(!strcmp(argv[ii], "-o"))
	outfile_g = argv[++ii];
      else if(!strcmp(argv[ii], "-t"))
	train_file_g = argv[++ii];
      else if(!strcmp(argv[ii], "-c"))
	cm_ratio_g=atof(argv[++ii]);
      else if(!strcmp(argv[ii], "-q"))
	quantile_g=atof(argv[++ii]);
      else if(!strcmp(argv[ii], "-z"))
	tile_rows_g = tile_cols_g = atoi(argv[++ii]);
      else
	break;
    }

  img_count_g = argc-ii;
  img_list_g = argv+ii;

  return 0;
}


int main(int argc, char *argv[])
{
  try {
  if(argc == 1)
    {
      help();
      return 1;
    }

  parse_opts(argc, argv);  

  if(train_file_g)
    do_training(tile_rows_g, tile_cols_g);

  if(appear_file_g && spatial_file_g)
    {
      // make sure output file doesn't exist already
      if(outfile_g)
	{
	  FILE *fp = fopen(outfile_g,"r");
	  if(fp)
	    {
	      cout << "output file " << outfile_g << " already exists! " 
		   << endl;
	      cout << "delete it first, then rerun" << endl;
	      fclose(fp);
	      exit(1);
	    }
	}

      AppearModel<TYPE> *appear;
      KFan<TYPE> *kfan;

      dout << "loading appear" << endl;
      appear= new Appear<TYPE>(appear_file_g);
      
      dout << "loading spatial" << endl;
      KFan<TYPE> fan(spatial_file_g);
    
      switch(fan.K)
	{
	case 0:
	  kfan=new K_0Fan<TYPE>(spatial_file_g);
	  break;
	case 1:
	  kfan=new K_1Fan<TYPE>(spatial_file_g);
	  break;
	default:
	  assert(0);
	}
      


      for(int i=0; i<img_count_g; i++)
	{
	  cout << img_list_g[i] << endl;
	  pair<TYPE, _DMatrix<TYPE> > result = 
	    match_process_image(img_list_g[i], appear, *kfan);
	  
	  if(outfile_g)
	    {
	      ofstream ofs;
	      ofs.open(outfile_g,ios::out | ios::app);
	      ofs << atoi(basename(img_list_g[i])) << endl;
	      ofs << img_list_g[i] << endl;
	      ofs << result.first << endl;
	      ofs << result.second << endl ;
	      ofs.close();
	    }
	  else
	    {
	      cout << atoi(basename(img_list_g[i])) << endl;
	      cout << img_list_g[i] << endl;
	      cout << result.first << endl;
	      cout << result.second << endl ;
	    }
	}
      
      delete appear;
      delete kfan;

      if(outfile_g)
	analyze(outfile_g);
    }

  }
  catch(string &str)
    {
      cout << str << endl;
      exit(1);
    }
}
