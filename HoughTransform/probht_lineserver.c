/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing, Lappeenranta, Finland
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki Kalviainen (Heikki.Kalviainen@lut.fi)
 *
 * For further information, please notice following papers:
 *
 *	Kiryati, N., Eldar, Y., and Bruckstein, A.M, A Probabilistic Hough
 *	Transform, EE PUB No. 746, Department of Eletrical Engineering,
 *	Technion Israel Institute of Technology, Haifa, Israel, February 1990.
 *
 *	Kiryati, N., Eldar, Y., and Bruckstein, A.M, A Probabilistic Hough
 *	Transform, Pattern Recognition, vol. 24, no. 4, 1991, pp. 303-316.
 *
 * File:    	probht_lineserver.c
 * Purpose: 	resolving options for the Probabilistic Hough Transform by Kiryati
 *          	et al. (ProbHT)
 * Date:        Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_formats.h"
#include <strings.h>

pic_type pic[MAX_SIZE][MAX_SIZE], pic1[MAX_SIZE][MAX_SIZE],
         rpic[MAX_SIZE][MAX_SIZE], gpic[MAX_SIZE][MAX_SIZE];
double real_params[MAX_LINES+2][2];
int StopFlag;
char proc_name[256];

void err_prnt(), warn_prnt(), show_usage(), print_options(),
     show_common_usage(), print_common_options();

/*****************************************************************************/
main(argc,argv)
int argc;
char *argv[];
{
  extern char *optarg;
  extern int optind;
  char *cptr, File[256], GrayFile[256], ResultFile[256], RealParamsFile[256],
       ParametersOut[256];
  pic_type Threshold=5;
  int c, dim_x, dim_y, dim_rho=ACC_MAX_SIZE, dim_theta=ACC_MAX_SIZE,
      gray_dim_x, gray_dim_y, TextInfo=0,
      Shrink=0, NumOfMaxs=999, MinSegLen=10, MaxSegWidth=2,
      MaxSegGap=5, FileFormat=PGM_PIC, StoreMaximas=0;
  double Noise=0.0, SampleLevel=20.0;

  File[0]=GrayFile[0]=ResultFile[0]=RealParamsFile[0]='\0';
  cptr=rindex(argv[0],'/');
  strcpy(proc_name,cptr?++cptr:argv[0]);

  while ((c = getopt(argc, argv, "f:G:R:r:F:sn:l:w:g:m:Mihd:t:S:P:u")) != -1)
    switch (c) {
      case 'f': strcpy(File,optarg); 
                break;
      case 'G': strcpy(GrayFile,optarg);
                break;
      case 'R': strcpy(RealParamsFile,optarg);
                break;
      case 'r': strcpy(ResultFile,optarg);
                break;
      case 'F': FileFormat=atoi(optarg);
                break;
      case 's': Shrink=1;
                break;
      case 'n': Noise=atof(optarg); 
                break;
      case 'l': MinSegLen=atoi(optarg); 
                break;
      case 'w': MaxSegWidth=atoi(optarg); 
                break;
      case 'g': MaxSegGap=atoi(optarg);
                break;
      case 'm': NumOfMaxs=atoi(optarg); 
                break;
      case 'M': StoreMaximas=1;
                break;
      case 'i': TextInfo=1;
                break;
      case 'h': print_options();
                exit(0);
      case 'd': dim_rho=atoi(optarg);
                dim_theta=atoi(argv[optind++]);
                break;
      case 't': Threshold=atoi(optarg); 
                break;
      case 'S': SampleLevel=atof(optarg); 
                break;
      case 'P': strcpy(ParametersOut, optarg); 
                break;
      case 'u':
      default : show_usage();
                exit(0);
    }

  if (make_object_from_image(File,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    if (strlen(GrayFile)) {
      pic_from_disk(GrayFile,gpic,Shrink,&gray_dim_x,&gray_dim_y);
      if (dim_x!=gray_dim_x || dim_y!=gray_dim_y)
        warn_prnt("Edge and graylevel pictures not same sized.\n");
      copy_pic(rpic,gpic,dim_x,dim_y);
    }
    if (SampleLevel>0.0 && SampleLevel<100.0)
      take_a_sample_of_pic(pic,pic1,dim_x,dim_y,SampleLevel);
    else
      copy_pic(pic1,pic,dim_x,dim_y);

    if (strlen(RealParamsFile))
      if (real_params_from_disk(RealParamsFile,real_params)!=real_params[0][0])
        warn_prnt("Something curious in real_params array.\n");

    if (dim_rho<3 || dim_rho>ACC_MAX_SIZE) {
      dim_rho=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator size, using default size.\n");
    }
    if (dim_theta<3 || dim_theta>ACC_MAX_SIZE) {
      dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator size, using default size.\n");
    }

    sht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_theta,dim_rho,
             real_params,NumOfMaxs,Threshold,MinSegLen,MaxSegWidth,MaxSegGap,
             TextInfo,ParametersOut);

    if (!(!strlen(ResultFile) && TextInfo))
      if (StoreMaximas)
        pic_to_disk(ResultFile,gpic,dim_x,dim_y,FileFormat,proc_name);
      else
        pic_to_disk(ResultFile,rpic,dim_x,dim_y,FileFormat,proc_name);
  }

}

/*****************************************************************************/
void show_usage()
{
  fprintf(stderr,"Usage: %s ",proc_name);
  show_common_usage(stderr);
  fprintf(stderr,"[-d dim_rho dim_theta] [-t Threshold] [-S SampleLevel]\n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"Probabilistic Hough Transform (ProbHT) for lines.\n\n");
  print_common_options(stdout);
  fprintf(stdout,"\t-d dim_rho dim_theta\t\t\t\tdef: %d %d\n",
          ACC_MAX_SIZE,ACC_MAX_SIZE);
  fprintf(stdout,"\t   Accumulator array size (max %d * %d).\n",
          ACC_MAX_SIZE,ACC_MAX_SIZE);
  fprintf(stdout,"\t-t Threshold\t\t\t\t\tdef: 5\n");
  fprintf(stdout,"\t   Minimum score accepted for a accumulator maximum.\n");
  fprintf(stdout,"\t-S SampleLevel\t\t\t\t\tdef: 20.0\n");
  fprintf(stdout,"\t   The number of edge points selected at random, \n");
  fprintf(stdout,"\t   as the percentage of all edge points (%%).\n");
}
