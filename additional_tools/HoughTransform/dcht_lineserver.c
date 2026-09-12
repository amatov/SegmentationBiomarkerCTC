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
 * File:    	dcht_lineserver.c
 * Purpose: 	resolving options for the Dynamic Combinatorial Hough Transform (DCHT)
 * Date:  	Jun 1, 1993
 * Last change:	Oct 10, 1995
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
  int c, dim_x, dim_y, dim_theta=ACC_MAX_SIZE,
      gray_dim_x, gray_dim_y, TextInfo=0, Shrink=0,
      NumOfTests=1, NumOfMaxs=999, MinSegLen=10, MaxSegWidth=2, MaxSegGap=5,
      Threshold=2, FileFormat=PGM_PIC, StoreMaximas=0;
  double Noise=0.0;

  File[0]=GrayFile[0]=ResultFile[0]=RealParamsFile[0]='\0';
  cptr=rindex(argv[0],'/');
  strcpy(proc_name,cptr?++cptr:argv[0]);

  while ((c = getopt(argc, argv, "f:G:R:r:F:sn:l:w:g:m:Mihd:t:c:uP:")) != -1)
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
      case 'd': dim_theta=atoi(optarg);
                break;
      case 't': Threshold=atoi(optarg);
                break;
      case 'c': NumOfTests=atoi(optarg);
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

    if (strlen(RealParamsFile))
      if (real_params_from_disk(RealParamsFile,real_params)!=real_params[0][0])
        warn_prnt("Something curious in real_params array.\n");

    if (dim_theta<3 || dim_theta>ACC_MAX_SIZE) {
      dim_theta=ACC_MAX_SIZE;
      warn_prnt("Illegal accumulator size, using default size.\n");
    }

    dcht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_theta,real_params,NumOfMaxs,
              MinSegLen,MaxSegWidth,MaxSegGap,Threshold,NumOfTests,TextInfo,
	      ParametersOut);

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
  fprintf(stderr,"[-d dim_theta] [-t Threshold] [-c NumOfTests]\n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"Dynamic Combinatorial Hough Transform (DCHT) for lines.\n\n");
  print_common_options(stdout);
  fprintf(stdout,"\t-d dim_theta\t\t\t\t\tdef: %d\n",ACC_MAX_SIZE);
  fprintf(stdout,"\t   Histogram vector length (max %d).\n",ACC_MAX_SIZE);
  fprintf(stdout,"\t-t Threshold\t\t\t\t\tdef: 2\n");
  fprintf(stdout,"\t   Minimum score accepted for a accumulator maximum.\n");
  fprintf(stdout,"\t-c NumOfTests\t\t\t\t\tdef: 1\n");
  fprintf(stdout,"\t   The number of tests to run.\n");
}
