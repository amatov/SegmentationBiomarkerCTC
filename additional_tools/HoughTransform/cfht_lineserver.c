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
 * File:         cfht_lineserver.c
 * Purpose:      resolving options for the Curve Fitting Hough Transform (CFHT)
 * Date:         Jun 1, 1993
 * Last change:  Oct 10, 1995
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
  int c, dim_x, dim_y, gray_dim_x, gray_dim_y, TextInfo=0, Shrink=0, Threshold=2,
      NumOfTests=1, NumOfMaxs=999, Msize=3, ThreshForPoints=0, DontRemovePoints=0,
      MinSegLen=10, MaxSegWidth=2, MaxSegGap=5, FileFormat=PGM_PIC,
      StoreMaximas=0, PostScript=0;
  double TolForFitting=2.0, Epsilon=1.0, Gamma=0.0, Noise=0.0;

  File[0]=GrayFile[0]=ResultFile[0]=RealParamsFile[0]='\0';
  cptr=rindex(argv[0],'/');
  strcpy(proc_name,cptr?++cptr:argv[0]);

  while ((c = getopt(argc,argv,"f:G:R:r:F:sn:l:w:g:m:Mihc:W:T:E:e:A:Dt:puP:"))!=-1)
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
      case 'c': NumOfTests=atoi(optarg); 
                break;
      case 'W': Msize=atoi(optarg); 
                break;
      case 'T': ThreshForPoints=atoi(optarg); 
                break;
      case 'E': TolForFitting=atof(optarg); 
                break;
      case 'e': Epsilon=atof(optarg); 
                break;
      case 'A': Gamma=atof(optarg); 
                break;
      case 'D': DontRemovePoints=1; 
                break;
      case 't': Threshold=atoi(optarg); 
                break;
      case 'p': PostScript=1;
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

    if (ThreshForPoints==0)
      ThreshForPoints=Msize;

    cfht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
              MaxSegWidth,MaxSegGap,Msize,ThreshForPoints,TolForFitting,Threshold,
              Epsilon,Gamma,DontRemovePoints,NumOfTests,TextInfo,PostScript,
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
  fprintf(stderr,"[-c NumOfTests] [-W Msize] [-T ThreshForPoints] ");
  fprintf(stderr,"[-E TolForFitting] [-e Epsilon] [-A Gamma] [-D] ");
  fprintf(stderr,"[-t Threshold] [-p]\n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"Curve Fitting Hough Transform (CFHT) for lines.\n\n");
  print_common_options(stdout);
  fprintf(stdout,"\t-c NumOfTests\t\t\t\t\tdef: 1\n");
  fprintf(stdout,"\t   The number of tests to run. This is only to get \n");
  fprintf(stdout,"\t   more reliable computation times.\n");
  fprintf(stdout,"\t-W Msize\t\t\t\t\tdef: 3\n");
  fprintf(stdout,"\t   The size of the window ((Msize*2+1)^2).\n");
  fprintf(stdout,"\t-T ThresholdForPoints\t\t\tdef: Msize\n");
  fprintf(stdout,"\t   Threshold for the number of points in the window.\n");
  fprintf(stdout,"\t   If there are less points in the window,\n");
  fprintf(stdout,"\t   the line fitting is not performed.\n");
  fprintf(stdout,"\t-E TolForFitting\t\t\t\tdef: 2.0\n");
  fprintf(stdout,"\t   Maximum mean square error for curve fitting.\n");
  fprintf(stdout,"\t-e Epsilon\t\t\t\t\tdef: 1.0\n");
  fprintf(stdout,"\t   Tolerance to existing accumulator cells.\n");
  fprintf(stdout,"\t   Tolerance of a new candidate cell to near \n");
  fprintf(stdout,"\t   accumulator cells when updating the accumulator.\n");
  fprintf(stdout,"\t   A new candidate cell is merged into an old one, \n"); 
  fprintf(stdout,"\t   if it satisfies the tolerance; otherwise, a new cell \n");
  fprintf(stdout,"\t   is created. The larger the tolerance, the less new cells \n");
  fprintf(stdout,"\t   are created. Tolerance is the Euclidean distance \n");
  fprintf(stdout,"\t   between positions of cells.\n");
  fprintf(stdout,"\t-A Gamma\t\t\t\t\tdef: 0.0\n");
  fprintf(stdout,"\t   Value for weighting accumulator cell position. Old ");
  fprintf(stdout,"position is\n\t   weighted by (1- Gamma). If Gamma is 0.0, ");
  fprintf(stdout,"the old position is\n\t   weighted by it's score, and ");
  fprintf(stdout,"the new one by one.\n");
  fprintf(stdout,"\t-D\t\t\t\t\t\tdef: not in use\n");
  fprintf(stdout,"\t   Do not remove the center pixel of the window from the ");
  fprintf(stdout,"picture.\n\t   Otherwise, if there is not ThresholdForPoints ");
  fprintf(stdout,"edge pixels\n\t   in the window, remove the center pixel.\n");
  fprintf(stdout,"\t-t Threshold\t\t\t\t\tdef: 2\n");
  fprintf(stdout,"\t   Minimun score accepted for a accumulator maximum.\n");
  fprintf(stdout,"\t-p Threshold\t\t\t\t\tdef: not in use\n");
  fprintf(stdout,"\t   PostScript output (houghn.ps; n=1...N) for dyn. accu.\n");
  fprintf(stdout,"\t   PS files illustrate the dynamic accu representation.\n");
}
