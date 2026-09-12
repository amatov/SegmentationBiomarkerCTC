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
 * File:    	rht_lineserver.c
 * Purpose: 	resolving options for the Randomized Hough Transform (RHT)
 * Date:    	Jun 1, 1993
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
  int c, dim_x, dim_y, gray_dim_x, gray_dim_y, TextInfo=0,
      Shrink=0, NumOfTests=1, NumOfMaxs=999, MinSegLen=10, MaxSegWidth=2,
      MaxSegGap=5, MinDist=2, MaxDist=20, Threshold=2, FileFormat=PGM_PIC,
      StoreMaximas=0, PostScript=0, MakeNoise=0;
  double Accuracy=0.01, Noise=0.0;

  File[0]=GrayFile[0]=ResultFile[0]=RealParamsFile[0]='\0';
  cptr=rindex(argv[0],'/');
  strcpy(proc_name,cptr?++cptr:argv[0]);

  while ((c = getopt(argc, argv, "f:G:R:r:F:sn:l:w:g:m:Mihc:d:t:a:pN:uP:")) != -1)
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
      case 'N': MakeNoise=1;
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
      case 'd': MinDist=atoi(optarg);
                MaxDist=atoi(argv[optind++]);
                break;
      case 't': Threshold=atoi(optarg);
                break;
      case 'a': Accuracy=atof(optarg);
                break;
      case 'p': PostScript=1;
                break;
      case 'P': strcpy(ParametersOut, optarg); 
                break;
      case 'u':
      default : show_usage();
                exit(0);
    }

  if (NumOfTests>1 && PostScript) {
    warn_prnt("No PS files with multiple tests.\n");
    PostScript=0;
  }

  if (make_object_from_image(File,pic1,Shrink,&dim_x,&dim_y,Noise)) {

    if (strlen(GrayFile)) {
      pic_from_disk(GrayFile,gpic,Shrink,&gray_dim_x,&gray_dim_y);
      if (dim_x!=gray_dim_x || dim_y!=gray_dim_y)
        warn_prnt("Edge and graylevel pictures not same sized.\n");
      copy_pic(rpic,gpic,dim_x,dim_y);
    }
    if (MakeNoise) {
      pic_to_disk(ResultFile,pic1,dim_x,dim_y,FileFormat,proc_name);
      exit(0);
    }

    if (strlen(RealParamsFile))
      if (real_params_from_disk(RealParamsFile,real_params)!=real_params[0][0])
        warn_prnt("Something curious in real_params array.\n");

    rht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
             MaxSegWidth,MaxSegGap,MinDist,MaxDist,Accuracy,Threshold,
             NumOfTests,TextInfo,PostScript,ParametersOut);

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
  fprintf(stderr,"[-c NumOfTests] [-d MinDist MaxDist] [-t Threshold] ");
  fprintf(stderr,"[-a Accuracy] [-p] [-N Noise] \n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"Randomized Hough Transform (RHT) for lines.\n\n");
  print_common_options(stdout);
  fprintf(stdout,"\t-c NumOfTests\t\t\t\t\tdef: 1\n");
  fprintf(stdout,"\t   The number of tests to run.\n");
  fprintf(stdout,"\t-d MinDist MaxDist\t\t\t\tdef: 2 20\n");
  fprintf(stdout,"\t   Distance limits for randomly selected points.\n");
  fprintf(stdout,"\t-t Threshold\t\t\t\t\tdef: 2\n");
  fprintf(stdout,"\t   Minimum score accepted for a accumulator maximum.\n");
  fprintf(stdout,"\t-a Accuracy\t\t\t\t\tdef: 0.01\n");
  fprintf(stdout,"\t   The accumulator accuracy.\n");
  fprintf(stdout,"\t-p\t\t\t\t\t\tdef: not in use\n");
  fprintf(stdout,"\t   PostScript output (houghn.ps; n=1...N).\n");
  fprintf(stdout,"\t   PS files illustrate the dynamic accu representation.\n");
  fprintf(stdout,"\t-N Noise\t\t\t\t\tdef: 0.0\n");
  fprintf(stdout,"\t   Randomly added noise to test image (%%).\n");
  fprintf(stdout,"\t   Write the test image to ResultFile and exit.\n");
}
