/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki Kalviainen (Heikki.Kalviainen@lut.fi)
 *
 * File:         aht_lineserver.c
 * Purpose:      resolving options for the Adaptive Hough Transform (AHT)
 * Date:         Jun 1,1993
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
  int c, acc_space_sel=1, dim_x, dim_y, dim_m=9, dim_c=9, gray_dim_x, gray_dim_y,
      TextInfo=0, Shrink=0, NumOfMaxs=999, MinSegLen=10, MaxSegWidth=2,
      MaxSegGap=5,FileFormat=PGM_PIC, StoreMaximas=0;
  double Accur_m=0.01,Accur_c=1.0,Noise=0.0,bin_level=0.9;

  File[0]=GrayFile[0]=ResultFile[0]=RealParamsFile[0]='\0';
  cptr=rindex(argv[0],'/');
  strcpy(proc_name,cptr?++cptr:argv[0]);

  while ((c = getopt(argc, argv, "f:G:R:r:F:sn:l:w:g:m:Mihd:S:L:a:uP:")) != -1)
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
      case 'd': dim_m=atoi(optarg); 
                dim_c=atoi(argv[optind++]);
                break;
      case 'S': acc_space_sel=atoi(optarg); 
                break;
      case 'L': bin_level=atof(optarg); 
                break;
      case 'a': Accur_m=atof(optarg); 
                Accur_c=atof(argv[optind++]);
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
    copy_pic(pic,pic1,dim_x,dim_y);

    if (strlen(RealParamsFile))
      if (real_params_from_disk(RealParamsFile,real_params)!=real_params[0][0])
        warn_prnt("Something curious in real_params array.\n");

    aht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_m,dim_c,real_params,NumOfMaxs,
             MinSegLen,MaxSegWidth,MaxSegGap,acc_space_sel,bin_level,Accur_m,
             Accur_c,TextInfo,ParametersOut);

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
  fprintf(stderr,"[-d dim_m dim_c] [-S acc_space_sel] [-L bin_level] ");
  fprintf(stderr,"[-a Accur_m Accur_c]\n");
}

/*****************************************************************************/
void print_options()
{
  fprintf(stdout,"Adaptive Hough Transform (AHT) for lines.\n\n"); 
  print_common_options(stdout);
  fprintf(stdout,"\t-d dim_m dim_c\t\t\t\t\tdef: 9 9\n");
  fprintf(stdout,"\t   Size of the accumulator array.\n");
  fprintf(stdout,"\t-S acc_space_sel\t\t\t\tdef: 1\n");
  fprintf(stdout,"\t   Accumulator range selector [1,2].\n");
  fprintf(stdout,"\t   1 for horizontal lines (0-45 or 135-180 degrees),\n");
  fprintf(stdout,"\t   2 for vertical lines (45-135 degrees).\n");
  fprintf(stdout,"\t-L bin_level\t\t\t\t\tdef: 0.9\n");
  fprintf(stdout,"\t   Level to binarize accumulator with respect to the maximum value.\n");
  fprintf(stdout,"\t-a Accur_m Accur_c\t\t\t\tdef: 0.01 1.0\n");
  fprintf(stdout,"\t   Relative accuracy of the parameters with respect to quantization.\n");
}
