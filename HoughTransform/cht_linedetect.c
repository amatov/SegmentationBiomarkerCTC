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
 * For further information, please notice following paper:
 *
 *	Ben-Tzvi, D. and Sandler, M.B., A Combinatorial Hough Transform,
 *	Pattern Recognition Letters, vol. 11, no. 3, 1990, pp. 167-174.
 *
 * File:         cht_linedetect.c
 * Purpose:      Combinatorial Hough Transform (CHT) for lines
 * Date:         Jun 1, 1993
 * Last change:  Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif

pic_type acc_space[ACC_MAX_SIZE][ACC_MAX_SIZE];

extern int StopFlag;

double maxs[MAX_LINES][2];

/*****************************************************************************/

cht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_the,dim_rho,real_params,NumOfMaxs,
         threshold,segments,overlap,mask_size,MinSegLen,MaxSegWidth,MaxSegGap,
         TextInfo,ParametersOut)
/*

  The Combinatorial Hough Transform (CHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y,dim_the,dim_rho - image and accu quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                threshold - minimum score accepted for an accumulator peak
                segments - numer of image segment in one dimension
                overlap - overlapped pixels
                mask_size - peak detection mask size
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                TextInfo - statistics output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE], pic1[][MAX_SIZE], gpic[][MAX_SIZE], rpic[][MAX_SIZE],
         threshold;
int dim_x,dim_y,dim_the,dim_rho,NumOfMaxs,segments,overlap,mask_size,MinSegLen,
    MaxSegWidth,MaxSegGap,TextInfo;
double real_params[][2];
char ParametersOut[256];
{
  line *start=NULL;
  int i, j, NumOfRealMaxs, false_alarms=0;
  float totcputime=0.0;
  struct tms cputimebase, cputime;
  StopFlag=0;

  NotifyDispatch_MACRO();

  for (i=0; i<dim_the; i++)
    for (j=0; j<dim_rho; j++)
      acc_space[i][j]=(pic_type)0;

  times(&cputimebase); /* clocking starts */

  /* the transform itself */
  combinatorial_hough_transform(pic,dim_x,dim_y,acc_space,dim_the,dim_rho,
                                segments,overlap);

  StopFlag=0;

  NotifyDispatch_MACRO();

  /* Draw Hough space */
  PutAccu_gray_MACRO(quads[3],acc_space,dim_the,dim_rho);
/*
  pic_to_disk("CHT_acc_space.pgm",acc_space,dim_the,dim_rho,2,"test");
*/

  /* Find maximas */
  if (!StopFlag) NumOfRealMaxs=FindMaxs_rho_theta_one(&start,acc_space,
                                       NumOfMaxs,threshold,dim_the,dim_rho,pic,
                                       dim_x,dim_y,maxs,MinSegLen,
                                       MaxSegWidth,MaxSegGap,&false_alarms,
                                       mask_size);

  times(&cputime); /* clocking ends */
  totcputime=1./60.*(cputime.tms_utime-cputimebase.tms_utime);

  NotifyDispatch_MACRO();

  if (TextInfo)
    TextOut_sht(NumOfRealMaxs,totcputime,maxs,real_params,false_alarms);

  store_all_line_segments(rpic,dim_x,dim_y,start,(pic_type)OBJECT_PIX_VAL);
  store_all_maxs_rho_theta(gpic,dim_x,dim_y,maxs,NumOfRealMaxs,
                           (pic_type)OBJECT_PIX_VAL,ParametersOut,start);

  StopFlag=0;
  NotifyDispatch_MACRO();
}

/*****************************************************************************/
combinatorial_hough_transform(pic,dim_x,dim_y,acc_space,dim_the,dim_rho,
                              segments,overlap)
/*
  Hough transformation from pic to acc_space, using line representation

          rho = x*cos(theta)+y*sin(theta).

  Original picture dimension is given by dim_x,dim_y and parameter
  space dimension is given by dim_the, dim_rho. Original picture is
  binary picture where pixel value 255 represent object, other values
  are background

  Parameters : pic - pic_type matrix, size [][MAX_SIZE]
               dim_x,dim_y - dimensions of the original picture
               acc_space - pic_type matrix , size [][ACC_MAX_SIZE]
               dim_the, dim_rho - dimensions of the accumulator space
               segments - number of image segments in one dimension
               overlap - number of overlapping pixels in each dimension

  Remarks : Treats accumalator space as continous (trough modulus indexing)
            in rho coordinate direction. No testing for accumulator space
            overflow.

*/


pic_type pic[][MAX_SIZE],acc_space[][ACC_MAX_SIZE];
int dim_x,dim_y,dim_the,dim_rho,segments,overlap;
{
  register i,j,k,l,m,n,n1,x_min,x_max,y_min,y_max;
  int rho_bin,theta_bin;
  double theta,rho,rho_min=-(double)dim_x,
         rho_max=sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));

  for (i=0; (i<segments) && !StopFlag; i++) { /* segment image */
    y_min=nint((double)i*((double)dim_y/(double)segments)-         /* divide */
               (double)overlap/2.0);                 /* overlapped pixels to */
    if (y_min<0)                                /* both sides of the segment */
      y_min=0;
    y_max=nint((double)(i+1)*((double)dim_y/(double)segments)+
               (double)overlap/2.0)+1;
    if (y_min>dim_y)
      y_min=dim_y;
    for (j=0; (j<segments) && !StopFlag; j++) {
      x_min=nint((double)j*((double)dim_x/(double)segments)-
                 (double)overlap/2.0);
      if (x_min<0)
        x_min=0;
      x_max=nint((double)(j+1)*((double)dim_x/(double)segments)+
               (double)overlap/2.0)+1;
      if (x_min>dim_x)
        x_min=dim_x;
      for (k=y_min; (k<y_max) && !StopFlag; k++) /* scan the segment */
        for (l=x_min; (l<x_max) && !StopFlag; l++) {
/*
          PutPixel_gray_MACRO(quads[1], l, k, MID_GRAY_LEVEL);
*/
          if (pic[k][l]==(pic_type)255) {
            n1=l+1;
            PutPixel_MACRO(quads[1], l, k);
            for (m=k; (m<y_max) && !StopFlag; m++) {
              for (n=n1; (n<x_max) && !StopFlag; n++)
                if (pic[m][n]==(pic_type)255) {
                  theta=atan2(1.0,(double)(m-k)/(double)(l-n));
                  theta_bin=nint(theta/PI*(double)dim_the);
                  if (theta_bin>(dim_the-1)) theta_bin=0;
                  rho=(double)l*cos(theta)+(double)k*sin(theta);
                  rho_bin=nint(((rho-rho_min)/(rho_max-rho_min))*
			       (double)(dim_rho-1));
                  acc_space[rho_bin][theta_bin]++; /* accmulate */
                }
              n1=x_min;
            }
            NotifyDispatch_MACRO();
          }
	}
      FlushDisp_MACRO();
    }
  }

  return(0);
}
