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
 *	Duda, R.O. and Hart, P.E., Use of the Hough Transform To Detect Lines
 *	and Curves in Pictures, Communications of the ACM, vol. 15, no. 1,
 *	1972, pp. 11-15.
 *
 *	Risse, T., Hough Transform for Line Recognition: Complexity of
 *	Evidence Accumulation and Cluster Detection, Computer Vision,
 *	Graphics, and Image Processing, vol. 46, no. 3, 1989, pp. 327-345.
 *
 * File:    	sht_linedetect.c
 * Purpose: 	Standard Hough Transform (SHT) for lines
 * Date:    	Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif

pic_type acc_space[ACC_MAX_SIZE][ACC_MAX_SIZE],
         hough_kernel[3][3]={ {0,0,0},  /* global variable for  */
			      {0,1,0},  /* Hough-space addition */
			      {0,0,0}}; /* and subtraction      */

extern int StopFlag;

double maxs[MAX_LINES][2];

/*****************************************************************************/
sht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_the,dim_rho,real_params,NumOfMaxs,
         threshold,MinSegLen,MaxSegWidth,MaxSegGap,TextInfo,ParametersOut)
/*

  Implement the Standard Hough Transform (SHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y,dim_the,dim_rho - image and accu quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                threshold - minimum score accepted for an accumulator maximum
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                TextInfo - statistics output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE], pic1[][MAX_SIZE], gpic[][MAX_SIZE], rpic[][MAX_SIZE],
         threshold;
int dim_x,dim_y,dim_the,dim_rho,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,
    TextInfo;
double real_params[][2];
char ParametersOut[256];
{
  line *start=NULL;
  int i, j, NumOfRealMaxs, false_alarms=0;
  float totcputime=0.0;
  double rho_max;
  struct tms cputimebase, cputime;

  StopFlag=0;

  NotifyDispatch_MACRO();

  for (i=0; i<dim_the; i++)
    for (j=0; j<dim_rho; j++)
      acc_space[i][j]=(pic_type)0;

  times(&cputimebase); /* clocking starts */

  standard_hough_transform(pic,dim_x,dim_y,acc_space,dim_the,dim_rho);

  StopFlag=0;

  NotifyDispatch_MACRO();

  /* Draw Hough space */
  PutAccu_gray_MACRO(quads[3],acc_space,dim_the,dim_rho);
/*
  pic_to_disk("STH_acc_space.pgm",acc_space,dim_the,dim_rho,2,"test");
*/
  /* Find maximas */
  if (!StopFlag) NumOfRealMaxs=FindMaxs_rho_theta(&start,acc_space,NumOfMaxs,
                                       threshold,dim_the,dim_rho,pic,pic1,
                                       dim_x,dim_y,maxs,MinSegLen,
                                       MaxSegWidth,MaxSegGap,&false_alarms);

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
standard_hough_transform(pic,dim_x,dim_y,acc_space,dim_the,dim_rho)
/*
  Hough transformation from pic to acc_space, using line representation

          rho = x*cos(theta)+y*sin(theta).

  Original picture dimension is given by dim_x,dim_y and parameter
  space dimension is given by dim_the, dim_rho. Original picture is
  binary picture where pixel value OBJECT_PIX_VAL represent object,
  other values are background

  Parameters : pic - pic_type matrix, size [][MAX_SIZE]
               dim_x,dim_y - size of the original picture
               acc_space - pic_type matrix , size [][ACC_MAX_SIZE]
               dim_the, dim_rho - size of the accumulator space

  Remarks : Treats accumalator space as continous (trough modulus indexing)
            in rho coordinate direction. No testing for accumulator space
            overflow.

*/
pic_type pic[][MAX_SIZE],acc_space[][ACC_MAX_SIZE];
int dim_x,dim_y,dim_the,dim_rho;
{
  register i,j,k,l;
  double theta,rho,rho_min=-(double)dim_x,
         rho_max=sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));
/*
         rho_max=sqrt(2.0)*sqrt((double)(dim_x*dim_x)+(double)(dim_y*dim_y));
*/

  for (i=0; (i<dim_y) && !StopFlag; i++) {
    for (j=0; (j<dim_x) && !StopFlag; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
        PutPixel_MACRO(quads[1], j, i);
        for (k=0; (k<dim_the) && !StopFlag; k++) {
/*
          theta=(double)(-PI/2.0+k*PI/dim_the);
*/
          theta=(double)k*PI/(double)dim_the; /*theta[0, PI]*/
          rho=(double)j*cos(theta)+(double)i*sin(theta);
/*
          l=(int)rint(rho*dim_rho/rho_max/2.0+dim_rho/2.0);
*/
          l=nint(((rho-rho_min)/(rho_max-rho_min))*(double)(dim_rho-1));
          add_acc_space(hough_kernel,acc_space,dim_the,dim_rho,l,k);
        }
        NotifyDispatch_MACRO();
      }
    FlushDisp_MACRO();
  }

  return(0);
}


/*****************************************************************************/
add_acc_space(hough_kernel,acc_space,dim_theta,dim_rho,l,k)

/*
  Adds Hough accumulator space by specified value or values.
  Addition mask is given by hough_kernel and accumulator space position
  by l and k.

  Parameters : hough_kernel - pic_type matrix, size [][3]
               acc_space - pic_type matrix, size [][ACC_MAX_SIZE]
               dim_the, dim_rho - size of the accumulator space
               l,k - accumulator space position (matrix indices)

*/
pic_type hough_kernel[][3],acc_space[][ACC_MAX_SIZE];
int dim_theta,dim_rho,l,k;
{
  register i,j;
  int x,y;

  for (i= -1; i<=1; i++)
    for (j= -1; j<=1; j++) {
      x = k+j;
      y = l+i;
      if (x==-1) {
        x=dim_theta-1;
        y=dim_rho-1-y;
      }
      if (x==dim_theta) {
        x=0;
        y=dim_rho-1-y;
      }
      if (y>=0 && y<dim_rho)
        acc_space[y][x] += hough_kernel[i+1][j+1];
    }
}


/*****************************************************************************/
sub_acc_space(hough_kernel,acc_space,dim_theta,dim_rho,l,k)

/*
  Subtracts Hough accumulator space by specified value or values.
  Subtraction mask is given by hough_kernel and accumulator space position
  by l and k.

  Parameters : hough_kernel - pic_type matrix, size [][3]
               acc_space - pic_type matrix, size [][ACC_MAX_SIZE]
               dim_the, dim_rho - size of the accumulator space
               l,k - accumulator space position (matrix indices)

*/
pic_type hough_kernel[][3],acc_space[][ACC_MAX_SIZE];
int dim_theta,dim_rho,l,k;
{
  register i,j;
  int x,y;

  for (i= -1; i<=1; i++)
    for (j= -1; j<=1; j++) {
      x = k+j;
      y = l+i;
      if (x==-1) {
        x=dim_theta-1;
        y=dim_rho-1-y;
      }
      if (x==dim_theta) {
        x=0;
        y=dim_rho-1-y;
      }
      if (y>=0 && y<dim_rho) {
        acc_space[y][x] -= hough_kernel[i+1][j+1];
        if (acc_space[y][x]<0)
          acc_space[y][x]=0;
      }
    }
}
