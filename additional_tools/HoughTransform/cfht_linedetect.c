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
 *	Liang, P., A New Transform for Curve Detection, Proceedings of Third
 *	International Conference on Computer Vision, Osaka, Japan, December
 *	4-7, 1990, pp. 748-751.
 *
 *	Liang, P., A New and Efficient Transform for Curve Detection, Journal
 *	of Robotic Systems, vol. 8, no. 6, 1991, pp. 841-847.
 *
 * File:         cfht_linedetect.c
 * Purpose:      the Curve Fitting Hough Transform (CFHT) for lines
 * Date:         Jun 1, 1993
 * Last change:  Oct 10,1995
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#define sqr(a) ((a)*(a))

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[], graph_quads[];
#endif

extern int StopFlag;

long curve_fitting_hough_transform();
double least_square_line_fit();

int size[MAX_TIME];
float f_size[MAX_TIME],f_ind[MAX_TIME];
double maxs[MAX_LINES][2];

/*****************************************************************************/
cfht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
          MaxSegWidth,MaxSegGap,Msize,ThreshForPoints,TolForFitting,Threshold,
          Epsilon,Gamma,DontRemovePoints,NumOfTests,TextInfo,PostScript,
	  ParametersOut)
/*

  The Curve Fitting Hough Transform (CFHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y - image quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                Msize - fitting window size (2*Msize+1)
                ThreshForPoints - minimum number of points in the window
                TolForFitting - maximum fitting error accepted
                Threshold - minimum score accepted for an accumulator maximum
                Epsilon - dist. between two accumulator cells when merging them 
                Gamma - score weight
                DontRemovePoints - any image points are not removed
                NumOfTests - number of tests (multible tests for accurate times)
                TextInfo - statistics output selector
                PostScript - PostScript output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE],pic1[][MAX_SIZE],gpic[][MAX_SIZE],rpic[][MAX_SIZE];
int dim_x,dim_y,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,Msize,
    ThreshForPoints,Threshold,DontRemovePoints,NumOfTests,TextInfo,PostScript;
double real_params[][2], TolForFitting, Epsilon, Gamma;
char ParametersOut[256];
{
  InfMat *acc_space;
  line *start=NULL;
  int i,accu_maxs,test,false_alarms;
  long timecount=0L;
  float totcputime=0L;
  struct tms cputimebase, cputime;

  StopFlag=0;

  NotifyDispatch_MACRO();

  times(&cputimebase); /*clocking starts */

  /* multible tests because of small run times */
  for (test=0; test<NumOfTests && !StopFlag; test++) {

    false_alarms=0;
    copy_pic(pic,pic1,dim_x,dim_y);
    if (test) {
      DeleteAllLineSegments(start);
      start=(line *)NULL;
      RemoveAccumulatorSpace(acc_space);
    }
    CreateInfMat(&acc_space);

    timecount=curve_fitting_hough_transform(acc_space,size,pic,dim_x,dim_y,
                                            Msize,ThreshForPoints,
                                            TolForFitting,Epsilon,Gamma,
                                            DontRemovePoints);

    /* what we have in the accumulator now */
    if (size[timecount]==0)
      warn_prnt("Something curious in accumulator space.\n");

    /* lets get all suitable parameters from the accumulator */

    accu_maxs=ReadAccuAndMarkLineSegments(&start,acc_space,maxs,Threshold,
                                          NumOfMaxs,pic,dim_x,dim_y,MinSegLen,
                                          MaxSegWidth,MaxSegGap,&false_alarms);
  }

  times(&cputime); /* clocking ends */
  totcputime=1./60.*(cputime.tms_utime-cputimebase.tms_utime);

  if (TextInfo)
    TextOut_cfht(size[timecount],accu_maxs,test,totcputime,maxs,real_params,
		 false_alarms);

  store_all_line_segments(rpic,dim_x,dim_y,start,(pic_type)OBJECT_PIX_VAL);
  store_all_maxs_line(gpic,dim_x,dim_y,maxs,accu_maxs,(pic_type)OBJECT_PIX_VAL,
	ParametersOut,1.0,start);

  for(i=0; i<(timecount+1); i++) {
    f_size[i]=(float)size[i];
    f_ind[i]=(float)i;
  }
  DrawLineData_graph_MACRO(graph_quads[1],f_ind,f_size,timecount+1,
                           "Accum. size","edges","cells");
  if (PostScript)
    DrawLineDataPS(f_ind,f_size,timecount+1,"Accum. size","edges","cells"); 

  DeleteAllLineSegments(start);
  RemoveAccumulatorSpace(acc_space);
  StopFlag=0;
}

/*****************************************************************************/
long curve_fitting_hough_transform(acc_space,size,pic,dim_x,dim_y,Msize,
                                   ThreshForPoints,TolForFitting,Epsilon,Gamma,
                                   DontRemovePoints)
/*

  Curve Fitting Hough Transform (CFHT) for line detection.

  Parameters :  acc_space - dynamic accumulator space
                size - accumulaotor size
                pic - input image
                dim_x,dim_y - image quantization
                Msize - fitting window size (2*Msize+1)
                ThreshForPoints - minimum number of points in the window
                TolForFitting - maximum fitting error accepted
                Epsilon - minimum dist. between two accumulator cells
                Gamma - score weight
                DontRemovePoints - any image points are not removed

*/
InfMat *acc_space;
int size[];
pic_type pic[][MAX_SIZE];
int dim_x, dim_y, Msize, ThreshForPoints, DontRemovePoints;
double TolForFitting, Epsilon, Gamma;
{
  int i,j,maxval;
  double a,b,amin,amax,bmin,bmax;
  long timecount=0L;

  /* scan the image (only edge points!!!) */
  for(i=0; i<dim_y && !StopFlag; i++)
    for(j=0; j<dim_x; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL &&
        least_square_line_fit(pic,j,i,Msize,ThreshForPoints,
                              dim_x,dim_y,&a,&b,DontRemovePoints)<TolForFitting) {

        /* accumulate parameter space */

        IncAccu_cfht(a,b,acc_space,Epsilon,Gamma,1);
        timecount++;

	size[timecount]=ExamineAccu(acc_space,&amin,&amax,&bmin,&bmax,&maxval);

        PutLongText_MACRO(quads[1],0,250,timecount,!(timecount % 100));
        DrawAccumulatorSpace_graph_MACRO(graph_quads[0],acc_space,"Accumulator",
                                         "a","b",amin,amax,bmin,bmax,maxval);

        FlushDisp_MACRO();
        NotifyDispatch_MACRO();
      }

  return timecount;
}

/*****************************************************************************/
double least_square_line_fit(pic,pj,pi,fit_block,num_of_points,dim_x,dim_y,a,b,
                             DontRemovePoints)
/*
  Curve fitting by the least square method.

  Parameters :	pic - edge image
                pj,pi - window center point
                fit_block - fitting window size (2*fit_block+1)
                num_of_points - minimum number of points in the window
                dim_x,dim_y - image quantization
                a,b - line parameters
                DontRemovePoints - any image points are not removed

*/
pic_type pic[][MAX_SIZE];
int pi, pj, fit_block, num_of_points, dim_x, dim_y, DontRemovePoints;
double *a, *b;
{
  int i, j, one_j,min_i=pi-fit_block,min_j=pj-fit_block,
                  max_i=pi+fit_block+1,max_j=pj+fit_block+1;
  double sx=0.0,sy=0.0,st2=0.0,ss=0.0,wt=1.0,sxoss,t,se=0.0;

  if (min_i<0) min_i=0;
  if (min_j<0) min_j=0;
  if (max_i>dim_y) max_i=dim_y;
  if (max_j>dim_x) max_j=dim_x;

  for(i=min_i; i<max_i; i++)
    for(j=min_j; j<max_j; j++) {
      PutPixel_gray_MACRO(quads[1],j,i,MID_GRAY_LEVEL);
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
        ss+=wt;
        sx+=(j*wt);
        sy+=(i*wt);
        PutPixel_MACRO(quads[1],j,i);
      }
    }

  if (ss<num_of_points) { /* is there enough points in the block for line fitting? */
    if (!DontRemovePoints) {
      pic[pi][pj]=(pic_type)0;
      PutPixel_gray_MACRO(quads[0],pj,pi,(pic_type)0);
    }
    return MAXDOUBLE;     /* if no, lets return far too big fitting error */
  }

  sxoss=sx/ss;
  *a=0.0;

  for(i=min_i; i<max_i; i++)
    for(j=min_j; j<max_j; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
        one_j=j;
        t=j-sxoss;
        st2+=t*t;
        *a+=t*i;
      }

  if (st2==0.0) {
    *a=INF_DOUBLE_VAL; /* line x = x0 -> Infinite value for a */
    *b=(double)one_j;
    return 0.0; /* no fitting error in this case */
  } else {
    *a/=st2;
    *b=(sy-sx*(*a))/ss;
  }

  for(i=min_i; i<max_i; i++) /* fitting error calculation */
    for(j=min_j; j<max_j; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL)
        se+=sqr((*a)*j-i+(*b))/((*a)*(*a)+1);

  return se/ss;
}
