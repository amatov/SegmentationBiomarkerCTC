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
 * The Randomized Hough Transform (RHT) algorithm is invented by Lei Xu, Erkki Oja, 
 * and Pekka Kultanen, and developed further by Heikki Kalviainen and Petri Hirvonen.
 * For further information, please notice following papers:
 *
 *          Kalviainen, H., Hirvonen, P., Xu, L., and Oja, E.,
 *          Probabilistic and Non-probabilistic Hough Transforms:
 *          Overview and Comparisons, vol. 13, no. 4, May 1995, pp. 239-252.
 *
 *          Xu, L., E. Oja, and Kultanen, P., A new curve detection method:
 *          Randomized Hough Transform (RHT),  Pattern Recognition Letters
 *          vol. 11, no. 5, 1990, pp. 331-338.
 *
 *          Kultanen, P., Xu, L., and  Oja, E., Randomized Hough Transform (RHT), 
 *          Proc. of the 10th International Conference on Pattern Recognition,
 *          Atlantic City, USA, June 16-21, 1990, pp. 631-635.
 *
 *          Xu, L. and Oja, E., Randomized Hough Transform (RHT): Basic
 *          Mechanisms, Algorithms, and Computational Complexities,  
 *          CVGIP: Image Understanding, vol. 57, no. 2, 1993, pp. 131-154.
 *
 *          Kalviainen, H., Hirvonen, P., Xu, L., and Oja, E.,
 *          Comparisons of Probabilistic and Non-probabilistic 
 *          Hough Transforms, 
 *          Proceedings of 3rd European Conference on Computer 
 *          Vision ECCV'94, Stockholm, Sweden, May 1994, pp. 351-360.
 *          
 *
 * File:    	drht_linedetect.c
 * Purpose: 	Dynamic Randomized Hough Transform (DRHT) for lines
 * Date:    	Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[], graph_quads[];
#endif

extern int StopFlag;

float f_size[MAX_TIME],f_ind[MAX_TIME];
int size[MAX_TIME],hits[MAX_LINES];
int maxs_histogram[MAX_LINES];
double maxs[MAX_LINES][2];

pic_vec_type vec[MAX_VEC_COORD];
pic_vec_type vec1[MAX_VEC_COORD]; /* could be smaller */

long rnd();

/*****************************************************************************/
drht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
          MaxSegWidth,MaxSegGap,MinDist,MaxDist,Accuracy,Accuracy1,Threshold,
          Threshold1,BlockWidth,MaxVar,NumOfTests,TextInfo,PostScript,
	  ParametersOut)
/*

  The Dynamic Randomized Hough Transform (DRHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y - image quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                MinDist,MaxDist - point distance limits
                Accuracy,Accuracy1 - accumulator resolutions for both iterations
                Threshold,Threshold1 - minimum score accepted for an accumulator
                                       maximum, both iterations
                BlockWidth - block width
                MaxVar - maximum variation of slope (a) in degrees
                NumOfTests - number of tests (multible tests for accurate times)
                TextInfo - statistics output selector
                PostScript - PostScript output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE],pic1[][MAX_SIZE],gpic[][MAX_SIZE],rpic[][MAX_SIZE];
int dim_x,dim_y,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,MinDist,MaxDist,
    Threshold,Threshold1,BlockWidth,NumOfTests,TextInfo,PostScript;
double real_params[][2],Accuracy,Accuracy1,MaxVar;
char ParametersOut[256];
{
  InfMat *acc_space;
  Accu *Cell;
  line *start=NULL;
  int test,x1,x2,y1,y2,maxval,i,k,NumOfAcc=0,max_size,false_alarms=0,found_params;
  long adv1, adv2, timecount=0L, n, n1, cumulative_size, num_of_real_params=0L,
       acquired_real_params=0L;
  float totcputime=0.0;
  double a,b,a1,b1,amin,amax,bmin,bmax,round=1./Accuracy,round1=1./Accuracy1;
  struct tms cputimebase,cputime;

  StopFlag=0;

  NotifyDispatch_MACRO();

  CreateInfMat(&acc_space);

  init_random_generator();

  times(&cputimebase); /* clocking starts */

  for (test=0; test<NumOfTests && !StopFlag; test++) {

    copy_pic(pic,pic1,dim_x,dim_y);
    convert_pic_vec(pic,dim_x,dim_y,vec,&n);
    timecount=0L; NumOfAcc=0;
    if (test) DeleteAllLineSegments(start);

    for (k=0; (k<NumOfMaxs) && (n>(2*MinSegLen)) && (NumOfAcc<=10) &&
              !StopFlag; ) {

      /* choose point pair randomly from binary picture */
      adv1=rnd(n-1);
      do {
        x1=vec[adv1].j;
        y1=vec[adv1].i;

        while (adv1==(adv2=rnd(n-1)));

        x2=vec[adv2].j;
        y2=vec[adv2].i; 
        adv1=adv2;

      } while (maxdist(x1,y1,x2,y2)>MaxDist || mindist(x1,y1,x2,y2)<MinDist);  

      /* calculate parameters & accumulate parameter space */

      a = (double)(y2-y1)/(double)(x2-x1);
      b = (double)y1-a*x1;

      /* Infinity value handling */
      if (!isinf(a))
        Cell=IncAccu(rint(round*a)/round,rint(round*b)/round,acc_space,1);
      else
        Cell=IncAccu(a,(double)x1,acc_space,1);

      /* show points which have been selected until now */

      timecount++;

      PutPixel_MACRO(quads[1], x1, y1);
      PutPixel_MACRO(quads[1], x2, y2);
      FlushDisp_MACRO();
      NotifyDispatch_MACRO();

      size[timecount]=ExamineAccu(acc_space,&amin,&amax,&bmin,&bmax,&maxval);

      PutLongText_MACRO(quads[1],0,250,timecount,!(timecount % 100));

      if (Cell->Accumulator>=Threshold) {

        if (size[timecount]==0)
          warn_prnt("Something curious in accumulator space.\n");

        DrawAccumulatorSpace_graph_MACRO(graph_quads[0],acc_space,
                                         "Accumulator","a","b",
                                         amin,amax,bmin,bmax,maxval);

        if (PostScript && k<1) 
          DrawAccumulatorSpacePS(acc_space,"Accumulator","a","b",
                                 amin,amax,bmin,bmax,maxval);

        RemoveAccumulatorSpace(acc_space);

        a=rint(round*a)/round;
        b=rint(round*b)/round;
        if (isinf(a)) b=(double)x1;

        ClearQuadrant_MACRO(quads[2]);

        /* Create a block around the line found and get blocked pixels */
            
	FollowAndGetBlockedPixels_Line(a,b,pic,dim_x,dim_y,BlockWidth,MaxVar,
                                       vec1,&n1);

        if (n1>1) { /* if not -> total false alarm */

          /* RHT using blocked pixels (2nd iteration) */
            
	  drht_line_2_iter(vec1,n1,MinDist,MaxDist,round1,Threshold1,
                           &a1,&b1,&timecount);

	  ShowDifferences(a,a1,b,b1);

          a=a1; /* new - hopefully more accurate - parameter values */
          b=b1;  
        }

        FollowAndMarkLineSegments_a_b(&start,a,b,pic,dim_x,dim_y,
                                      MinSegLen,MaxSegWidth,MaxSegGap);

        if (RemoveEdgePoints(start,pic,dim_x,dim_y)) {
          convert_pic_vec(pic,dim_x,dim_y,vec,&n);
          maxs[k][0]=a;
          maxs[k][1]=b;
          k++;
          hits[k]=timecount;
          NumOfAcc=0;
        } else
	  false_alarms++;

        CreateInfMat(&acc_space);

        ClearQuadrant_MACRO(quads[1]);
        PutLongText_MACRO(quads[1],0,250,timecount,TRUE);

        NumOfAcc++;
      }

      NotifyDispatch_MACRO();
    }
    maxs_histogram[k]++;
    if ((int)real_params[0][0]>0) {
      found_params=test_found_params(maxs,k,real_params);
      num_of_real_params+=found_params;
      if (found_params==NumOfMaxs)
        acquired_real_params++;
    }
  }

  times(&cputime); /*clocking ends */
  totcputime=1./60.*(cputime.tms_utime-cputimebase.tms_utime);

  if (TextInfo) {
    if ((int)real_params[0][0]<1)
      num_of_real_params=-1;
    TextOut_rht(size,timecount,maxs_histogram,test,totcputime,NumOfMaxs,
		num_of_real_params,acquired_real_params,false_alarms);
  }

  store_all_line_segments(rpic,dim_x,dim_y,start,(pic_type)OBJECT_PIX_VAL);
  store_all_maxs_line(gpic,dim_x,dim_y,maxs,k,(pic_type)OBJECT_PIX_VAL,
	ParametersOut,round,start);

  for(i=0; i<(timecount+1); i++) {
    f_size[i]=(float)size[i];
    f_ind[i]=(float)i;
  }
  DrawLineData_graph_MACRO(graph_quads[1],f_ind,f_size,timecount+1,
                           "Accum. size","picks","cells");
  if (PostScript)
    DrawLineDataPS(f_ind,f_size,timecount+1,"Accum. size","picks","cells"); 

  for(i=0; i<(k+1); i++) {
    f_size[i]=(float)hits[i];
    f_ind[i]=(float)i;
  }
  DrawLineData_graph_MACRO(graph_quads[2],f_size,f_ind,k+1,"      Hits",
                           "picks","maxs");
  if (PostScript)
    DrawLineDataPS(f_size,f_ind,k+1,"Hits","picks","maxs");

  DeleteAllLineSegments(start);
  RemoveAccumulatorSpace(acc_space);
  StopFlag=0;
}

/*****************************************************************************/

drht_line_2_iter(vec1,n1,min_dist,max_dist,round1,Threshold1,a1,b1,time) 
/*

  Dynamic Randomized Hough Transform (DRHT) for line detection, second iteration

  Parameters :  vec1 - edge point vector
                n - vector length
                min_dist,max_dist - point distance limits
                round1 - accumulator resolutions for 2nd iteration
                Threshold1 - minimum score accepted for an accumulator peak
                a1 - new slope
                b1 - new intercept
                time - timecount

*/
pic_vec_type vec1[];
long n1, *time;
int min_dist, max_dist, Threshold1;
double *a1, *b1, round1;
{
  InfMat *acc_space1;
  Accu *Cell;
  int x1,x2,y1,y2;
  long adv1, adv2, timer=0L;   
  int LineFound=FALSE,max_val;
  double a_min, a_max, b_min, b_max;

  CreateInfMat(&acc_space1);

  while (!LineFound) {
    adv1=rnd(n1-1);
    do {
      x1=vec1[adv1].j;
      y1=vec1[adv1].i;

      while (adv1==(adv2=rnd(n1-1)));

      x2=vec1[adv2].j;
      y2=vec1[adv2].i; 
      adv1=adv2;

    } while (maxdist(x1,y1,x2,y2)>max_dist || mindist(x1,y1,x2,y2)<min_dist);

    PutPixel_MACRO(quads[2], x1, y1);
    PutPixel_MACRO(quads[2], x2, y2);

    (*a1) = (double)(y2-y1)/(double)(x2-x1);
    (*b1) = (double)y1-(*a1)*x1;

    if (!isinf(*a1))
      Cell=IncAccu(rint(round1*(*a1))/round1,rint(round1*(*b1))/round1,
                   acc_space1,1);
    else
      Cell=IncAccu(*a1,(double)x1,acc_space1,1);

    (*time)++; timer++;

    size[*time]=ExamineAccu(acc_space1,&a_min,&a_max,&b_min,&b_max,&max_val);

    PutLongText_MACRO(quads[2],0,11,timer,!(timer % 10));

    if (Cell->Accumulator>=Threshold1)
      LineFound=TRUE;

    NotifyDispatch_MACRO();
  }

  (*a1)=rint(round1*(*a1))/round1; 
  (*b1)=rint(round1*(*b1))/round1;	    
  if (isinf(*a1)) (*b1)=(double)x1;

  PutLongText_MACRO(quads[2],0,11,timer,TRUE);

  RemoveAccumulatorSpace(acc_space1);
}


/*****************************************************************************/

ShowDifferences(a,a1,b,b1)
/*

  Show old and new lines

  Parameters :  a - old slope
                b - old intercept
                a1 - new slope
                b1 - new intercept

*/
double a,a1,b,b1;
{
  int x,y;

  if (isinf(a)) {
    x = (int)rint(b);
    for (y=5;y<250;y++)
      if (x>5 && x<250)
        PutPixel_gray_MACRO(quads[2],x,y,MID_GRAY_LEVEL);
  }
  if (fabs(a)<0.001) { 
    y = (int)rint(b);
    for (x=5;x<250;x++)
      if (y>5 && y<250)
        PutPixel_gray_MACRO(quads[2],x,y,MID_GRAY_LEVEL);
  } else
    for (x=5;x<250;x++) {
      y = (int)rint(a*x+b);
      if (y>5 && y<250)
        PutPixel_gray_MACRO(quads[2],x,y,MID_GRAY_LEVEL);
    }

  if (isinf(a1)) {
    x = (int)rint(b1);
    for (y=5;y<250;y++)
      if (x>30 && x<200)
        PutPixel_MACRO(quads[2],x,y);
  }
  if (fabs(a1)<0.001) {
    y = (int)rint(b1);
    for (x=5;x<250;x++)
      if (y>30 && y<200)
        PutPixel_MACRO(quads[2],x,y);
  } else
    for (x=5;x<250;x++) {
      y = (int)rint(a1*x+b1);
      if (y>30 && y<200)
        PutPixel_MACRO(quads[2],x,y);
  }
}
