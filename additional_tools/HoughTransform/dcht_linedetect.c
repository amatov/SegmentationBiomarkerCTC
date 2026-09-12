/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing, Lappeenranta, Finland
 *
 * For further information, please notice following papers:
 *
 *	Leavers, V.F., Ben-Tzvi, D. and Sandler, M.B., A Dynamic Combinatorial
 *	Hough Transform for Straight Lines and Circles, Proceedings of 5th
 *	Alvey Vision Conference, Reading, UK, September 25-28, 1989,
 *	pp. 163-168.
 *
 *	Ben-Tzvi, D., Leavers, V.F. and Sandler, M.B., A Dynamic Combinatorial
 *	Hough Transform, Proceedings of the 5th International Conference on
 *	Image Analysis and Processing, Positano, Italy, September 20-22,
 *	1989, pp. 152-159.
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki Kalviainen (Heikki.Kalviainen@lut.fi)
 *
 * File:    	dcht_linedetect.c
 * Purpose: 	Dynamic Combinatorial Hough Transform (DCHT) for lines
 * Date:    	Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[], graph_quads[];
#endif

extern int StopFlag;

pic_type theta_histogram[ACC_MAX_SIZE];
int maxs_histogram[MAX_LINES];
double maxs[MAX_LINES][2];

pic_vec_type vec[MAX_VEC_COORD];

long rnd();

/*****************************************************************************/
dcht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_theta,real_params,NumOfMaxs,
          MinSegLen,MaxSegWidth,MaxSegGap,Threshold,NumOfTests,TextInfo,
	  ParametersOut)
/*

  The Dynamic Combinatorial Hough Transform (DCHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y,dim_theta - image and accu quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                Threshold - minimum score accepted for a accumulator peak
                NumOfTests - number of tests to run
                TextInfo - statistics output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE],pic1[][MAX_SIZE],gpic[][MAX_SIZE],rpic[][MAX_SIZE];
int dim_x,dim_y,dim_theta,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,Threshold,
    NumOfTests,TextInfo;
double real_params[][2];
char ParametersOut[256];
{
  line *start=NULL;
  int test, x1, x2, y1, y2, i, k, false_alarms=0, max_theta, theta_bin,
      found_params;
  long adv1, adv2, n, num_of_real_params=0L, acquired_real_params=0L;
  float totcputime=0.0;
  double a, b, rho, theta;
  struct tms cputimebase, cputime;

  StopFlag=0;

  NotifyDispatch_MACRO();

  init_random_generator();

  times(&cputimebase); /*clocking starts */

  for (test=0; test<NumOfTests && !StopFlag; test++) {

    copy_pic(pic,pic1,dim_x,dim_y);
    n=999;
    if (test) DeleteAllLineSegments(start);

    for (k=0; (k<NumOfMaxs) && (n>(2*MinSegLen)) && !StopFlag; ) {

      ClearQuadrant_MACRO(quads[1]);
      convert_pic_vec(pic,dim_x,dim_y,vec,&n);
      for (i=0; i<dim_theta; i++) /* zero the accumulator */
        theta_histogram[i]=(pic_type)0;

     /* choose the seed point randomly from binary picture */
      adv1=rnd(n-1);
      x1=vec[adv1].j;
      y1=vec[adv1].i;
      PutPixel_MACRO(quads[1], x1, y1);

      /* calculate parameters & accumulate parameter space */
      for (adv2=0; adv2<adv1 && !StopFlag; adv2++) { /* before the seed point */
        x2=vec[adv2].j;
        y2=vec[adv2].i;
        PutPixel_MACRO(quads[1], x2, y2);

        theta=atan2(1.0,(double)(y2-y1)/(double)(x1-x2));
        theta_bin=nint(theta/PI*(double)dim_theta);
        if (theta_bin>(dim_theta-1))
          theta_bin=0;
        theta_histogram[theta_bin]++;

        NotifyDispatch_MACRO();
      }
      FlushDisp_MACRO();

      for (adv2=adv1+1; adv2<n && !StopFlag; adv2++) { /* after the seed point */
        x2=vec[adv2].j;
        y2=vec[adv2].i;
        PutPixel_MACRO(quads[1], x2, y2);

        theta=atan2(1.0,(double)(y2-y1)/(double)(x1-x2));
        theta_bin=nint(theta/PI*(double)dim_theta);
        if (theta_bin>(dim_theta-1))
          theta_bin=0;
        theta_histogram[theta_bin]++;

        NotifyDispatch_MACRO();
      }
      FlushDisp_MACRO();

      for (i=max_theta=0; i<dim_theta && !StopFlag; i++) /* find maximum */
        if (theta_histogram[i]>max_theta) {
          max_theta=theta_histogram[i];
          theta_bin=i;
        }

      if (max_theta>=Threshold) {
          
        theta=(double)theta_bin*PI/(double)dim_theta;
        rho=x1*cos(theta)+y1*sin(theta);

        if (!fabs(sin(theta))) { /* (rho, theta) --> (a, b) */
          a=INF_DOUBLE_VAL;
          b=rint(rho);
        } else {
          a=-cos(theta)/sin(theta);
          if (isinf(a))
            b=rint(rho);
          else
            b=rho/sin(theta);
        }
        FollowAndMarkLineSegments_a_b(&start,a,b,pic,dim_x,dim_y,
                                      MinSegLen,MaxSegWidth,MaxSegGap);

        if (RemoveEdgePoints(start,pic,dim_x,dim_y)) {
          maxs[k][0]=rho;
          maxs[k][1]=theta;
          k++;
        } else
          false_alarms++;
      }
      pic[y1][x1]=(pic_type)0;

      NotifyDispatch_MACRO();
    }
    maxs_histogram[k]++;
    if ((int)real_params[0][0]>0) { /* check are the lines found 'real' lines */
      found_params=test_found_params_rho_theta(maxs,k,real_params);
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
    TextOut_dcht(maxs_histogram,test,totcputime,NumOfMaxs,
                 num_of_real_params,acquired_real_params,false_alarms);
  }

  store_all_line_segments(rpic,dim_x,dim_y,start,(pic_type)OBJECT_PIX_VAL);
  store_all_maxs_rho_theta(gpic,dim_x,dim_y,maxs,k,(pic_type)OBJECT_PIX_VAL,
	ParametersOut,start);

  DeleteAllLineSegments(start);
  StopFlag=0;
}
