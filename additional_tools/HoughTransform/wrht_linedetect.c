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
 * File:    	wrht_linedetect.c
 * Purpose: 	Window Randomized Hough Transform (WRHT) for lines
 *              (including the Connective Randomized Hough Transform (CRHT))
 * Date:    	Jun 1, 1993
 * Last change: Oct 10, 1995
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#define sqr(a) ((a)*(a))

long rnd();
double WindowLeastSqrLineFit();

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
pic_type window[MAX_FIT_WIN_SIZE][MAX_FIT_WIN_SIZE];

/*****************************************************************************/
wrht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
          MaxSegWidth,MaxSegGap,Accuracy,Threshold,ThreshForPoints,
          TolForFitting,WinSize,DontRmAccu,ConnectiveCheck,Sectoring,
          NumOfTests,TextInfo,PostScript,ParametersOut)
/*

  The Window RHT (WRHT) including the Connective RHT (CRHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y - image quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                Accuracy - accmultor resolution
                Threshold - minimum score accepted for an accumulator maximum
                ThreshForPoints - minimum number of points in the window
                TolForFitting - maximum fitting error accepted
                WinSize - fitting window size (2*Msize+1)
                DontRmAccu - accumulator is not removed after line is found
                ConnectiveCheck - do the connective component check for window
                Sectoring - do the connective component check as sectored
                NumOfTests - number of tests
                TextInfo - statistics output selector
                PostScript - PostScript output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE],pic1[][MAX_SIZE],gpic[][MAX_SIZE],rpic[][MAX_SIZE];
int dim_x,dim_y,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,Threshold,
    ThreshForPoints,WinSize,DontRmAccu,ConnectiveCheck,Sectoring,NumOfTests,
    TextInfo,PostScript;
double real_params[][2],Accuracy,TolForFitting;
char ParametersOut[256];
{
  InfMat *acc_space;
  Accu *Cell;
  line *start=NULL;
  int test,x1,y1,maxval,i,k,NumOfAcc=0,LineFound=FALSE,WindowFits,
      max_size,false_alarms=0,found_params;
  long adv1, timecount=0L, n, cumulative_size, num_of_real_params=0L,
       acquired_real_params=0L;
  float totcputime=0.0;
  double a,b,amin,amax,bmin,bmax,round=1./Accuracy;
  struct tms cputimebase, cputime;

  StopFlag=0;

  NotifyDispatch_MACRO();

  init_random_generator();

  times(&cputimebase); /*clocking starts */

  for (test=0; test<NumOfTests && !StopFlag; test++) {

    if (Threshold>1) {
      if (test) RemoveAccumulatorSpace(acc_space);
      CreateInfMat(&acc_space);
    }

    copy_pic(pic,pic1,dim_x,dim_y);
    convert_pic_vec(pic,dim_x,dim_y,vec,&n);
    timecount=0L; NumOfAcc=0;
    if (test) DeleteAllLineSegments(start);

    for (k=0; (k<NumOfMaxs) && (n>(MinSegLen*2)) && (NumOfAcc<=10) &&
              !StopFlag; ) {

      /* choose one point randomly from binary picture */

      WindowFits=0;
      do {
        adv1=rnd(n-1);
        x1=vec[adv1].j;
        y1=vec[adv1].i;
        WindowFits++;
        NotifyDispatch_MACRO();

      } while (WindowLeastSqrLineFit(pic,x1,y1,window,WinSize,ThreshForPoints,
                                     dim_x,dim_y,&a,&b,ConnectiveCheck,
                                     Sectoring)>TolForFitting &&
               WindowFits<n && !StopFlag);

      timecount++;

      if (WindowFits>=n || StopFlag) { /* no more lines in image ? */
        size[timecount]=0;
        break;
      }

      if (Threshold>1) {

        /* Infinity value handling */

        if (!isinf(a))
          Cell=IncAccu(rint(round*a)/round,rint(round*b)/round,
                       acc_space,1);
        else
          Cell=IncAccu(a,b,acc_space,1);

        size[timecount]=ExamineAccu(acc_space,&amin,&amax,&bmin,&bmax,&maxval);

        if (Cell->Accumulator>=Threshold) {

          LineFound=TRUE;

          if (size[timecount]==0)
            warn_prnt("Something curious in accumulator space.\n");

          DrawAccumulatorSpace_graph_MACRO(graph_quads[0],acc_space,
                                           "Accumulator","a","b",
                                           amin,amax,bmin,bmax,maxval);
          if (PostScript && k<1) 
            DrawAccumulatorSpacePS(acc_space,"Accumulator","a","b",
                                   amin,amax,bmin,bmax,maxval);

          if (DontRmAccu)
            Cell->Accumulator=0; /* 'Removes' only one cell */
          else {
	    RemoveAccumulatorSpace(acc_space);
            CreateInfMat(&acc_space);
          }
        }

      } else {
        size[timecount]=0;
        LineFound=TRUE;
      }

      if (LineFound) {

        a=rint(round*a)/round;
        b=rint(round*b)/round;

        FollowAndMarkLineSegments_a_b(&start,a,b,pic,dim_x,dim_y,MinSegLen,
                                      MaxSegWidth,MaxSegGap);

        if (RemoveEdgePoints(start,pic,dim_x,dim_y)) {
          convert_pic_vec(pic,dim_x,dim_y,vec,&n);
          maxs[k][0]=a;
          maxs[k][1]=b;
          k++;
          hits[k]=timecount;
          NumOfAcc=0;
        } else
	  false_alarms++;


        ClearQuadrant_MACRO(quads[1]);
        PutLongText_MACRO(quads[1],0,250,timecount,TRUE);

        NumOfAcc++;
        LineFound=FALSE;
      }

      NotifyDispatch_MACRO();
    }
    maxs_histogram[k]++;
    if ((int)real_params[0][0]>0) { /* real params testing */
      found_params=test_found_params(maxs,k,real_params);
      num_of_real_params+=found_params;
      if (found_params==NumOfMaxs)
        acquired_real_params++;
    }
  }

  times(&cputime);
  totcputime=1./60.*(cputime.tms_utime-cputimebase.tms_utime);

  if (TextInfo) {
    if (real_params[0][0]<1)
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
  if (Threshold>1)
    RemoveAccumulatorSpace(acc_space);
  StopFlag=0;
}

/*****************************************************************************/
double WindowLeastSqrLineFit(pic,pj,pi,window,fit_window,num_of_points,
                             dim_x,dim_y,a,b,connective_check,sectoring)
/*
  Curve fitting by the least square method.

  Parameters :  pic - edge image
                pj,pi - window center point
                window - window array
                fit_window - fitting window size (2*fit_block+1)
                num_of_points - minimum number of points in the window
                dim_x,dim_y - image quantization
                a,b - line parameters
                connective_check - do the connective component check for window
                sectoring - do the connective component check as sectored

*/
pic_type pic[][MAX_SIZE],window[][MAX_FIT_WIN_SIZE];
int pi, pj, fit_window, num_of_points, dim_x, dim_y, connective_check,
    sectoring;
double *a, *b;
{
  int i,j,min_i=pi-fit_window,min_j=pj-fit_window,
          max_i=pi+fit_window+1,max_j=pj+fit_window+1,win_ci,win_cj;
  double sx=0.0,sy=0.0,st2=0.0,ss=0.0,wt=1.0,sxoss,t,se=0.0;

  if (min_i<0) min_i=0;                 /* i=y */
  if (min_j<0) min_j=0;                 /* j=x */
  if (max_i>dim_y) max_i=dim_y;
  if (max_j>dim_x) max_j=dim_x;
  win_ci=pi-min_i;
  win_cj=pj-min_j;

  if (connective_check) { /* use only center point connected points */
  fit_window=fit_window*2+1;
  for(i=0; i<fit_window; i++)
    for(j=0; j<fit_window; j++)
      window[i][j]=(pic_type)0; /* zero the window */
    window[win_ci][win_cj]=(pic_type)OBJECT_PIX_VAL; /* obviously an edge */
    get_connected_neighbors(window,win_ci,win_cj,pic,pi,pj,min_i+1,min_j+1,
                            max_i-2,max_j-2,sectoring?0:-1);
    for(i=min_i; i<max_i; i++)
      for(j=min_j; j<max_j; j++) {
        PutPixel_gray_MACRO(quads[1],j,i,MID_GRAY_LEVEL);
        if (window[i-min_i][j-min_j]==(pic_type)OBJECT_PIX_VAL) {
          ss+=wt;
          sx+=(j*wt);
          sy+=(i*wt);
          PutPixel_MACRO(quads[1],j,i);
        }
      }

  } else /* use all windowed points */
    for(i=min_i; i<max_i; i++)
      for(j=min_j; j<max_j; j++) {
        PutPixel_gray_MACRO(quads[1],j,i,MID_GRAY_LEVEL);
        if ((window[i-min_i][j-min_j]=pic[i][j])==(pic_type)OBJECT_PIX_VAL) {
          ss+=wt;
          sx+=(j*wt);
          sy+=(i*wt);
          PutPixel_MACRO(quads[1],j,i);
        }
      }

  if (ss<num_of_points) /* is there enough points in window for line fit ? */
    return MAXDOUBLE;   /* if not, lets return far too big fitting error */

  sxoss=sx/ss;
  *a=0.0;

  for(i=min_i; i<max_i; i++)
    for(j=min_j; j<max_j; j++)
      if (window[i-min_i][j-min_j]==(pic_type)OBJECT_PIX_VAL) {
        t=j-sxoss;
        st2+=t*t;
        *a+=t*i;
      }

  if (st2==0.0) {
    *a=INF_DOUBLE_VAL; /* line x = x0 -> Infinite value for 'a' */
    *b=(double)pj;
    return 0.0; /* no fitting error in this case */
  } else {
    *a/=st2;
    *b=(sy-sx*(*a))/ss;
  }

  for(i=min_i; i<max_i; i++)
    for(j=min_j; j<max_j; j++)
      if (window[i-min_i][j-min_j]==(pic_type)OBJECT_PIX_VAL)
        se+=sqr((*a)*j-i+(*b))/((*a)*(*a)+1);

  return se/ss; /* mean square line fitting error */
}

/*****************************************************************************/
get_connected_neighbors(win,win_ci,win_cj,pic,pic_ci,pic_cj,min_i,min_j,
                        max_i,max_j,sector)
/*
  Curve fitting by the least square method.

  Parameters :  win - window array
                win_ci,win_cj - window center point in window
                pic - edge image
                pic_ci,pic_cj - window center point in image
                fit_window - fitting window size (2*fit_block+1)
                min_i,min_j,max_i,max_j - window borders in image
                sector- current sector (-1: no sectoring; 0: all sectors)

*/
pic_type pic[][MAX_SIZE],win[][MAX_FIT_WIN_SIZE];
int win_ci,win_cj,pic_ci,pic_cj,min_i,min_j,max_i,max_j,sector;
{
  int checked_sectors=0;

  if (pic_ci<min_i || pic_ci>max_i || pic_cj<min_j || pic_cj>max_j ||
      win_ci>(MAX_FIT_WIN_SIZE-2) || win_cj>(MAX_FIT_WIN_SIZE-2))
    return; /* the window edge (vertex) has been reached */
                                             /*            8 7 6  */
  switch (sector) {                          /* Sectors:   1 0 5  */
    case -1:/* no sectoring */               /*            2 3 4  */
    case 0: checked_sectors=-5; /* starting from the center point */
    case 1: if (pic[pic_ci+1][pic_cj-1] && !win[win_ci+1][win_cj-1]) {
              win[win_ci+1][win_cj-1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci+1,win_cj-1,pic,pic_ci+1,pic_cj-1,
                                      min_i,min_j,max_i,max_j,!sector?8:sector);
            }
            checked_sectors++;
    case 2: if (pic[pic_ci][pic_cj-1] && !win[win_ci][win_cj-1]) {
              win[win_ci][win_cj-1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci,win_cj-1,pic,pic_ci,pic_cj-1,
                                      min_i,min_j,max_i,max_j,!sector?1:sector);
            }
            checked_sectors++;
    case 3: if (pic[pic_ci-1][pic_cj-1] && !win[win_ci-1][win_cj-1]) {
              win[win_ci-1][win_cj-1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci-1,win_cj-1,pic,pic_ci-1,pic_cj-1,
                                      min_i,min_j,max_i,max_j,!sector?2:sector);
            }
            if (++checked_sectors==3) return;
    case 4: if (pic[pic_ci-1][pic_cj] && !win[win_ci-1][win_cj]) {
              win[win_ci-1][win_cj]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci-1,win_cj,pic,pic_ci-1,pic_cj,
                                      min_i,min_j,max_i,max_j,!sector?3:sector);
            }
            if (++checked_sectors==3) return;
    case 5: if (pic[pic_ci-1][pic_cj+1] && !win[win_ci-1][win_cj+1]) {
              win[win_ci-1][win_cj+1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci-1,win_cj+1,pic,pic_ci-1,pic_cj+1,
                                      min_i,min_j,max_i,max_j,!sector?4:sector);
            }
            if (++checked_sectors==3) return;
    case 6: if (pic[pic_ci][pic_cj+1] && !win[win_ci][win_cj+1]) {
              win[win_ci][win_cj+1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci,win_cj+1,pic,pic_ci,pic_cj+1,
                                      min_i,min_j,max_i,max_j,!sector?5:sector);
            }
            if (++checked_sectors==3) return;
    case 7: if (pic[pic_ci+1][pic_cj+1] && !win[win_ci+1][win_cj+1]) {
              win[win_ci+1][win_cj+1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci+1,win_cj+1,pic,pic_ci+1,pic_cj+1,
                                      min_i,min_j,max_i,max_j,!sector?6:sector);
            }
            if (++checked_sectors==3) return;
    case 8: if (pic[pic_ci+1][pic_cj] && !win[win_ci+1][win_cj]) {
              win[win_ci+1][win_cj]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci+1,win_cj,pic,pic_ci+1,pic_cj,
                                      min_i,min_j,max_i,max_j,!sector?7:sector);
            }
            if (++checked_sectors==3) return;
            if (pic[pic_ci+1][pic_cj-1] && !win[win_ci+1][win_cj-1]) {
              win[win_ci+1][win_cj-1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci+1,win_cj-1,pic,pic_ci+1,pic_cj-1,
                                      min_i,min_j,max_i,max_j,!sector?8:sector);
            }
            if (++checked_sectors==3) return;
            if (pic[pic_ci][pic_cj-1] && !win[win_ci][win_cj-1]) {
              win[win_ci][win_cj-1]=(pic_type)OBJECT_PIX_VAL;
              get_connected_neighbors(win,win_ci,win_cj-1,pic,pic_ci,pic_cj-1,
                                      min_i,min_j,max_i,max_j,!sector?1:sector);
            }
    default:return;
  }
}
