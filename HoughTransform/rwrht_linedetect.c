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
 * File:    	rwrht_linedetect.c
 * Purpose: 	Random Window Randomized Hough Transform (RWRHT) for lines
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

pic_vec_type vec[MAX_VEC_COORD], winvec[MAX_VEC_COORD];

long rnd();

/*****************************************************************************/
rwrht_line(pic,pic1,gpic,rpic,dim_x,dim_y,real_params,NumOfMaxs,MinSegLen,
           MaxSegWidth,MaxSegGap,MinDist,MaxDist,Accuracy,Threshold,MinWinSize,
           MaxWinSize,AccumulationsLimit,NumOfTests,TextInfo,PostScript,
	   ParametersOut)
/*

  Implement the Random Window RHT (RWRHT) for line detection.

  Parameters :  pic,pic1,gpic,rpic - input and output images
                dim_x,dim_y - image quantization
                real_params - real line parameter list
                NumOfMaxs - number of lines to find
                MinSegLen,MaxSegWidth,MaxSegGap - line definitions
                MinDist,MaxDist - point distance limits
                Accuracy - accumulator resolution
                Threshold - minimum score accepted for an accumulator maximum
                MinWinSize,MaxWinSize - window size limits
                AccumulationsLimit - are accumulations from the same window
                                     limited (on/off)
                NumOfTests - number of tests
                TextInfo - statistics output selector
                PostScript - PostScript output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE],pic1[][MAX_SIZE],gpic[][MAX_SIZE],rpic[][MAX_SIZE];
int dim_x,dim_y,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,MinDist,MaxDist,
    Threshold,MinWinSize,MaxWinSize,AccumulationsLimit,NumOfTests,TextInfo,
    PostScript;
double real_params[][2],Accuracy;
char ParametersOut[256];
{
  InfMat *acc_space;
  Accu *Cell;
  line *start=NULL;
  int test,cx,cy,x1,x2,y1,y2,maxval,i,k,NumOfAcc=0,max_size,WinSize,
      NumOfSel,Selections,LineFound,false_alarms=0,found_params,
      NumOfAccumulations,Accumulations,WinSelections;
  long adv1, adv2, timecount=0L, n, win_n, cumulative_size,
       num_of_real_params=0L, acquired_real_params=0L;
  float totcputime=0.0;
  double a,b,amin,amax,bmin,bmax,round=1./Accuracy;
  struct tms cputimebase,cputime;

  if (AccumulationsLimit)
    NumOfAccumulations=20*Threshold;
  else
    NumOfAccumulations=MAXINT;

  StopFlag=0;

  NotifyDispatch_MACRO();

  CreateInfMat(&acc_space);

  init_random_generator();

  times(&cputimebase);

  for (test=0; test<NumOfTests && !StopFlag; test++) {

    copy_pic(pic,pic1,dim_x,dim_y);
    convert_pic_vec(pic,dim_x,dim_y,vec,&n);
    timecount=0L; NumOfAcc=0;
    if (test) DeleteAllLineSegments(start);

    for (k=0; (k<NumOfMaxs) && (n>(2*MinSegLen)) && (NumOfAcc<=10) &&
              !StopFlag; ) {

      LineFound=FALSE;
      WinSelections=Selections=Accumulations=0;

      /* choose window center point */

      do {

        adv1=rnd(n-1);
        cx=vec[adv1].j;
        cy=vec[adv1].i;

        /* choose window size */

        WinSize=(int)rnd(MaxWinSize-MinWinSize)+MinWinSize;

        /* create pixel vector using windowed edge points */

        convert_win_pic_vec(pic,dim_x,dim_y,cx,cy,WinSize,winvec,&win_n);

        WinSelections++;

        NotifyDispatch_MACRO();

      } while (/*win_n<=WinSize*/win_n<MinSegLen && /* select new if too empty */
               WinSelections<n && !StopFlag);

      if (win_n<MinSegLen && WinSelections>=n) {
        if (Accumulations) {
          RemoveAccumulatorSpace(acc_space);
          CreateInfMat(&acc_space);
        }
        break; /* only noise points left (?) in the image */
      }

      /* choose the number of maximum RHT selections */

      /* NumOfSel = f(win_n) */
      NumOfSel=win_n*win_n;

      /* do the RHT for the selected window */

      do { /* pick and accumulate */
        if (AccumulationsLimit)
          Selections=0;

        /* choose point pair randomly from binary picture (window) */

        adv1=rnd(win_n-1);

        do { /* pick */
          x1=winvec[adv1].j;
          y1=winvec[adv1].i;
	       
          while (adv1==(adv2=rnd(win_n-1)));

          x2=winvec[adv2].j;
          y2=winvec[adv2].i; 
          adv1=adv2;
          Selections++;

          NotifyDispatch_MACRO();

        } while ((maxdist(x1,y1,x2,y2)>MaxDist || mindist(x1,y1,x2,y2)<MinDist)
                 && Selections<NumOfSel && !StopFlag);

        if (Selections>=NumOfSel && (maxdist(x1,y1,x2,y2)>MaxDist ||
                                     mindist(x1,y1,x2,y2)<MinDist) ) {
          if (Accumulations) {
            RemoveAccumulatorSpace(acc_space);
            CreateInfMat(&acc_space);
          }
          NumOfAcc++;
          break; /* no suitable points in the window */
        }

        /* calculate parameters & accumulate parameter space */

        a = (double)(y2-y1)/(double)(x2-x1);
        b = (double)y1-a*x1;

        /* Infinity value handling */
        if (!isinf(a))
          Cell=IncAccu(rint(round*a)/round,rint(round*b)/round,acc_space,1);
        else
          Cell=IncAccu(a,(double)x1,acc_space,1);

        Accumulations++;
        timecount++;

        /* show points which have been selected until now */

        PutPixel_MACRO(quads[1], x1, y1);
        PutPixel_MACRO(quads[1], x2, y2);
        FlushDisp_MACRO();
        NotifyDispatch_MACRO();

        size[timecount]=ExamineAccu(acc_space,&amin,&amax,&bmin,&bmax,&maxval);

        PutLongText_MACRO(quads[1],0,250,timecount,!(timecount % 100));

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

          if (isinf(a))
            b=(double)x1;
          else {
            a=rint(round*a)/round;
            b=rint(round*b)/round;
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

          RemoveAccumulatorSpace(acc_space);
          CreateInfMat(&acc_space);
          NumOfAcc++;

          ClearQuadrant_MACRO(quads[1]);
          PutLongText_MACRO(quads[1],0,250,timecount,TRUE);
        }

      } while (!LineFound && Selections<NumOfSel &&
               Accumulations<NumOfAccumulations && !StopFlag);
      NotifyDispatch_MACRO();
    } /* endfor(maxs) */
    maxs_histogram[k]++;
    if ((int)real_params[0][0]>0) {
      found_params=test_found_params(maxs,k,real_params);
      num_of_real_params+=found_params;
      if (found_params==NumOfMaxs)
        acquired_real_params++;
    }
  } /* endfor(tests) */

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
