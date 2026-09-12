/*****************************************************************************
 * Houghtool - The Software Package for efficiency measuring and visualization
 * of the HT and it's variants for line detection
 *
 * Lappeenranta University of Technology, Department of Information Technology
 * Laboratory of Information Processing
 *
 * Authors: Petri Hirvonen, Jouni Ikonen (Jouni.Ikonen@lut.fi)
 *	    Pekka Kultanen, Heikki K{lvi{inen (Heikki.Kalviainen@lut.fi)
 *
 * File:    cfht_fndmaxline.c
 * Purpose: verifying lines found
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"

extern int StopFlag;

/*****************************************************************************/
ReadAccuAndMarkLineSegments(start,Mat,Params,Threshold,NumOfMaxs,pic,
                            dim_x,dim_y,MinSegLen,MaxSegWidth,MaxSegGap,
			    false_alarms)
/*

  Follow line defined by Hough accumulator space maxima location and
  form line segments according to given minimum segment length and 
  given deviation from absolute line.

  Parameters :	start - data structure for line segments
		Mat - dynamic accumulator
		Params - real line parameter list
		NumOfMaxs - number of lines to find
		pic - binary image containing  edge points, 
		       pic_type matrix size [][MAX_SIZE]
		dim_x, dim_y - edge image size
		MinSegLen - minimum legth for line segment
		MaxSegWidth - line scanning width
		MaxSegGap - maximum gap between pixels within line segment
		false_alarms - number of false alarms

*/
line **start;
InfMat *Mat;
double Params[][2];
int Threshold,NumOfMaxs,dim_x,dim_y,MinSegLen,MaxSegWidth,MaxSegGap,
    *false_alarms;
pic_type pic[][MAX_SIZE];
{
  InfIndex *p;
  Accu *s, *s_max;
  int n, i=0, thr_max=0, val=MAXINT, k, l;

  do {
    p=Mat;
    n=0;

    if (p->Accu) { /* infinite slope */
      s=p->Accu;
      do {
        if (s->Accumulator>=Threshold && s->Accumulator>thr_max) {
          thr_max=s->Accumulator;
          Params[i][0]=INF_DOUBLE_VAL;
          Params[i][1]=s->Y_coord;
          s_max=s;
          n++;
        }
        s=s->Up;
      } while (s);
    }

    if (p->NextIndex) { /* finite slope */
      p=p->NextIndex;

      do {
        s=p->Accu;
        do {
          if (s->Accumulator>=Threshold && s->Accumulator>thr_max) {
            thr_max=s->Accumulator;
            Params[i][0]=p->X_coord;
            Params[i][1]=s->Y_coord;
            s_max=s;
            n++;
          }
          s=s->Up;
        } while (s);
        p=p->NextIndex;
      } while (p);
    }

    if (n) { /* at least one cell was found */
/*
      printf("a: %lf b: %lf val: %d i: %d n: %d\n",
             Params[i][0],Params[i][1],s_max->Accumulator,i,n);
*/
      FollowAndMarkLineSegments_a_b(start,Params[i][0],Params[i][1],
                                    pic,dim_x,dim_y,MinSegLen,MaxSegWidth,
				    MaxSegGap);

      if (RemoveEdgePoints(*start,pic,dim_x,dim_y)) {
        i++;
        for (k=val=0; k<dim_y; k++)
          for (l=0; l<dim_x; l++)
            if (pic[k][l]==(pic_type)OBJECT_PIX_VAL)
              val++;
      } else /* there is no matching edges -> false alarm */
	(*false_alarms)++;

      s_max->Accumulator=thr_max=0;

    }

  } while(n && i<NumOfMaxs && val>MinSegLen && !StopFlag);

  return i;
}
