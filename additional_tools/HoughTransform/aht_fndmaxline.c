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
 * File:    aht_fndmaxline.c
 * Purpose: verifying lines found
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

/*****************************************************************************/
void FollowAndMarkLineSegments_m_c(start,m,c,range,edge,dim_x,dim_y,min_seg_length,
                              max_seg_width,max_seg_gap)
/*

  Follow line defined by Hough accumulator space maxima location and
  form line segments according to given minimum segment length and 
  given deviation from absolute line.

  Parameters :	start - data structure for line segments
		m, c - Hough accumulator space coordinates
		range - Hough accumulator space range selector
		edge - binary image containing  edge points, 
		       pic_type matrix size [][MAX_SIZE]
		dim_x, dim_y - edge image size
		min_seg_length - minimum legth for line segment
		max_seg_width  - maximum width for line segment
		max_seg_gap    - maximum gap between pixels within
                                 line segment

  Remarks : ...

*/

line **start;
double *m, *c;
int range, dim_x, dim_y, min_seg_length, max_seg_width, max_seg_gap;
pic_type edge[][MAX_SIZE];
{
  double a=*m, b=*c;

  if (range!=1) {               /* (m, c)  --->  (a, b) */
    if (fabs(*m)<0.00001)
      a=((*m)/fabs(*m))*INF_DOUBLE_VAL;
    else {
      a=1.0/(*m);
      b=-(*c)/(*m);
    }
    *m=a;
    *c=b;
  }
  /* do the actual job elsewhere */
  FollowAndMarkLineSegments_a_b(start,a,b,edge,dim_x,dim_y,min_seg_length,
                                max_seg_width,max_seg_gap);
}
