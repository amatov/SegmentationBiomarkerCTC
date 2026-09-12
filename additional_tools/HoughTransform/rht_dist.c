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
 * File:    rht_dist.c
 * Purpose: point distance calculation
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

/*****************************************************************************/
maxdist(x1,y1,x2,y2)
/*

  Calculate x-axis and y-axis distances and reuturn longer.

  Parameters :  x1,y1,x2,y2 - points coordinates

*/
int x1,y1,x2,y2;
{
  int x,y;

  x=abs(x1-x2);
  y=abs(y1-y2);

  if (x<y)
    return y;

  return x;
}

/*****************************************************************************/
mindist(x1,y1,x2,y2)
/*

  Calculate x-axis and y-axis distances and reuturn sum.

  Parameters :  x1,y1,x2,y2 - points coordinates

*/
int x1,y1,x2,y2;
{
  return abs(x1-x2)+abs(y1-y2);
}

/*****************************************************************************/
euclid_dist(x1,y1,x2,y2)
/*

  Calculate euclide distance and reuturn it.

  Parameters :  x1,y1,x2,y2 - points coordinates

*/
int x1,y1,x2,y2;
{
  int x,y;

  x=x1-x2;
  y=y1-y2;

  return nint(sqrt((double)(x*x+y*y)));
}
