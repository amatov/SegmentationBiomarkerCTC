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
 * File:    ht_nonansi.c
 * Purpose: non-ansi features
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

#ifdef NON_ANSI_RECOVERY
/*****************************************************************************/
int nint(x)
/*

  Round the input value to the nearest integer.

  Parameters :  x - input floating point value

*/
double x;
{
  if ((fabs(x)>=MAXINT) || isinf(x) || isnan(x))
    if (x>0.0)
      return MAXINT;
    else
      return -(MAXINT);
  return (int)(x+0.5);
}

/*****************************************************************************/
double rint(x)
/*

  Round the input value to the nearest integer but preserve value type.

  Parameters :  x - input floating point value

*/
double x;
{
  return (double)floor(x+0.5);
}
#endif
