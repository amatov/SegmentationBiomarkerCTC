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
 * File:    rwrht_imgutils.c
 * Purpose: image conversion
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"
#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif

/*****************************************************************************/
convert_win_pic_vec(pic, dim_x, dim_y, cx, cy, size, vec, n)
/*

  Convert window to edge point vector.

  Parameters :  pic - image
                dim_x,dim_y - image quantization
                cx,cy - window center point
                size - window size ((2*size+1)^2)
                vec - vector
                n - vector length

*/
pic_type pic[][MAX_SIZE];
int dim_x, dim_y, cx, cy, size;
pic_vec_type vec[];
long *n;
{
  int i,j,min_i=cy-size,min_j=cx-size,
          max_i=cy+size+1,max_j=cx+size+1;
  long k=0L;

  if (min_i<0) min_i=0;
  if (min_j<0) min_j=0;
  if (max_i>dim_y) max_i=dim_y;
  if (max_j>dim_x) max_j=dim_x;

  for(i=min_i; i<max_i; i++)
    for(j=min_j; j<max_j; j++) {
      PutPixel_gray_MACRO(quads[1],j,i,MID_GRAY_LEVEL);
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
        vec[k].i=i;
        vec[k].j=j;
        k++;
        PutPixel_MACRO(quads[1],j,i);
      }
    }

  *n=k;
}
