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
 * File:    rht_imgutils.c
 * Purpose: image conversion
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

/*****************************************************************************/
convert_pic_vec(pic, dim_x, dim_y, vec, n)
/*

  Convert image to edge point vector.

  Parameters :  pic - image
                dim_x,dim_y - image quantization
                vec - vector
                n - vector length

*/
pic_type pic[][MAX_SIZE];
int dim_x, dim_y;
pic_vec_type vec[];
long *n;
{
  int i,j;
  long k=0L;

  for (i=0; i<dim_y; i++)
    for (j=0; j<dim_x; j++)
      if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
        vec[k].i=i;
        vec[k].j=j;
        k++;
      }

  *n=k;
}
