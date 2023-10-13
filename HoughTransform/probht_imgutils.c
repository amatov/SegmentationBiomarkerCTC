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
 * File:    probht_imgutils.c
 * Purpose: sampling image
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

long rnd();

/*****************************************************************************/
take_a_sample_of_pic(pic, pic1, dim_x, dim_y, level)
/*

  Copy an amount of points of pic1 defined by level to pic.

  Parameters :  pic,pic1 - images
                dim_x,dim_y - image quantization
                level - sample level (percentage)

*/
pic_type pic[][MAX_SIZE], pic1[][MAX_SIZE];
int dim_x, dim_y;
double level;
{
  int i, j, k, n=0;

  for (i=0; i<dim_x; i++)
    for (j=0; j<dim_y; j++) {
      pic[i][j]=(pic_type)0;
      if (pic1[i][j]==(pic_type)OBJECT_PIX_VAL)
        n++; /* count of edges */
    }

  n=nint((double)n*(level/100.0)); /* new number of edges */

  for (k=0; k<n; ) {
    i=(int)rnd(dim_x-1); /* random coordinates */
    j=(int)rnd(dim_y-1);
    if (pic1[i][j]==(pic_type)OBJECT_PIX_VAL &&
        pic[i][j]!=(pic_type)OBJECT_PIX_VAL) {
      pic[i][j]=(pic_type)OBJECT_PIX_VAL;
      k++;
    }
  }
}
