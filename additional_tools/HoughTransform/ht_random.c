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
 * File:    	ht_random.c
 * Purpose: 	random value generator
 * Date:    	Jun 1, 1993
 * Last change: Nov 30,1995
 *****************************************************************************/

#include "ht_hough.h"

#include <sys/timeb.h>

/*****************************************************************************/
long rnd(dev)
/*

  Give a random long integer value between 0 and dev.

  Parameters :  dev - max value

*/
long dev;
{
  double drand48();

  return((long)floor((dev+1)*drand48()));
}

/*****************************************************************************/
init_random_generator()
/*

  Initialize the random generator.

  Parameters :

*/
{
  struct timeb seedvalue;

  if (time(0) == -1) {
     ftime(&seedvalue);
     srand48(seedvalue.time);
     
  }
  else
     srand48(time(0)); 
}
