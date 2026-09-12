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
 * File:    cfht_infmat.c
 * Purpose: functions for manipulating the accumulator of the CFHT
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"

#define sqr(a) ((a)*(a))

int timecount;

/*****************************************************************************/
Accu *IncAccu_cfht(X,Y,Mat,Epsilon,Gamma,val)
/*

  Increment an accumulator cell by 'val'.

  Parameters :	X,Y - parameter values
		Mat - dynamic accumulator space
		range - Hough accumulator space range selector
		Epsilon - Euclidean distance between cells when merging them 
		Gamma - score weigth
		val - increment

*/
double X,Y,Epsilon,Gamma;
InfMat *Mat;
int val;
{
  InfIndex *p=Mat,*prev_p=Mat,*near_p,*near_prev_p;
  Accu *s,*New,*prev_s,*near_prev_s,*near_s;
  double new_X,new_Y,near_X,near_Y,gap,near_gap=MAXDOUBLE;
  int new_Accuval,found=0;

  if (isinf(X)) { /* infinite slope */

    if (p->Accu==(Accu *)NULL) {

      New=AllocAccu();
      New->Up=(Accu *)NULL;
      p->Accu=New;
      New->Y_coord=Y;
      New->Accumulator=val; 
      New->Time=timecount;

      return New;

    } else {

      s=prev_s=p->Accu;

      while (s->Up && s->Y_coord<=(Y-Epsilon)) {
        if (s!=prev_s)
          prev_s=s;
        s=s->Up;
      }

      while (s && (gap=fabs(s->Y_coord-Y))<Epsilon) {
        if (gap<near_gap) {
          near_gap=gap;
          near_s=s;
          near_prev_s=prev_s;
          near_Y=s->Y_coord;
          found=1;
        }

        if (s!=prev_s)
          prev_s=s;
        s=s->Up;
      }

      if (found) {
        new_Accuval=near_s->Accumulator+val;
        if (Gamma)
          new_Y=(1.0-Gamma)*near_s->Y_coord+Gamma*Y;
        else
          new_Y=(double)(near_s->Accumulator*near_s->Y_coord+val*Y)/
                (double)(new_Accuval);
        if (near_s==p->Accu) {
          p->Accu=near_s->Up;
          free((Accu *)near_s);
        } else
          DeleteAccu(near_prev_s);
        return IncAccu(X,new_Y,Mat,new_Accuval); /* incrementing one cell */
      } else
	return IncAccu(X,Y,Mat,val); /* new cell */

    }

  } else { /* finite slope */

    if (p->NextIndex==(InfIndex *)NULL) {

      p=AddInfIndex(p, X);
  
      New=AllocAccu();
      New->Up=(Accu *)NULL;
      p->Accu=New;
      New->Y_coord=Y;
      New->Accumulator=val;
      New->Time=timecount; 
 
      return New;

    } else {

      p=p->NextIndex;

      while (p->NextIndex && p->X_coord<=(X-Epsilon)) {
        if (p!=prev_p)
          prev_p=p;
        p=p->NextIndex;
      }

      while (p && p->X_coord<(X+Epsilon)) {

        s=prev_s=p->Accu;

        while (s->Up && s->Y_coord<=(Y-Epsilon)) {
          if (s!=prev_s)
            prev_s=s;
          s=s->Up;
        }

        while (s && s->Y_coord<(Y+Epsilon)) {
          if (((gap=sqrt(sqr(p->X_coord-X)+sqr(s->Y_coord-Y)))<Epsilon) &&
              (gap<near_gap)) {
            near_gap=gap;
            near_p=p;
            near_prev_p=prev_p;
            near_s=s;
            near_prev_s=prev_s;
            near_X=p->X_coord;
            near_Y=s->Y_coord;
            found=1;
          }

          if (s!=prev_s)
            prev_s=s;
          s=s->Up;
        }

        if (p!=prev_p)
          prev_p=p;
        p=p->NextIndex;
      }

      if (found) {
        new_Accuval=near_s->Accumulator+val;
        if (Gamma) {
          new_X=(1.0-Gamma)*near_p->X_coord+Gamma*X;
          new_Y=(1.0-Gamma)*near_s->Y_coord+Gamma*Y;
        } else {
          new_X=(double)(near_s->Accumulator*near_p->X_coord+val*X)/
                (double)(new_Accuval);
          new_Y=(double)(near_s->Accumulator*near_s->Y_coord+val*Y)/
                (double)(new_Accuval);
        }
        if (near_s==near_p->Accu) {
          if (near_s->Up==(Accu *)NULL)
            DeleteInfIndex(near_prev_p);
          else {
            near_p->Accu=near_s->Up;
            free((Accu *)near_s);
          }
        } else
          DeleteAccu(near_prev_s);
        return IncAccu(new_X,new_Y,Mat,new_Accuval); /* increment one cell */
      } else
	return IncAccu(X,Y,Mat,val); /* new cell */
    }
  }
}
