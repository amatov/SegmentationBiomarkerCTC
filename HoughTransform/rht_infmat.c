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
 * File:    rht_infmat.c
 * Purpose: manipulating the accumulator of the RHT
 * Date:    Jun 1, 1993
 *****************************************************************************/

/*****************************************************************************/
/*                                Accumulator                                */
/*                                                                           */
/*  Accu       Accu                    Accu        Accu                      */
/*    ^          ^                       ^           ^                       */
/*    |          |                       |           |                       */
/*  Accu       Accu        Accu        Accu        Accu        Accu          */
/*    ^          ^           ^           ^           ^           ^           */
/*    |          |           |           |           |           |           */
/* InfMat -> infIndex -> infIndex -> infIndex -> infIndex -> infIndex        */
/*                                                                           */
/*  (inf)  |                        (finite)                                 */
/*****************************************************************************/

#include "ht_hough.h"
#include "rht_infmat.h"

int timecount;

/*****************************************************************************/
InfMat *AllocInfMat() /*** alloc space for accumulator root ***/
{
  InfMat *p=(InfMat *)malloc(sizeof(InfMat));
  if (p==(InfMat *)NULL) {
    err_prnt("Cannot allocate InfMat.");
    exit(-1);
  }
  p->NextIndex=(InfIndex *)NULL;
  p->Accu=(Accu *)NULL;
  p->X_coord=MAXDOUBLE;

  return p;
}

/*****************************************************************************/
InfIndex *AllocInfIndex() /*** alloc space for `first' index ***/
{
  InfIndex *p=(InfIndex *)malloc(sizeof(InfIndex));
  if (p==(InfIndex *)NULL) {
    err_prnt("Cannot allocate InfIndex.");
    exit(-1);
  }

  return p;
}

/*****************************************************************************/
Accu *AllocAccu() /*** alloc space for accumulator cell ***/
{
  Accu *p=(Accu *)malloc(sizeof(Accu));
  if (p==(Accu *)NULL) {
    err_prnt("Cannot allocate Accu.");
    exit(-1);
  }

  return p;
}

/*****************************************************************************/
CreateInfMat(Mat) /*** create accumulator ***/
InfMat **Mat;
{
  *Mat=AllocInfMat();
}

/*****************************************************************************/
InfIndex *AddInfIndex(Previous, Value) /*** add new `first' index ***/
InfIndex *Previous;
double Value;
{
  InfIndex *New;

  New=AllocInfIndex();

  New->NextIndex=Previous->NextIndex;
  Previous->NextIndex=New;
  New->Accu=(Accu *)NULL;
  New->X_coord=Value;

  return New;
}

/*****************************************************************************/
DeleteInfIndex(Previous) /*** delete `first' index ***/
InfIndex *Previous;
{
  InfIndex *Item=Previous->NextIndex;
  Previous->NextIndex=Item->NextIndex;

  free((char *)(Item->Accu));
  free((char *)Item);
}

/*****************************************************************************/
Accu *AddAccu(Previous,Coord,val) /*** add new cell ***/
Accu *Previous;
double Coord;
int val;
{
  Accu *New;

  New=AllocAccu();
  New->Up=Previous->Up;
  Previous->Up=New;
  New->Y_coord=Coord;
  New->Accumulator=val;
  New->Time=timecount; 

  return New;
}

/*****************************************************************************/
DeleteAccu(Previous) /*** delete cell ***/
Accu *Previous;
{
  Accu *Item=Previous->Up;
  Previous->Up=Item->Up;
  free((char *)Item);
}

/*****************************************************************************/
Accu *IncAccu(X,Y,Mat,val) /*** increment cell score ***/
double X,Y;
InfMat *Mat;
int val;
{
  InfIndex *p=Mat;
  Accu *s,*New;

  if (isinf(X)) { /* infinite `first' index */
    if (p->Accu==(Accu *)NULL) {
      New=AllocAccu();
      New->Up=(Accu *)NULL;
      p->Accu=New;
      New->Y_coord=Y; /* Notice : this space is not for Y-coordinate. */
      New->Accumulator=val; 
      New->Time=timecount;
 
      return New; 
    }

    s=p->Accu;
    do {
      if (s->Y_coord==Y) {
        s->Accumulator+=val;
        s->Time=timecount;
        return s;
      } else
        if ((s->Up==(Accu *)NULL) || (s->Up->Y_coord>Y))
          return AddAccu(s,Y,val);
      s=s->Up;
    } while (1);

  }

  if (p->NextIndex==(InfIndex *)NULL) {
    p=AddInfIndex(p, X);
  
    New=AllocAccu();
    New->Up=(Accu *)NULL;
    p->Accu=New;
    New->Y_coord=Y;
    New->Accumulator=val;
    New->Time=timecount; 
 
    return New; 
  }

  p=p->NextIndex;

  do {
    if (p->X_coord==X) { 
      s=p->Accu;
      do {
        if (s->Y_coord==Y) {
          s->Accumulator+=val;
          s->Time=timecount;
          return s;
        } else
          if ((s->Up==(Accu *)NULL) || (s->Up->Y_coord>Y))
            return AddAccu(s,Y,val);
        s=s->Up;
      } while (1);
    } else
      if ((p->NextIndex==(InfIndex *)NULL) || (p->NextIndex->X_coord>X)) { 
        p=AddInfIndex(p, X);
  
        New=AllocAccu();
        New->Up=(Accu *)NULL;
        p->Accu=New;
        New->Y_coord=Y;
        New->Accumulator=val; 
        New->Time=timecount;
 
        return New; 
      }
      p=p->NextIndex;
      
  } while (1);
}

/*****************************************************************************/
ExamineAccu(Mat,XMin,XMax,YMin,YMax,AccuMax) /*** count cells ***/
InfMat *Mat;
double *XMin,*XMax,*YMin,*YMax;
int *AccuMax;
{
  InfIndex *p=Mat;
  Accu *s;
  int n=0;

  if ((p->NextIndex==(InfIndex *)NULL) && (p->Accu==(Accu *)NULL))
    return 0;

  if (p->Accu) {
    s=p->Accu;
    if (s->Y_coord<*YMin)
      *YMin=s->Y_coord;
    do {
      if (s->Accumulator>*AccuMax) *AccuMax=s->Accumulator;
      if (s->Y_coord>*YMax) *YMax=s->Y_coord;

      s=s->Up; n++;
    } while (s);

  }

  if (p->NextIndex==(InfIndex *)NULL)
    return n;

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X_coord;
  *XMax=p->X_coord;

  s=p->Accu;
  *YMin=s->Y_coord;
  *YMax=s->Y_coord;

  do {
    s=p->Accu;
    if (s->Y_coord<*YMin)
      *YMin=s->Y_coord;

    do {
      if (s->Accumulator>*AccuMax)
        *AccuMax=s->Accumulator;
      if (s->Y_coord>*YMax)
        *YMax=s->Y_coord;

      s=s->Up;
      n++;
    } while (s);

    if (*XMin>p->X_coord)
      *XMin=p->X_coord;
    if (*XMax<p->X_coord)
      *XMax=p->X_coord;

    p=p->NextIndex;
  } while (p);

  return n;
}

/*****************************************************************************/
ExamineAccu2(Mat,XMin,XMax,YMin,YMax,AccuMax,XaMax,YaMax) /*** count cells ***/
InfMat *Mat;
double *XMin,*XMax,*YMin,*YMax,*XaMax,*YaMax;
int *AccuMax;
{
  InfIndex *p=Mat;
  Accu *s;
  int n=0;

  if (p->NextIndex==(InfIndex *)NULL)
    return 0;

  /* Infinity value checking here ??? */

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X_coord;
  *XMax=p->X_coord;

  s=p->Accu;
  *YMin=s->Y_coord;
  *YMax=s->Y_coord;

  do {
    s=p->Accu;
    if (s->Y_coord<*YMin)
      *YMin=s->Y_coord;

    do {
      if (s->Accumulator>*AccuMax) {
        *AccuMax=s->Accumulator;
        *YaMax=s->Y_coord;
        *XaMax=p->X_coord;
      }
      if (s->Y_coord>*YMax)
        *YMax=s->Y_coord;

      s=s->Up;
      n++;
    } while (s);

    if (*XMin>p->X_coord)
      *XMin=p->X_coord;
    if (*XMax<p->X_coord)
      *XMax=p->X_coord;

    p=p->NextIndex;
  } while (p);

  return n;
}


/*****************************************************************************/
PrintAccu(Mat) /*** print accumulator contents ***/
InfMat *Mat;
{
  InfIndex *p=Mat;
  Accu *s;

  /* Infinity value printing */
  if (p->Accu!=(Accu *)NULL) {
    s=p->Accu;
    do {
      printf("(x,y) = %lf, %lf val = %d\n ",
             p->X_coord,s->Y_coord,s->Accumulator);
	    
      s=s->Up;
    } while (s);
  }

  if (p->NextIndex==(InfIndex *)NULL)
    return;

  p=p->NextIndex;

  s=p->Accu; 

  do {
    s=p->Accu;
    do {
      printf("(x,y) = %lf, %lf val = %d\n ",p->X_coord,s->Y_coord,s->Accumulator);
	    
      s=s->Up;;
    } while (s);

    p=p->NextIndex;
  } while (p);

}

/*****************************************************************************/
PrintAccu2(Mat,nhits) /*** print accumulator contents ***/
InfMat *Mat;
int nhits;
{
  InfIndex *p=Mat;
  Accu *s;

  if (p->NextIndex==(InfIndex *)NULL)
    return 0;

  /* Infinity value printing */
  if (p->Accu!=(Accu *)NULL) {
    s=p->Accu;
    do {
      if (nhits == s->Accumulator)
        printf("(%ld,%ld) ",(int)rint(p->X_coord),(int)rint(s->Y_coord));
      s=s->Up;;
    } while (s);
  }

  p=p->NextIndex;

  s=p->Accu; 

  do {
    s=p->Accu;
    do {
      if (nhits == s->Accumulator)
            printf("(%ld,%ld) ",(int)rint(p->X_coord),(int)rint(s->Y_coord));

      s=s->Up;;
    } while (s);

    p=p->NextIndex;
  } while (p);

}

/*****************************************************************************/
RemoveAccumulatorSpace(Mat) /*** remove the whole accumulator ***/
InfMat *Mat;
{
  InfIndex *p=Mat, *tmp_p;
  Accu *s, *tmp_s;

  if (p->NextIndex==(InfIndex *)NULL)
    return 0;

  if (p->Accu!=(Accu *)NULL) {
    s=p->Accu;
    do {
      tmp_s=s->Up;
      free(s);
      s=tmp_s;
    } while (s);
  }

  free(p);
  p=p->NextIndex;

  s=p->Accu; 

  do {
    s=p->Accu;
    do {
      tmp_s=s->Up;
      free(s);
      s=tmp_s;
    } while (s);
    tmp_p=p->NextIndex;
    free(p);
    p=tmp_p;
  } while (p);
}


/*****************************************************************************/
InfMat3D *AllocInfMat3D() /*** alloc space for 3D accumulator root ***/
{
  InfMat3D *p=(InfMat3D *)malloc(sizeof(InfMat3D));
  if (p==(InfMat3D *)NULL) {
    err_prnt("Cannot allocate InfMat3D.");
    exit(-1);
  }
  p->NextIndex=(InfIndex3DX *)NULL;
  p->YIndex=(InfIndex3DY *)NULL;
  p->X=MAXDOUBLE;

  return p;
}

/*****************************************************************************/
InfIndex3DX *AllocInfIndex3DX() /*** alloc space for first index ***/
{
  InfIndex3DX *p=(InfIndex3DX *)malloc(sizeof(InfIndex3DX));
  if (p==(InfIndex3DX *)NULL) {
    err_prnt("Cannot allocate InfIndex3DX.");
    exit(-1);
  }

  return p;
}

/*****************************************************************************/
InfIndex3DY *AllocInfIndex3DY() /*** alloc space for second index ***/
{
  InfIndex3DY *p=(InfIndex3DY *)malloc(sizeof(InfIndex3DY));
  if (p==(InfIndex3DY *)NULL) {
    err_prnt("Cannot allocate InfIndex3DY.");
    exit(-1);
  }

  return p;
}

/*****************************************************************************/
InfIndex3DZ *AllocInfIndex3DZ() /*** alloc space for third index (cell) ***/
{
  InfIndex3DZ *p=(InfIndex3DZ *)malloc(sizeof(InfIndex3DZ));
  if (p==(InfIndex3DZ *)NULL) {
    err_prnt("Cannot allocate InfIndex3DZ.");
    exit(-1);
  }

  return p;
}

/*****************************************************************************/
CreateInfMat3D(Mat) /*** create 3D accumulator ***/
InfMat3D **Mat;
{
  *Mat=AllocInfMat3D();
}

/*****************************************************************************/
InfIndex3DX *AddInfIndex3DX(Previous, Value) /*** add new `first' index ***/
InfIndex3DX *Previous;
double Value;
{
  InfIndex3DX *New;

  New=AllocInfIndex3DX();

  New->NextIndex=Previous->NextIndex;
  Previous->NextIndex=New;
  New->YIndex=(InfIndex3DY *)NULL;
  New->X=Value;

  return New;
}

/*****************************************************************************/
InfIndex3DY *AddInfIndex3DY(Previous, Value) /*** add new `second' index ***/
InfIndex3DY *Previous;
double Value;
{
  InfIndex3DY *New;

  New=AllocInfIndex3DY();

  New->NextIndex=Previous->NextIndex;
  Previous->NextIndex=New;
  New->ZIndex=(InfIndex3DZ *)NULL;
  New->Y=Value;

  return New;
}

/*****************************************************************************/
InfIndex3DZ *AddInfIndex3DZ(Previous, Coord, Value) /*** add new cell ***/
InfIndex3DZ *Previous;
double Coord, Value;
{
  InfIndex3DZ *New;

  New=AllocInfIndex3DZ();

  New->NextIndex=Previous->NextIndex;
  Previous->NextIndex=New;
  New->Z=Coord;
  New->Accumulator=Value;

  return New;
}

/*****************************************************************************/
InfIndex3DZ *AddInfIndex3DZ2(Previous, Coord, Value, XVal1,YVal1,XVal,YVal)
InfIndex3DZ *Previous;                                   /*** add new cell ***/
double Coord,XVal1,YVal1,XVal,YVal; 
int Value;
{
  InfIndex3DZ *New;

  New=AllocInfIndex3DZ();

  New->NextIndex=Previous->NextIndex;
  Previous->NextIndex=New;
  New->Z=Coord;
  New->Accumulator=Value; New->Hit=1;
  New->XValue=XVal; New->YValue=YVal;
  New->XValue1=XVal1; New->YValue1=YVal1;

  return New;
}

/*****************************************************************************/
InfIndex3DZ *IncAccu3D(X,Y,Z,Mat,val) /*** increment cell score ***/
double X,Y,Z;
InfMat3D *Mat;
int val;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s,*NewY;
  InfIndex3DZ *r,*NewZ;

  if (p->NextIndex==(InfIndex3DX *)NULL) {
    p=AddInfIndex3DX(p, X);
  
    NewY=AllocInfIndex3DY();
    NewY->NextIndex=(InfIndex3DY *)NULL;
    p->YIndex=NewY;
    NewY->Y=Y;

    NewZ=AllocInfIndex3DZ();
    NewZ->NextIndex=(InfIndex3DZ *)NULL;
    NewY->ZIndex=NewZ;
    NewZ->Z=Z;
    NewZ->Accumulator=val;

    return NewZ; 
  }

  p=p->NextIndex;

  do {
    if (p->X==X) { 
      s=p->YIndex;
      do {
        if (s->Y==Y) {  
          if (s->ZIndex==(InfIndex3DZ *)NULL) { 
            NewZ=AllocInfIndex3DZ();
            NewZ->NextIndex=(InfIndex3DZ *)NULL;
            NewZ->Z=Z;
            NewZ->Accumulator=val;
            s->ZIndex=NewZ;

            return NewZ;
          } else {
            r=s->ZIndex;

            do {
              if (r->Z==Z) {
                r->Accumulator+=val;
                return r;
              } else 
                if ((r->NextIndex==(InfIndex3DZ *)NULL) || (r->NextIndex->Z>Z))
                   return AddInfIndex3DZ(r,Z,val);

              r=r->NextIndex;

            } while (1);
          }
        } else
          if ((s->NextIndex==(InfIndex3DY *)NULL) || (s->NextIndex->Y>Y)) { 
            NewY=AddInfIndex3DY(s, Y);

            NewZ=AllocInfIndex3DZ();
            NewZ->NextIndex=(InfIndex3DZ *)NULL;
            NewY->ZIndex=NewZ;
            NewZ->Z=Z;
            NewZ->Accumulator=val;

            return NewZ; 
          }
        s=s->NextIndex;
      } while (1);
    } else
      if ((p->NextIndex==(InfIndex3DX *)NULL) || (p->NextIndex->X>X)) { 
        p=AddInfIndex3DX(p, X);
  
        NewY=AllocInfIndex3DY();
        NewY->NextIndex=(InfIndex3DY *)NULL;
        p->YIndex=NewY;
        NewY->Y=Y;

        NewZ=AllocInfIndex3DZ();
        NewZ->NextIndex=(InfIndex3DZ *)NULL;
        NewY->ZIndex=NewZ;
        NewZ->Z=Z;
        NewZ->Accumulator=val;

        return NewZ; 
      }
    p=p->NextIndex;
  } while (1);
}

/*****************************************************************************/
InfIndex3DZ *IncAccu3D2(X,Y,Z,Mat,val,xval1,yval1,xval,yval)
double X,Y,Z,xval1,yval1,xval,yval;              /*** increment cell score ***/
InfMat3D *Mat;
int val;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s,*NewY;
  InfIndex3DZ *r,*NewZ;

  if (p->NextIndex==(InfIndex3DX *)NULL){
    p=AddInfIndex3DX(p, X);
  
    NewY=AllocInfIndex3DY();
    NewY->NextIndex=(InfIndex3DY *)NULL;
    p->YIndex=NewY;
    NewY->Y=Y;

    NewZ=AllocInfIndex3DZ();
    NewZ->NextIndex=(InfIndex3DZ *)NULL;
    NewY->ZIndex=NewZ;
    NewZ->Z=Z;
    NewZ->Accumulator=val;
    NewZ->Hit=1; 
    NewZ->XValue=xval;
    NewZ->YValue=yval;
    NewZ->XValue1=xval1;
    NewZ->YValue1=yval1;

    return NewZ; 
  }

  p=p->NextIndex;

  do {
    if (p->X==X) { 
      s=p->YIndex;
      do {
        if (s->Y==Y) {  
          if (s->ZIndex==(InfIndex3DZ *)NULL) { 
            NewZ=AllocInfIndex3DZ();
            NewZ->NextIndex=(InfIndex3DZ *)NULL;
            NewZ->Z=Z;
            NewZ->Accumulator=val;
            NewZ->Hit=1;
            NewZ->XValue=xval;
            NewZ->YValue=yval;
            NewZ->XValue1=xval1;
            NewZ->YValue1=yval1;
            s->ZIndex=NewZ;

            return NewZ;
          } else {
            r=s->ZIndex;

            do {
              if (r->Z==Z) {
                r->Accumulator+=val;
                r->Hit++;
                r->XValue=((r->XValue)*(r->Hit -1)+xval)/r->Hit;
                r->XValue1=((r->XValue1)*(r->Hit -1)+xval1)/r->Hit; 
                r->YValue=((r->YValue)*(r->Hit -1)+yval)/r->Hit;
                r->YValue1=((r->YValue1)*(r->Hit -1)+yval1)/r->Hit;

                return r; 
              } else 
                if ((r->NextIndex==(InfIndex3DZ *)NULL) || (r->NextIndex->Z>Z))
                  return AddInfIndex3DZ2(r,Z,val,xval1,yval1,xval,yval);

              r=r->NextIndex;

            } while (1);
          }
        } else
          if ((s->NextIndex==(InfIndex3DY *)NULL) || (s->NextIndex->Y>Y)) { 
            NewY=AddInfIndex3DY(s,Y);

            NewZ=AllocInfIndex3DZ();
            NewZ->NextIndex=(InfIndex3DZ *)NULL;
            NewY->ZIndex=NewZ;
            NewZ->Z=Z;
            NewZ->Accumulator=val;
            NewZ->Hit=1;
            NewZ->XValue=xval;
            NewZ->YValue=yval;
            NewZ->XValue1=xval1;
            NewZ->YValue1=yval1;

            return NewZ; 
          }
        s=s->NextIndex;
      } while (1);
    } else
      if ((p->NextIndex==(InfIndex3DX *)NULL) || (p->NextIndex->X>X)) { 
        p=AddInfIndex3DX(p, X);
  
        NewY=AllocInfIndex3DY();
        NewY->NextIndex=(InfIndex3DY *)NULL;
        p->YIndex=NewY;
        NewY->Y=Y;

        NewZ=AllocInfIndex3DZ();
        NewZ->NextIndex=(InfIndex3DZ *)NULL;
        NewY->ZIndex=NewZ;
        NewZ->Z=Z;
        NewZ->Accumulator=val;
        NewZ->Hit=1;
        NewZ->XValue=xval;
        NewZ->YValue=yval;
        NewZ->XValue1=xval1;
        NewZ->YValue1=yval1;

        return NewZ; 
      }
    p=p->NextIndex;

  } while (1);
}

/*****************************************************************************/
ExamineAccu3D(Mat,XMin,XMax,YMin,YMax,ZMin,ZMax,AccuMax) /*** count cells ***/
InfMat3D *Mat;
double *XMin,*XMax,*YMin,*YMax,*ZMin,*ZMax;
int *AccuMax;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  /* Infinity value checking here ??? */

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X;
  *XMax=p->X;

  s=p->YIndex;
  *YMin=s->Y;
  *YMax=s->Y;

  r=s->ZIndex;
  *ZMin=r->Z;
  *ZMax=r->Z;

  do {
    s=p->YIndex;
    if (s->Y<*YMin)
      *YMin=s->Y;

    do {
      if (s->Y>*YMax)
        *YMax=s->Y;

      r=s->ZIndex;
      if (r->Z<*ZMin)
        *ZMin=r->Z;
      do {
        if (r->Z>*ZMax)
          *ZMax=r->Z;
        if (r->Accumulator>*AccuMax)
          *AccuMax=r->Accumulator;

         r=r->NextIndex;
         n++;
       } while (r);

       s=s->NextIndex;
     } while (s);

     if (*XMin>p->X)
       *XMin=p->X;
     if (*XMax<p->X)
       *XMax=p->X;

     p=p->NextIndex;
  } while (p);

  return n;
}

/*****************************************************************************/
ExamineAccu3D2(Mat,XMin,XMax,YMin,YMax,ZMin,ZMax,AccuMax,XaMax,YaMax,ZaMax)
InfMat3D *Mat;                                            /*** count cells ***/
double *XMin,*XMax,*YMin,*YMax,*ZMin,*ZMax,*XaMax,*YaMax,*ZaMax;
int *AccuMax;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  /* Infinity value checking here ??? */

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X;
  *XMax=p->X;

  s=p->YIndex;
  *YMin=s->Y;
  *YMax=s->Y;

  r=s->ZIndex;
  *ZMin=r->Z;
  *ZMax=r->Z;

  do {
    s=p->YIndex;
    if (s->Y<*YMin)
      *YMin=s->Y;

    do {
      if (s->Y>*YMax)
        *YMax=s->Y;

      r=s->ZIndex;
      if (r->Z<*ZMin)
        *ZMin=r->Z;
      do {
        if (r->Z>*ZMax)
          *ZMax=r->Z;
        if (r->Accumulator>*AccuMax) {
          *AccuMax=r->Accumulator;
          *ZaMax=r->Z;
           *YaMax=s->Y;
          *XaMax=p->X;
        }

        r=r->NextIndex; n++;
      } while (r);
 
      s=s->NextIndex;
    } while (s);

    if (*XMin>p->X)
      *XMin=p->X;
    if (*XMax<p->X)
      *XMax=p->X;

     p=p->NextIndex;
  } while (p);

  return n;
}

/*****************************************************************************/
ExamineAccu3D3(Mat,XMin,XMax,YMin,YMax,ZMin,ZMax,AccuMax,XaMax,YaMax,ZaMax,
               XVal1,YVal1,XVal,YVal)                     /*** count cells ***/
InfMat3D *Mat;
double *XMin,*XMax,*YMin,*YMax,*ZMin,*ZMax,*XaMax,*YaMax,*ZaMax,
       *XVal1,*YVal1,*XVal,*YVal;
int *AccuMax;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  /* Infinity value checking here ??? */

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X;
  *XMax=p->X;

  s=p->YIndex;
  *YMin=s->Y;
  *YMax=s->Y;

  r=s->ZIndex;
  *ZMin=r->Z;
  *ZMax=r->Z;

  do {
    s=p->YIndex;
    if (s->Y<*YMin)
      *YMin=s->Y;

    do {
      if (s->Y>*YMax)
        *YMax=s->Y;

      r=s->ZIndex;
      if (r->Z<*ZMin)
        *ZMin=r->Z;
      do {
        if (r->Z>*ZMax)
          *ZMax=r->Z;
        if (r->Accumulator>*AccuMax) {
          *AccuMax=r->Accumulator;
          *XVal=r->XValue;
          *XVal1=r->XValue1;
          *YVal=r->YValue;
          *YVal1=r->YValue1;
          *ZaMax=r->Z;
          *YaMax=s->Y;
          *XaMax=p->X;
        }

        r=r->NextIndex; n++;
      } while (r);

      s=s->NextIndex;
    } while (s);

    if (*XMin>p->X)
      *XMin=p->X;
    if (*XMax<p->X)
      *XMax=p->X;

    p=p->NextIndex;
  } while (p);

  return n;
}


/*****************************************************************************/
ExamineAccu3D4(Mat,XMin,XMax,YMin,YMax,ZMin,ZMax,AccuMax,XaMax,YaMax,ZaMax,
               XVal1,YVal1,XVal,YVal,AccuHit)             /*** count cells ***/
InfMat3D *Mat;
double *XMin,*XMax,*YMin,*YMax,*ZMin,*ZMax,*XaMax,*YaMax,*ZaMax,
       *XVal1,*YVal1,*XVal,*YVal;
int *AccuMax,*AccuHit;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  /* Infinity value checking here ??? */

  p=p->NextIndex;

  *AccuMax=0;

  *XMin=p->X;
  *XMax=p->X;

  s=p->YIndex;
  *YMin=s->Y;
  *YMax=s->Y;

  r=s->ZIndex;
  *ZMin=r->Z;
  *ZMax=r->Z;

  do {
    s=p->YIndex;
    if (s->Y<*YMin)
      *YMin=s->Y;

    do {
      if (s->Y>*YMax)
        *YMax=s->Y;

      r=s->ZIndex;
      if (r->Z<*ZMin)
        *ZMin=r->Z;
      do {
        if (r->Z>*ZMax)
          *ZMax=r->Z;
        if (r->Accumulator>*AccuMax) {
          *AccuMax=r->Accumulator;
          *AccuHit=r->Hit;
          *XVal=r->XValue;
          *XVal1=r->XValue1;
          *YVal=r->YValue;
          *YVal1=r->YValue1;
          *ZaMax=r->Z;
          *YaMax=s->Y;
          *XaMax=p->X;
        }

        r=r->NextIndex; n++;
      } while (r);

      s=s->NextIndex;
    } while (s);

    if (*XMin>p->X)
      *XMin=p->X;
    if (*XMax<p->X)
      *XMax=p->X;

    p=p->NextIndex;
 } while (p);

  return n;
}

/*****************************************************************************/
PrintAccu3D(Mat) /*** print accumulator contents ***/
InfMat3D *Mat;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  p=p->NextIndex;

  do {
    s=p->YIndex;
    do {
      r=s->ZIndex;
      do {
        printf("(%lf,%lf,%lf)=%d\n",p->X,s->Y,r->Z,r->Accumulator);

        r=r->NextIndex;
      } while (r);

      s=s->NextIndex; n++;
    } while (s);

    p=p->NextIndex;
  } while (p);

  return n;
}

/*****************************************************************************/
PrintAccu3D2(Mat,nhits) /*** print accumulator contents ***/
InfMat3D *Mat;
int nhits;
{
  InfIndex3DX *p=Mat;
  InfIndex3DY *s;
  InfIndex3DZ *r;
  int n=0;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return 0;

  p=p->NextIndex;

  do {
    s=p->YIndex;
    do {
      r=s->ZIndex;
      do {
        if (nhits == r->Accumulator)
          printf("(%ld,%ld,%ld) ",nint(p->X),nint(s->Y),nint(r->Z));

        r=r->NextIndex;
      } while (r);

      s=s->NextIndex; n++;
    } while (s);

    p=p->NextIndex;
  } while (p);

  return n;
}

/*****************************************************************************/
RemoveAccumulatorSpace3D(Mat) /*** remove the whole accumulator ***/
InfMat3D *Mat;
{
  InfIndex3DX *p=Mat, *tmp_p;
  InfIndex3DY *s, *tmp_s;
  InfIndex3DZ *r, *tmp_r;

  if (p->NextIndex==(InfIndex3DX *)NULL)
    return -1;   /* 0  -> -1 */

  free(p);
  p=p->NextIndex;

  do {
    s=p->YIndex;
    tmp_p=p->NextIndex;
    do {
      r=s->ZIndex;
      tmp_s=s->NextIndex;
      do {
        tmp_r=r->NextIndex;
        free(r);
        r=tmp_r;
      } while (r);
      free(s);
      s=tmp_s;
    } while (s);
    free(p);
    p=tmp_p;
  } while (p);

}
