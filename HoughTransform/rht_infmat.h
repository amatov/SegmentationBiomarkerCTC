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
 * File:    ht_infmat.h
 * Purpose: header include for the RHT dynamic accumulator
 * Date:    Jun 1, 1993
 *****************************************************************************/

#ifndef TRUE
#define TRUE                    1
#endif
#ifndef FALSE
#define FALSE                   0
#endif
#ifndef NULL
#define NULL                    0
#endif

typedef union double_t
  {
    unsigned char bytes[8];
    double value;
  } Double;
     

typedef struct  InfIndex_t InfMat;

typedef struct InfIndex_t
  {
    struct Accu_t *Accu;
    struct InfIndex_t *NextIndex;
    double X_coord;
  } InfIndex;

typedef struct Accu_t
  {
    struct Accu_t *Up;
    double Y_coord;
    int Accumulator,Time;
  } Accu;


typedef struct  InfIndex3DX_t InfMat3D;

typedef struct InfIndex3DX_t
  {
    struct InfIndex3DX_t *NextIndex;
    struct InfIndex3DY_t *YIndex;
    double X;
  } InfIndex3DX;

typedef struct InfIndex3DY_t
  {
    struct InfIndex3DY_t *NextIndex;
    struct InfIndex3DZ_t *ZIndex;
    double Y;
  } InfIndex3DY;

typedef struct InfIndex3DZ_t
  {
    struct InfIndex3DZ_t *NextIndex;
    double Z;
    int Accumulator,Hit;
    double XValue1,YValue1,XValue,YValue;
  } InfIndex3DZ;


InfMat *AllocInfMat();
InfIndex *AllocInfIndex();
Accu *AllocAccu();
InfIndex *AddInfIndex();
Accu *AddAccu();
Accu *IncAccu();

/* for cfth method */
Accu *IncAccu_cfht();

/* for circle finding */
InfMat3D *AllocInfMat3D();
InfIndex3DX *AllocInfIndex3DX();
InfIndex3DY *AllocInfIndex3DY();
InfIndex3DZ *AllocInfIndex3DZ();
InfIndex3DX *AddInfIndex3DX();
InfIndex3DY *AddInfIndex3DY();
InfIndex3DZ *AddInfIndex3DZ();
InfIndex3DZ *AddInfIndex3DZ2();
InfIndex3DZ *IncAccu3D();
InfIndex3DZ *IncAccu3D2();
