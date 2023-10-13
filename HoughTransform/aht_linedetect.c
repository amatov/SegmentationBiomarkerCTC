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
 * For further information, please notice following papers:
 *
 *	Illingworth, J. and Kittler J., The Adaptive Hough Transform,
 *	IEEE Transactions on Pattern Analysis and Machine Intelligence, vol. 9,
 *	no. 5, 1987, pp. 690-698.	
 *
 *	Princen, J., Yuen, H.K., Illingworth, J. and Kittler, J., Properties of
 *	the Adaptive Hough Transform, Proceedings of 6th Scandinavian
 *	Conference on Image Analysis, Oulu, Finland, June 19-22, 1989, 
 *	pp. 613-620. 
 *
 * File:         aht_linedetect.c
 * Purpose:      The Adaptive Hough Transform (AHT) algorithm for lines
 * Date:         Jun 1, 1993
 * Last change:  Oct 10,1995
 *****************************************************************************/

#include "ht_hough.h"
#include "ht_graphmacros.h"

#include <sys/types.h>
#include <sys/times.h>

#ifdef VISUAL_PACK
#include "xhoughtool.h"

extern quadrant quads[];
#endif
extern int StopFlag;

pic_type acc_space[ACC_MAX_SIZE][ACC_MAX_SIZE];
pic_type tmp_acc_space[ACC_MAX_SIZE][ACC_MAX_SIZE];

double maxs[MAX_LINES][2];

#define LOW_M_RANGE -1.0001 /* parameter ranges for m and c */
#define HIGH_M_RANGE 1.0001
#define LOW_C_RANGE -256.0
#define HIGH_C_RANGE 256.0

/*****************************************************************************/

aht_line(pic,pic1,gpic,rpic,dim_x,dim_y,dim_m,dim_c,real_params,NumOfMaxs,
         MinSegLen,MaxSegWidth,MaxSegGap,acc_space_sel,bin_level,Accur_m,
         Accur_c,TextInfo,ParametersOut)
/*

  The Adapative Hough Transform (AHT) for line detection.

  Parameters :	pic,pic1,gpic,rpic - input and output images
		dim_x,dim_y,dim_m,dim_c - image and accu quantization
		real_params - real line parameter list
		NumOfMaxs - number of lines to find
		MinSegLen,MaxSegWidth,MaxSegGap - line definitions
		acc_space_sel - accumulator space selector [1,2]
		bin_level - level to which the accu is binarized with respect
                            to the maximum value
		Accur_m,Accur_c - maximum accumulator resolution
		TextInfo - statistics output selector
		ParametersOut - output parameter file

*/
pic_type pic[][MAX_SIZE], pic1[][MAX_SIZE], gpic[][MAX_SIZE], rpic[][MAX_SIZE];
int dim_x,dim_y,dim_m,dim_c,NumOfMaxs,MinSegLen,MaxSegWidth,MaxSegGap,
    acc_space_sel,TextInfo;
double real_params[][2],bin_level,Accur_m,Accur_c;
char ParametersOut[256];
{
  line *start=NULL;
  int i, j, k, max_val, NoMoreLinesToFind=FALSE, false_alarms=0, iteration;
  float totcputime=0.0;
  double m, c, m_max, c_max, small_m, big_m, small_c, big_c;
  struct tms cputimebase, cputime;
    
  StopFlag=0;

  NotifyDispatch_MACRO();

  times(&cputimebase); /* clocking starts */

  for (k=0; k<NumOfMaxs && !false_alarms && !StopFlag; ) {

    small_m=LOW_M_RANGE; /* set default ranges */
    big_m=HIGH_M_RANGE;
    small_c=LOW_C_RANGE;
    big_c=HIGH_C_RANGE;
    iteration=0;

    do {
      if (adaptive_hough_transform(pic,dim_x,dim_y,acc_space,dim_m,dim_c,
                                   acc_space_sel,small_m,big_m,small_c,big_c)
                                                               <(MinSegLen*2))
        NoMoreLinesToFind=TRUE;

      /* Draw Hough space */
      DrawAccumulatorSpace_aht_MACRO(quads[3],acc_space,dim_m,dim_c,pic1,
                                     MAX_SIZE,MAX_SIZE,
                                     small_m,big_m,small_c,big_c,
                                     LOW_M_RANGE,HIGH_M_RANGE,
                                     LOW_C_RANGE,HIGH_C_RANGE);
/*
      PutAccu_gray_MACRO(quads[2],acc_space,dim_c,dim_m);
      PrintAccu_aht(acc_space,dim_m,dim_c,small_m,big_m,small_c,big_c);
*/
      /* set new parameter ranges */
      find_significant_maxima_and_analyze_acc(acc_space,dim_m,dim_c,
                           bin_level,&small_m,&big_m,&small_c,&big_c);
/*
      PrintAccu_aht(acc_space,dim_m,dim_c,small_m,big_m,small_c,big_c);
*/
      NotifyDispatch_MACRO();

      /* Draw Hough space */
/*
      PutAccu_gray_MACRO(quads[3],acc_space,dim_c,dim_m);
*/
    iteration++;

    } while ( ((fabs(big_m-small_m)>(dim_m*Accur_m)) ||
               (fabs(big_c-small_c)>(dim_c*Accur_c)))  && !NoMoreLinesToFind
              && !StopFlag && iteration<20); /* stop iterating ? */

    if (NoMoreLinesToFind || StopFlag)
      break;

    /* accmulate once again */
    adaptive_hough_transform(pic,dim_x,dim_y,acc_space,dim_m,dim_c,
                             acc_space_sel,small_m,big_m,small_c,big_c);

    /* find maximum */
    find_best_parameter_values(acc_space,&m_max,&c_max,dim_m,dim_c,
                               small_m,big_m,small_c,big_c);
/*
    printf("%lf m: %lf c: %lf acc_space_sel: %d\n",
           m_max,c_max,acc_space_sel);
    PrintAccu_aht(acc_space,dim_m,dim_c,small_m,big_m,small_c,big_c);
*/

    m=m_max;
    c=c_max;
    FollowAndMarkLineSegments_m_c(&start,&m,&c,acc_space_sel,pic,dim_x,dim_y,
                                  MinSegLen,MaxSegWidth,MaxSegGap);

    if (RemoveEdgePoints(start,pic,dim_x,dim_y)) {
      maxs[k][0]=m; /* here really a and b */
      maxs[k][1]=c;
      k++;
    } else {
      false_alarms++;
      if (k==0)
        warn_prnt("Probably wrongly selected accumulator range\n");
      break;
    }
/*
    printf("m:%lf c:%lf m_max:%lf c_max:%lf %d\n",m,c,m_max,c_max,false_alarms);
*/
    NotifyDispatch_MACRO();
  }

  times(&cputime); /*clocking ends */
  totcputime=1./60.*(cputime.tms_utime-cputimebase.tms_utime);

  if (TextInfo)
    TextOut_aht(totcputime,maxs,k,real_params,false_alarms);

  store_all_line_segments(rpic,dim_x,dim_y,start,(pic_type)OBJECT_PIX_VAL);
  store_all_maxs_line(gpic,dim_x,dim_y,maxs,k,(pic_type)OBJECT_PIX_VAL,
	ParametersOut,1.0,start);

  DeleteAllLineSegments(start);
  StopFlag=0;
}

/*****************************************************************************/
adaptive_hough_transform(pic,dim_x,dim_y,acc_space,dim_m,dim_c,acc_space_sel,
                         small_m,big_m,small_c,big_c)
/*
  Hough transformation from pic to acc_space, using line representation

          y = m*x + c.

  Original picture dimension is given by dim_x,dim_y and parameter
  space dimension is given by dim_m, dim_c. Original picture is
  binary picture where pixel value OBJECT_PIX_VAL represent object,
  other values are background

  Parameters : pic - pic_type matrix, size [][MAX_SIZE]
               dim_x,dim_y - dimensions of the original picture
               acc_space - pic_type matrix , size [][ACC_MAX_SIZE]
               dim_m,dim_c - dimensions of the accumulator space
               acc_space_sel - parameter values selector
               small_m,big_m,small_c,big_c - parameter value ranges


*/

pic_type pic[][MAX_SIZE],acc_space[][ACC_MAX_SIZE];
int dim_x,dim_y,dim_m,dim_c,acc_space_sel;
double small_m,big_m,small_c,big_c;
{
  register i,j,k,l,pix;
  double m,c,m_bin_l,m_bin_h,c_bin_l,c_bin_h;

  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      acc_space[i][j]=tmp_acc_space[i][j]=0;

  if (acc_space_sel==1) { /* using (m, c) */
    for (i=pix=0; (i<dim_y) && !StopFlag; i++)
      for (j=0; (j<dim_x) && !StopFlag; j++)
        if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
          pix++;
	  for (k=0; k<dim_m+1; k++) { /* one c per each m */
            m=(((big_m-small_m)*k)/(double)dim_m)+small_m;
            c=(double)i-(double)j*m;
	    for (l=0; l<dim_c; l++) {
              c_bin_l=(((big_c-small_c)*l)/(double)dim_c)+small_c;
              c_bin_h=(((big_c-small_c)*(l+1))/(double)dim_c)+small_c;
              if (c>=c_bin_l && c<=c_bin_h) {
                if (k<dim_m) tmp_acc_space[k][l]=1;
                if (k>0) tmp_acc_space[k-1][l]=1;
              }
            }
          }
	  for (l=0; l<dim_c+1; l++) { /* one m per each c */
            c=(((big_c-small_c)*l)/(double)dim_c)+small_c;
            m=((double)i-c)/(double)j;
	    for (k=0; k<dim_m; k++) {
              m_bin_l=(((big_m-small_m)*k)/(double)dim_m)+small_m;
              m_bin_h=(((big_m-small_m)*(k+1))/(double)dim_m)+small_m;
              if (m>=m_bin_l && m<=m_bin_h) {
                if (l<dim_c) tmp_acc_space[k][l]=1;
                if (l>0) tmp_acc_space[k][l-1]=1;
              }
            }
          }
          for (k=0; k<dim_m; k++) /* actual accumulation */
            for (l=0; l<dim_c; l++) {
              acc_space[k][l]+=tmp_acc_space[k][l];
              tmp_acc_space[k][l]=0;
            }
        }
  } else /* using (m', c') */
    for (i=pix=0; (i<dim_y) && !StopFlag; i++)
      for (j=0; (j<dim_x) && !StopFlag; j++)
        if (pic[i][j]==(pic_type)OBJECT_PIX_VAL) {
          pix++;
	  for (k=0; k<dim_m+1; k++) { /* one c per each m */
            m=(((big_m-small_m)*k)/(double)dim_m)+small_m;
            c=(double)j-(double)i*m;
	    for (l=0; l<dim_c; l++) {
              c_bin_l=(((big_c-small_c)*l)/(double)dim_c)+small_c;
              c_bin_h=(((big_c-small_c)*(l+1))/(double)dim_c)+small_c;
              if (c>=c_bin_l && c<=c_bin_h) {
                if (k<dim_m) tmp_acc_space[k][l]=1;
                if (k>0) tmp_acc_space[k-1][l]=1;
              }
            }
          }
          for (l=0; l<dim_c+1; l++) { /* one m per each c */
            c=(((big_c-small_c)*l)/(double)dim_c)+small_c;
            m=((double)j-c)/(double)i;
	    for (k=0; k<dim_m; k++) {
              m_bin_l=(((big_m-small_m)*k)/(double)dim_m)+small_m;
              m_bin_h=(((big_m-small_m)*(k+1))/(double)dim_m)+small_m;
              if (m>=m_bin_l && m<=m_bin_h) {
                if (l<dim_c) tmp_acc_space[k][l]=1;
                if (l>0) tmp_acc_space[k][l-1]=1;
              }
            }
          }
          for (k=0; k<dim_m; k++) /* actual accumulation */
            for (l=0; l<dim_c; l++) {
              acc_space[k][l]+=tmp_acc_space[k][l];
              tmp_acc_space[k][l]=0;
            }
        }

  return pix;
}

/*****************************************************************************/
find_significant_maxima_and_analyze_acc(acc_space,dim_m,dim_c,
                                        level,small_m,big_m,small_c,big_c)
/*

  Find the most significant maximum value of the accumulator and analyze
  the accumulator to set the new parameter ranges.

  Parameters :	acc_space - accumulator array
		dim_m,dim_c - accumulator space quantization
		level - binarizing level
		small_m,big_m,small_c,big - parameter ranges

*/

pic_type acc_space[][ACC_MAX_SIZE];
int dim_m,dim_c;
double level,*small_m,*big_m,*small_c,*big_c;
{
  int i,j,k,max=0,label=1,label_list[1000][3],low_m=MAXINT,high_m=-1,low_c=MAXINT,
      high_c=-1,count=0;
  double center_m,center_c,new_small_m,new_big_m,new_small_c,new_big_c,
         hi_score=0.0;


  /* find max value */

  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j]>max) {
        max=acc_space[i][j];
        tmp_acc_space[i][j]=acc_space[i][j];
      }


  /* binarize accumulator */

  level*=max;
  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j]>=level)
        acc_space[i][j]=1;
      else
        acc_space[i][j]=0;


  /* label accumulator cells */

  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j])
        if ((i-1)>=0 && acc_space[i-1][j])
          if ((j-1)>=0 && acc_space[i][j-1]) {
            acc_space[i][j]=acc_space[i-1][j];
            if (acc_space[i-1][j]!=acc_space[i][j-1]) {
              label_list[count][0]=acc_space[i-1][j];
              label_list[count][1]=acc_space[i][j-1];
              count++;
            }
          } else
            acc_space[i][j]=acc_space[i-1][j];
        else
          if ((j-1)>=0 && acc_space[i][j-1])
            acc_space[i][j]=acc_space[i-1][j];
          else
            acc_space[i][j]=label++;


  /* resolve equivalent labels */

  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j])
        for (k=count-1; k>=0; k--)
          if (acc_space[i][j]==label_list[k][1])
            acc_space[i][j]=label_list[k][0];
/*
  PrintAccu_aht(acc_space,dim_m,dim_c,*small_m,*big_m,*small_c,*big_c);
*/
  /* find the most significant maxima */

  max=dim_m*dim_c+1;
  for (k=count=0; k<max; k++) /* zero label list */
    label_list[k][0]=label_list[k][1]=label_list[k][2]=0;

  for (i=count=0; i<dim_m; i++) /* make histogram of labels */
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j]) {
        for (k=label=0; k<count; k++)
          if (acc_space[i][j]==label_list[k][0]) {
            label_list[k][1]++;
            label_list[k][2]+=tmp_acc_space[i][j];
            label=1;
          }
        if (!label) {
          label_list[count][0]=acc_space[i][j];
          label_list[count][1]=1;
          label_list[count][2]=tmp_acc_space[i][j];
          count++;
        }
      }

  for (k=0; k<count; k++) /* find the most occurrent label */
    if (((double)label_list[k][2]/(double)label_list[k][1])>hi_score) {
      label=label_list[k][0];
      hi_score=(double)label_list[k][2]/(double)label_list[k][1];
    }

  for (i=0; i<dim_m; i++) /* find the cells of the best max */
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j]==label) {
        if (i<low_m) low_m=i;
        if (j<low_c) low_c=j;
        if (i>high_m) high_m=i;
        if (j>high_c) high_c=j;
      }

  center_m=((double)(high_m-low_m)/2.0)+low_m; /* find center of the max */
  center_c=((double)(high_c-low_c)/2.0)+low_c;

  if ((high_m-low_m)<2) /* special conditions */
    if ((high_m-low_m)==0) {
      low_m--; high_m++;
    } else
      if (low_m>(dim_m-high_m)) low_m--;
      else high_m++;

  if ((high_c-low_c)<2)
    if ((high_c-low_c)==0) {
      low_c--; high_c++;
    } else
      if (low_c>(dim_c-high_c)) low_c--;
      else high_c++;


  /* set new parameter dimensions */

  if (low_m>0 && high_m<(dim_m-1)) {
    new_small_m=((((*big_m)-(*small_m))*low_m)/(double)dim_m)+(*small_m);
    new_big_m=((((*big_m)-(*small_m))*(high_m+1))/(double)dim_m)+(*small_m);
  } else
    if (low_m<=0 && high_m>=(dim_m-1)) {
      new_small_m=*small_m;
      new_big_m=*big_m;
    } else {
      new_small_m=((((*big_m)-(*small_m))*(center_m-((double)dim_m/2.0)))/
                  (double)dim_m)+(*small_m);
      new_big_m=((((*big_m)-(*small_m))*(center_m-((double)dim_m/2.0)))/
                (double)dim_m)+(*big_m);
    }

  if (low_c>0 && high_c<(dim_c-1)) {
    new_small_c=((((*big_c)-(*small_c))*low_c)/(double)dim_c)+(*small_c);
    new_big_c=((((*big_c)-(*small_c))*(high_c+1))/(double)dim_c)+(*small_c);
  } else
    if (low_c<=0 && high_c>=(dim_c-1)) {
      new_small_c=*small_c;
      new_big_c=*big_c;
    } else {
      new_small_c=((((*big_c)-(*small_c))*(center_c-((double)dim_c/2.0)))/
                  (double)dim_c)+(*small_c);
      new_big_c=((((*big_c)-(*small_c))*(center_c-((double)dim_c/2.0)))/
                (double)dim_c)+(*big_c);
    }

  *small_m=new_small_m<LOW_M_RANGE?LOW_M_RANGE:new_small_m;
  *big_m=new_big_m>HIGH_M_RANGE?HIGH_M_RANGE:new_big_m;
  *small_c=new_small_c<LOW_C_RANGE?LOW_C_RANGE:new_small_c;
  *big_c=new_big_c>HIGH_C_RANGE?HIGH_C_RANGE:new_big_c;

}

/*****************************************************************************/
find_best_parameter_values(acc_space,m_max,c_max,dim_m,dim_c,
                           small_m,big_m,small_c,big_c)
/*

  Find the maximum value of the accumulator array and set
  the corresponding parameter values.

  Parameters :	acc_space - accumulator array
		m_max,c_max - parameter values
		dim_m,dim_c - accumulator space quantization
		small_m,big_m,small_c,big - parameter ranges

*/
pic_type acc_space[][ACC_MAX_SIZE];
double *m_max,*c_max,small_m,big_m,small_c,big_c;
int dim_m,dim_c;
{
  int i,j,max=0,m,c;

  for (i=0; i<dim_m; i++)
    for (j=0; j<dim_c; j++)
      if (acc_space[i][j]>max) {
        max=acc_space[i][j];
        m=i;
        c=j;
      }

  *m_max=(((big_m-small_m)*((double)m+.5))/(double)dim_m)+small_m;
  *c_max=(((big_c-small_c)*((double)c+.5))/(double)dim_c)+small_c;
}

/*****************************************************************************/
PrintAccu_aht(acc_space,dim_m,dim_c,small_m,big_m,small_c,big_c)
/*

  Print accumulator array.

  Parameters :	acc_space - accumulator array
		dim_m,dim_c - accumulator space quantization
		small_m,big_m,small_c,big - parameter ranges

*/
pic_type acc_space[][ACC_MAX_SIZE];
int dim_m,dim_c;
double small_m,big_m,small_c,big_c;
{
  int i,j;

  printf("%lf\t",big_c);

  for (j=dim_c-1; j>=0; j--) {
    for (i=0; i<dim_m; i++)
      printf("%3d ",acc_space[i][j]);
    if (j==1)
      printf("\n%lf\t",small_c);
    else
      printf("\n\t\t");
  }

  printf("%lf\t\t\t%lf\n\n",small_m,big_m);
}
