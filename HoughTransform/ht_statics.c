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
 * File:    ht_statics.c
 * Purpose: detection statistics output
 * Date:    Jun 1, 1993
 *****************************************************************************/

#include "ht_hough.h"

char text[20],text_out[256];
int flag_vec[MAX_LINES];

/*****************************************************************************/
test_found_params(maxs,num_of_maxs,real_params)
/*

  Check do the maxs array include an parameters real_params array.

  Parameters :  maxs - array of accumulator maximum coordinates
                num_of_maxs - size of the maxs array
                real_params - real line parameter list

*/
double maxs[][2], real_params[][2];
int num_of_maxs;
{
  int i, j, num_of_real_params=0, size_of_params=(int)real_params[0][0];
  double rho, theta;

  for (i=0; i<size_of_params; i++)
    flag_vec[i]=0;
  size_of_params+=2;

  for(i=0; i<num_of_maxs; i++) { /* maxs list */
    rho=isinf(maxs[i][0]) ? maxs[i][1] : /* (rho, theta) --> (a, b) */
                            maxs[i][1]/sqrt(maxs[i][0]*maxs[i][0]+1.0);
    theta=atan2(1.0,-maxs[i][0]);
    for(j=2; j<size_of_params; j++) /* real params list */
      if (!flag_vec[j-2] && /* check whether the parameters are close enough */
          fabs(rho-real_params[j][0])<=real_params[1][0] &&   /* to the real */
          fabs(theta-real_params[j][1])<=real_params[1][1]) {        /* ones */
        num_of_real_params++;
        flag_vec[j-2]=1;
        break;
      }
  }
  return num_of_real_params;
}

/*****************************************************************************/
test_found_params_rho_theta(maxs,num_of_maxs,real_params)
/*

  Check whether the maxs array include parameters of the real_params array.

  Parameters :  maxs - array of accumulator maximum coordinates
                num_of_maxs - size of the maxs array
                real_params - real line parameter list

*/
double maxs[][2], real_params[][2];
int num_of_maxs;
{
  int i, j, num_of_real_params=0, size_of_params=(int)real_params[0][0];

  for (i=0; i<size_of_params; i++)
    flag_vec[i]=0;
  size_of_params+=2;

  for(i=0; i<num_of_maxs; i++) { /* maxs list */
    for(j=2; j<size_of_params; j++) /* real params list */
      if (!flag_vec[j-2] && /* test parameters */
          fabs(maxs[i][0]-real_params[j][0])<=real_params[1][0] &&
          fabs(maxs[i][1]-real_params[j][1])<=real_params[1][1]) {
        num_of_real_params++;
        flag_vec[j-2]=1;
        break;
      }
  }
  return num_of_real_params;
}

/*****************************************************************************/
TextOut_rht(size,timecount,maxs_histogram,test,totcputime,num_of_maxs,
	    num_of_real_params,acquired_real_params,false_alarms)
/*

  Print detection statistics.

  Parameters :  size - array of accumulator sizes
                timecount - index for the size array
                maxs_histogram - histogram of maxima of multible tests
                test - number of tests
                totcputime - CPU time spent in the test sequence
                num_of_maxs - total number of maxima found
                num_of_real_params - number of real lines found
                acquired_real_params - number of real lines in the image
                false_alarms - number of false alarms

*/
int size[], maxs_histogram[], test, num_of_maxs, false_alarms;
long timecount, num_of_real_params, acquired_real_params;
float totcputime;
{
  int i, j, max_size=0;
  double average_found_maxs=0.0, cumulative_size=0.0;

  sprintf(text_out,"%sCPU time: %7.2f s\n",test>1?"Average ":"",
          totcputime/test);
  out_prnt(text_out);

  for (i=0; i<MAX_LINES; i++)
    average_found_maxs+=maxs_histogram[i]*i;
  if (test>1) {
    sprintf(text_out,"Average found maxs: %7.2f\n",average_found_maxs/test);
    out_prnt(text_out);
  }

  sprintf(text_out,"Found maxs:  ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4d ",i);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"Tests/found: ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4d ",maxs_histogram[i]);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"   --\"--   %%: ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4.1lf ",
              100.0*(double)maxs_histogram[i]/(double)test);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"Found acquired num. of maxs: %4.1lf %%\n",
          100.0*(double)maxs_histogram[num_of_maxs]/(double)test);
  out_prnt(text_out);

  if (average_found_maxs>0.0 || false_alarms>0)
    sprintf(text_out,"%salse alarm rate: %4.1lf %%\n",test>1?"Average f":"F",
            (double)false_alarms*100.0/
            (average_found_maxs+(double)false_alarms));
  else
    sprintf(text_out,"%salse alarm rate: 0.0 %%\n",test>1?"Average f":"F");
  out_prnt(text_out);

  if (num_of_real_params>=0) {
    sprintf(text_out,"%sum. of real params: %7.2lf\n",test>1?"Average n":"N",
            (double)num_of_real_params/(double)test);
    out_prnt(text_out);
    sprintf(text_out,"Found acquired num. of real params: %4.1lf %%\n",
            100.0*(double)acquired_real_params/(double)test);
    out_prnt(text_out);
  }

  if (timecount>0) {
    for(i=0; i<timecount; i++) {
      cumulative_size+=size[i];
      if (size[i]>max_size)
        max_size=size[i];
    }
    sprintf(text_out,"Average accu size: ");
    sprintf(text,"%d\n",nint(cumulative_size/(double)timecount));
    strcat(text_out,text);
    out_prnt(text_out);
    sprintf(text_out,"Maximum accu size: %d\n",max_size);
    out_prnt(text_out);
  }
}

/*****************************************************************************/
TextOut_cfht(size,accu_maxs,test,totcputime,maxs,params,false_alarms)
/*

  Print detection statistics.

  Parameters :  size - array of accumulator sizes
                accu_maxs - number of accumulator maxima
                test - number of tests
                totcputime - CPU time spent in the test sequence
                maxs - array of accumulator maxima
                params - array of real line parameters
                false_alarms - number of false alarms

*/
int size, accu_maxs, test, false_alarms;
float totcputime;
double maxs[][2], params[][2];
{
  int i, j, num_of_real_params;

  sprintf(text_out,"%sCPU time: %7.2f s\n",test>1?"Average ":"",
          totcputime/test);
  out_prnt(text_out);

  sprintf(text_out,"Found maxs: %d\n",accu_maxs);
  out_prnt(text_out);

  if (accu_maxs>0 || false_alarms>0)
    sprintf(text_out,"False alarm rate: %4.1lf %%\n",
            (double)false_alarms*100.0/(double)(accu_maxs+false_alarms));
  else
    sprintf(text_out,"False alarm rate: 0.0 %%\n");
  out_prnt(text_out);

  if ((int)params[0][0]>0) {
    num_of_real_params=test_found_params(maxs,accu_maxs,params);
    sprintf(text_out,"Num. of real params: %d\n",num_of_real_params);
    out_prnt(text_out);
  }

  sprintf(text_out,"Maximum accu size: %d\n",size);
  out_prnt(text_out);
}

/*****************************************************************************/
TextOut_sht(accu_maxs,totcputime,maxs,params,false_alarms)
/*

  Print detection statistics.

  Parameters :  accu_maxs - number of accumulator maxima
                totcputime - CPU time spent in the test sequence
                maxs - array of accumulator maxima
                params - array of real line parameters
                false_alarms - number of false alarms

*/
int accu_maxs, false_alarms;
float totcputime;
double maxs[][2], params[][2];
{
  int i, j, num_of_real_params=0, size_of_params=(int)params[0][0]+2;

  sprintf(text_out,"CPU time: %7.2f s\n",totcputime);
  out_prnt(text_out);

  sprintf(text_out,"Found maxs: %d\n",accu_maxs);
  out_prnt(text_out);

  if (accu_maxs>0 || false_alarms>0)
    sprintf(text_out,"False alarm rate: %4.1lf %%\n",
            (double)false_alarms*100.0/(double)(accu_maxs+false_alarms));
  else
    sprintf(text_out,"False alarm rate: 0.0 %%\n");
  out_prnt(text_out);

  if ((int)params[0][0]>0) {
    num_of_real_params=test_found_params_rho_theta(maxs,accu_maxs,params);
    sprintf(text_out,"Num. of real params: %d\n",num_of_real_params);
    out_prnt(text_out);
  }

}

/*****************************************************************************/
TextOut_aht(totcputime,maxs,num_of_maxs,params,false_alarms)
/*

  Print detection statistics.

  Parameters :  totcputime - CPU time spent in the test sequence
                maxs - array of accumulator maxima
                num_of_maxs - number of accumulator maxima
                params - array of real line parameters
                false_alarms - number of false alarms

*/
float totcputime;
int num_of_maxs, false_alarms;
double maxs[][2], params[][2];
{
  int i, j, num_of_real_params;

  sprintf(text_out,"CPU time: %7.2f s\n",totcputime);
  out_prnt(text_out);

  sprintf(text_out,"Found maxs: %d\n",num_of_maxs);
  out_prnt(text_out);

  if (num_of_maxs>0 || false_alarms>0)
    sprintf(text_out,"False alarm rate: %4.1lf %%\n",
            (double)false_alarms*100.0/(double)(num_of_maxs+false_alarms));
  else
    sprintf(text_out,"False alarm rate: 0.0 %%\n");
  out_prnt(text_out);

  if ((int)params[0][0]>0) {
    num_of_real_params=test_found_params(maxs,num_of_maxs,params);
    sprintf(text_out,"Num. of real params: %d\n",num_of_real_params);
    out_prnt(text_out);
  }
}

/*****************************************************************************/
TextOut_dcht(maxs_histogram,test,totcputime,num_of_maxs,
	     num_of_real_params,acquired_real_params,false_alarms)
/*

  Print detection statistics.

  Parameters :  maxs_histogram - histogram of maxima of multible tests
                test - number of tests
                totcputime - CPU time spent in the test sequence
                num_of_maxs - total number of maxima found
                num_of_real_params - number of real lines found
                acquired_real_params - number of real lines in the image
                false_alarms - number of false alarms

*/
int maxs_histogram[], test, num_of_maxs, false_alarms;
long num_of_real_params,acquired_real_params;
float totcputime;
{
  int i, j;
  double average_found_maxs=0.0;

  sprintf(text_out,"%sCPU time: %7.2f s\n",test>1?"Average ":"",
          totcputime/test);
  out_prnt(text_out);

  for (i=0; i<MAX_LINES; i++)
    average_found_maxs+=maxs_histogram[i]*i;
  if (test>1) {
    sprintf(text_out,"Average found maxs: %7.2f\n",average_found_maxs/test);
    out_prnt(text_out);
  }

  sprintf(text_out,"Found maxs:  ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4d ",i);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"Tests/found: ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4d ",maxs_histogram[i]);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"   --\"--   %%: ");
  for(i=0; i<MAX_LINES; i++)
    if (maxs_histogram[i]) {
      sprintf(text,"%4.1lf ",
              100.0*(double)maxs_histogram[i]/(double)test);
      strcat(text_out,text);
    }
  strcat(text_out,"\n");
  out_prnt(text_out);

  sprintf(text_out,"Found acquired num. of maxs: %4.1lf %%\n",
          100.0*(double)maxs_histogram[num_of_maxs]/(double)test);
  out_prnt(text_out);

  if (average_found_maxs>0.0 || false_alarms>0)
    sprintf(text_out,"%salse alarm rate: %4.1lf %%\n",test>1?"Average f":"F",
            (double)false_alarms*100.0/
            (average_found_maxs+(double)false_alarms));
  else
    sprintf(text_out,"%salse alarm rate: 0.0 %%\n",test>1?"Average f":"F");
  out_prnt(text_out);

  if (num_of_real_params>=0) {
    sprintf(text_out,"%sum. of real params: %7.2lf\n",test>1?"Average n":"N",
            (double)num_of_real_params/(double)test);
    out_prnt(text_out);
    sprintf(text_out,"Found acquired num. of real params: %4.1lf %%\n",
            100.0*(double)acquired_real_params/(double)test);
    out_prnt(text_out);
  }
}
