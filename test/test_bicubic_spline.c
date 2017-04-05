#include "common.h"
#include "arrays.h"

#define _FILEVALUESIZE_ 64

int fill_spline_info(char *filename,
                     double **xknots,
                     int *nx,
                     double **yknots,
                     int *ny,
                     double **coeffs,
                     ErrorMsg error_message);

int fill_eval_info(char *filename,
                   double **xeval,
                   int nx,
                   double **yeval,
                   int ny,
                   ErrorMsg error_message);

int main(int argc, char const *argv[])
{

  char spline_file[_FILEVALUESIZE_];
  char bound_file[_FILEVALUESIZE_];
  char out_file[_FILEVALUESIZE_];
  FILE * ofp;

  int i;

  double * xknots;
  int n_xknots;
  double * yknots;
  int n_yknots;
  double * coeffs;

  double * x_eval;
  int mx = 21;

  double * y_eval;
  int my = 21;

  double * z;
  ErrorMsg error_message;

  // Read in command line arguments
  if (argc > 1) {
    sscanf(argv[1], "%d", &mx);
    sscanf(argv[2], "%d", &my);
  }

  if (argc > 3) {
    strcpy(bound_file, argv[3]);
  }
  else {
    strcpy(bound_file, "spline.bnds");
  }

  if (argc > 4) {
    strcpy(spline_file, argv[4]);
  }
  else {
    strcpy(spline_file, "spline.info");
  }

  
  // Allocate new evaluation points
  printf("Reading in bounds file '%s'...\n", bound_file);
  class_call(fill_eval_info(bound_file,
                            &x_eval,
                            mx,
                            &y_eval,
                            my,
                            error_message),
             error_message,
             error_message);

  // Allocate results array
  z = malloc((mx*my) * sizeof(double));
  if (z == NULL) {
    printf("%s(L:%d) Cannot allocate z\n",__func__,__LINE__);
    return _FAILURE_;
  }

  // Read in spline file
  printf("Reading in spline file '%s'...\n", spline_file);
  class_call(fill_spline_info(spline_file,
                              &xknots,
                              &n_xknots,
                              &yknots,
                              &n_yknots,
                              &coeffs,
                              error_message),
             error_message,
             error_message);

  // Evaluate spline along x_eval and y_eval
  printf("Evaluating spline...\n");
  class_call(array_eval_bicubic_bspline(xknots,
                                        n_xknots,
                                        yknots,
                                        n_yknots,
                                        coeffs,
                                        x_eval,
                                        mx,
                                        y_eval,
                                        my,
                                        z,
                                        error_message),
             error_message,
             error_message);

  free(x_eval);
  free(y_eval);
  free(xknots);
  free(yknots);
  free(coeffs);

  // Write results to an output file
  strcpy(out_file, "spline.out");
  printf("Writing results to '%s'...\n", out_file);
  ofp = fopen(out_file, "w");

  for (i = 0; i < mx*my; ++i) {
    if (i != mx*my-1) {
      fprintf(ofp, "%.24e,", z[i]);
    }
    else {
      fprintf(ofp, "%.24e", z[i]);
    }
  }
  fprintf(ofp, "\n");

  fclose(ofp);
  free(z);

  return _SUCCESS_;
}


int fill_spline_info(
                     char *filename,
                     double **xknots,
                     int *nx,
                     double **yknots,
                     int *ny,
                     double **coeffs,
                     ErrorMsg error_message
                     ) {
  FILE *fp;
  int n_coeffs;
  int col, char_ind;
  int file_line_size;
  char * line_buf;
  char value_buf[_FILEVALUESIZE_];
  char char_buf[2];
  int deg = 3;
  
  // Set terminating C string character
  char_buf[1] = '\0';

  class_open(fp, filename, "r", error_message);

  // The first line contains the number of xknots, nx
  fscanf(fp, "%d\n", nx);
  *xknots = malloc(*nx * sizeof(double));
  class_test((*xknots == NULL),
       error_message,
       "Cannot allocate xknots\n");

  // Allocate line buffer
  file_line_size = *nx * _FILEVALUESIZE_;
  line_buf = malloc(file_line_size * sizeof(char));
  class_test((line_buf == NULL),
       error_message,
       "Cannot allocate line_buf\n");
  
  // Get second line read into line buffer
  fgets(line_buf, file_line_size, fp);
  class_test((line_buf == NULL),
       error_message,
       "unexpected end of file");

  // Read second line into xknots
  col = 0;
  char_ind = 0;
  while (line_buf[char_ind] != '\n') {
    if (line_buf[char_ind] != ',') { // continue reading in characters
      char_buf[0] = line_buf[char_ind];
      strncat(value_buf, char_buf, 2);
    }
    else { // save value to array of doubles
      sscanf(value_buf, "%lf", *xknots+col);
      strncpy(value_buf, "", 2);
      col++;
    }
    char_ind++;
  }
  sscanf(value_buf, "%lf", *xknots+col);
  strncpy(value_buf, "", 2);
  col++;

  free(line_buf);
  // printf("Read %d/%d xknots\n", col, *nx);


  // The third line contains the number of yknots, ny
  fscanf(fp, "%d\n", ny);
  *yknots = malloc(*ny * sizeof(double));
  class_test((*yknots == NULL),
       error_message,
       "Cannot allocate yknots");

  // Allocate line buffer
  file_line_size = *ny * _FILEVALUESIZE_;
  line_buf = malloc(file_line_size * sizeof(char));
  class_test((line_buf == NULL),
       error_message,
       "Cannot allocate line_buf\n");

  // Get fourth line read into line buffer
  fgets(line_buf, file_line_size, fp);
  class_test((line_buf == NULL),
       error_message,
       "unexpected end of file");

  // Read fourth line into yknots
  col = 0;
  char_ind = 0;
  while (line_buf[char_ind] != '\n') {
    if (line_buf[char_ind] != ',') { // continue reading in characters
      char_buf[0] = line_buf[char_ind];
      strncat(value_buf, char_buf, 2);
    }
    else { // save value to array of doubles
      sscanf(value_buf, "%lf", *yknots+col);
      strncpy(value_buf, "", 2);
      col++;
    }
    char_ind++;
  }
  sscanf(value_buf, "%lf", *yknots+col);
  strncpy(value_buf, "", 2);
  col++;

  free(line_buf);
  // printf("Read %d/%d yknots\n", col, *ny);


  // The fifth line contains the number of coeffs, n_coeffs
  fscanf(fp, "%d\n", &n_coeffs);
  class_test((n_coeffs != (*nx-deg-1)*(*ny-deg-1)),
       error_message,
       "number of supplied coefficients is incorrect");

  *coeffs = malloc(n_coeffs * sizeof(double));
  class_test((*coeffs == NULL),
       error_message,
       "Cannot allocate coeffs");

  // Allocate line buffer
  file_line_size = n_coeffs * _FILEVALUESIZE_;
  line_buf = malloc(file_line_size * sizeof(char));
  class_test((line_buf == NULL),
       error_message,
       "Cannot allocate line_buf\n");

  // Get sixth line read into line buffer
  fgets(line_buf, file_line_size, fp);
  class_test((line_buf == NULL),
       error_message,
       "unexpected end of file");

  // Read sixth line into coeffs
  col = 0;
  char_ind = 0;
  while (line_buf[char_ind] != '\n') {
    if (line_buf[char_ind] != ',') { // continue reading in characters
      char_buf[0] = line_buf[char_ind];
      strncat(value_buf, char_buf, 2);
    }
    else { // save value to array of doubles
      sscanf(value_buf, "%lf", *coeffs+col);
      strncpy(value_buf, "", 2);
      col++;
    }
    char_ind++;
  }
  sscanf(value_buf, "%lf", *coeffs+col);
  strncpy(value_buf, "", 2);
  col++;

  free(line_buf);
  // printf("Read %d/%d coeffs\n", col, n_coeffs);

  fclose(fp);

  return _SUCCESS_;
}


int fill_eval_info(char *filename,
                   double **xeval,
                   int nx,
                   double **yeval,
                   int ny,
                   ErrorMsg error_message
                   ) {
  FILE *fp;
  int i;

  double x_min, x_max;
  double y_min, y_max;

  class_open(fp, filename, "r", error_message);

  // The first line contains x_min and x_max
  fscanf(fp, "%lf,%lf\n", &x_min, &x_max);
  *xeval = malloc(nx * sizeof(double));
  class_test((*xeval == NULL),
       error_message,
       "Cannot allocate xeval\n");


  for (i = 0; i < nx; ++i) {
    *(*xeval+i) = x_min + (double) (i)/(double) (nx-1) * (x_max - x_min);
  }

  // The second line contains y_min and y_max
  fscanf(fp, "%lf,%lf\n", &y_min, &y_max);
  *yeval = malloc(nx * sizeof(double));
  class_test((*yeval == NULL),
       error_message,
       "Cannot allocate yeval\n");

  for (i = 0; i < ny; ++i) {
    *(*yeval+i) = y_min + (double) (i)/(double) (ny-1) * (y_max - y_min);
  }

  fclose(fp);

  return _SUCCESS_;
}
