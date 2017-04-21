#include "common.h"
#include "arrays.h"

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

  char spline_file[_CSVVALUESIZE_];
  char bound_file[_CSVVALUESIZE_];
  char out_file[_CSVVALUESIZE_];
  FILE * ofp;

  int i;

  // double * xknots;
  // int n_xknots;
  // double * yknots;
  // int n_yknots;
  // double * coeffs;

  FILE * fp;

  struct bspline_2d bsp;
  struct bspline_2d * pbsp;

  double * x_eval;
  int mx = 21;

  double * y_eval;
  int my = 21;

  double * z;
  ErrorMsg error_message;

  pbsp = &bsp;

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
  class_open(fp, spline_file, "r", error_message);
  class_call(class_read_bicubic_bspline(fp,
                                        pbsp,
                                        error_message),
             error_message,
             error_message);
  fclose(fp);

  // Evaluate spline along x_eval and y_eval
  printf("Evaluating spline...\n");
  class_call(array_eval_bicubic_bspline(pbsp,
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
  free(pbsp->xknots);
  free(pbsp->yknots);
  free(pbsp->coeffs);

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
  int deg = 3;
  
  class_open(fp, filename, "r", error_message);

  class_call(class_read_1d_array(fp,
                                 xknots,
                                 nx,
                                 error_message),
             error_message,
             error_message);


  class_call(class_read_1d_array(fp,
                                 yknots,
                                 ny,
                                 error_message),
             error_message,
             error_message);


  class_call(class_read_1d_array(fp,
                                 coeffs,
                                 &n_coeffs,
                                 error_message),
             error_message,
             error_message);

  class_test((n_coeffs != (*nx-deg-1)*(*ny-deg-1)),
       error_message,
       "number of supplied coefficients is incorrect");

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
  *yeval = malloc(ny * sizeof(double));
  class_test((*yeval == NULL),
       error_message,
       "Cannot allocate yeval\n");

  for (i = 0; i < ny; ++i) {
    *(*yeval+i) = y_min + (double) (i)/(double) (ny-1) * (y_max - y_min);
  }

  fclose(fp);

  return _SUCCESS_;
}
