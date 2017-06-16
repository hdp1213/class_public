#include "common.h"

void class_protect_sprintf(char* dest, char* tpl,...) {
  va_list args;
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
}

void class_protect_fprintf(FILE* stream, char* tpl,...) {
  va_list args;
  char dest[6000];
  va_start(args,tpl);
  vsnprintf(dest, 2048,tpl,args);
  va_end(args);
  fprintf(stream,"%s",dest);
}

void* class_protect_memcpy(void* dest, void* from, size_t sz) {
  return memcpy(dest, from,sz);
}

int get_number_of_titles(char * titlestring){
  int i;
  int number_of_titles=0;

  for (i=0; i<strlen(titlestring); i++){
    if (titlestring[i] == '\t')
      number_of_titles++;
  }
  return number_of_titles;
}

/**
 * Read in a 2D array given in a .csv file
 *
 *
 * @param filename      Input: name of file to read in
 * @param array         Input/Output: already allocated array to store file values in of size x_size*y_size, indexed like array[x_size*row+col]
 * @param x_size        Input: size of x-axis of array
 * @param y_size        Input: size of y-axis of array
 * @param error_message Output: error message
 * @return the error status
 */

int class_read_2d_array(
                        FILE * fp,
                        double * array,
                        int x_size,
                        int y_size,
                        ErrorMsg error_message
                        ) {
  int row,col,char_ind;
  char line_buf[y_size*_CSVVALUESIZE_];
  char value_buf[_CSVVALUESIZE_];
  char char_buf[2];

  // Set terminating C string character
  char_buf[1] = '\0';

  // First read in each row of data from file
  row=0;
  while (fgets(line_buf, y_size*_ERRORMSGSIZE_, fp) != NULL) {
    // Then collect each value in the row
    col=0;
    char_ind=0;
    while (line_buf[char_ind] != '\n') {
      if (line_buf[char_ind] != ',') { // continue reading in value until it hits delimiter
        char_buf[0] = line_buf[char_ind];
        strncat(value_buf, char_buf, 2);
      }
      else { // save value to array of doubles
        sscanf(value_buf, "%lf", array+x_size*row+col);
        strncpy(value_buf, "", 2);
        col++;
      }
      char_ind++;
    }

    // Collect the last value at the end of the row
    sscanf(value_buf, "%lf", array+x_size*row+col);
    strncpy(value_buf, "", 2);
    row++;
  }

  return _SUCCESS_;
}

/**
 * Read in a 1D array given in a .csv file
 *
 * The .csv file must contain as its first line the number of elements
 * in the array, and as its second line the elements of the array
 * separated by commas.
 *
 *
 * @param fp            Input: file pointer
 * @param array         Output: array of values
 * @param array_size    Output: size of array
 * @param error_message Output: error message
 * @return the error status
 */

int class_read_1d_array(
                        FILE * fp,
                        double ** array,
                        int * array_size,
                        ErrorMsg error_message
                        ) {
  char * line_buf;
  char value_buf[_CSVVALUESIZE_];
  char char_buf[2];
  int col, char_ind;
  int file_line_size;

  // Set terminating C string character and clear value_buf buffer
  // WARNING: It is sometimes full!!
  char_buf[1] = '\0';
  strncpy(value_buf, "", 2);

  // The first line contains the array size
  fscanf(fp, "%d\n", array_size);
  *array = malloc(*array_size * sizeof(double));
  class_test((*array == NULL),
       error_message,
       "Cannot allocate array\n");

  // Allocate line buffer
  file_line_size = *array_size * _CSVVALUESIZE_;
  line_buf = malloc(file_line_size * sizeof(char));
  class_test((line_buf == NULL),
       error_message,
       "Cannot allocate line_buf\n");
  
  // Get second line read into line buffer
  fgets(line_buf, file_line_size, fp);
  class_test((line_buf == NULL),
       error_message,
       "unexpected end of file");

  // Read line buffer into array
  col = 0;
  char_ind = 0;
  while (line_buf[char_ind] != '\n') {
    if (line_buf[char_ind] != ',') { // continue reading in characters
      char_buf[0] = line_buf[char_ind];
      strncat(value_buf, char_buf, 2);
    }
    else { // save value to array of doubles
      sscanf(value_buf, "%lf", *array+col);
      strncpy(value_buf, "", 2);
      col++;
    }
    char_ind++;
  }
  sscanf(value_buf, "%lf", *array+col);
  strncpy(value_buf, "", 2);
  col++;

  // printf("Read %d/%d array points\n", col, *array_size);
  free(line_buf);

  return _SUCCESS_;
}

/**
 * Read in a bicubic b-spline given in a .csv file
 *
 * The .csv file must contain the following information on separate
 * lines in the following order:
 *   * number of knots in x direction
 *   * x-coordinates of the knots, in CSV format
 *   * number of knots in y direction
 *   * y-coordinates of the knots, in CSV format
 *   * spline coefficients
 *
 *
 * @param fp            Input: file pointer
 * @param pbsp          Input/Output: pointer to bspline_2d structure
 * @param error_message Output: error message
 * @return the error status
 */

int class_read_bicubic_bspline(
                               FILE * fp,
                               struct bspline_2d * pbsp,
                               ErrorMsg error_message
                               ) {
  int n_coeffs;
  int deg = 3;

  pbsp->degree = deg;

  class_call(class_read_1d_array(fp,
                                 &(pbsp->xknots),
                                 &(pbsp->nxknots),
                                 error_message),
             error_message,
             error_message);


  class_call(class_read_1d_array(fp,
                                 &(pbsp->yknots),
                                 &(pbsp->nyknots),
                                 error_message),
             error_message,
             error_message);


  class_call(class_read_1d_array(fp,
                                 &(pbsp->coeffs),
                                 &n_coeffs,
                                 error_message),
             error_message,
             error_message);

  class_test((n_coeffs != (pbsp->nxknots-deg-1)*(pbsp->nyknots-deg-1)),
       error_message,
       "number of supplied coefficients does not match degree of spline requested");

  return _SUCCESS_;
}
