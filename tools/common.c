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
 * Read in a 2D array given in a .csv-like file with space delimiter
 *
 *
 * @param filename      Input: name of file to read in
 * @param array         Input/Output: already allocated array to store file values in of size x_size*y_size, indexed like array[x_size*row+col]
 * @param x_size        Input: size of x-axis of array
 * @param y_size        Input: size of y-axis of array
 * @param error_message Output: error message
 * @return the error status
 */

int class_read_2d_array_from_file(
                                  char * filename,
                                  double * array,
                                  int x_size,
                                  int y_size,
                                  ErrorMsg error_message
                                  ) {
  int row,col,char_ind;
  char line_buf[_ERRORMSGSIZE_];
  char value_buf[_ERRORMSGSIZE_];
  char char_buf[2];
  FILE *fp;

  class_open(fp, filename, "r", error_message);

  // First read in each row of data from file
  row=0;
  while (fgets(line_buf, _ERRORMSGSIZE_, fp)!=NULL) {
    // Then collect each value in the row
    col=0;
    char_ind=0;
    while (line_buf[char_ind] != '\n') {
      if (line_buf[char_ind] != ' ') { // continue reading in value until it hits delimiter
        char_buf[0] = line_buf[char_ind];
        char_buf[1] = '\0';
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

  fclose(fp);

  return _SUCCESS_;
}
