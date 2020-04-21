/**
* @file asciireader.hpp
* @brief Class that reads an ascii file into memory to enable parallel parsing.
* @author Aurel Neic
* @version
* @date 2018-09-18
*/

#include <stdio.h>
#include <errno.h>
#include <string.h>

#include <iostream>
#include <string>
#include <vector>



/// Class that reads an ascii file into memory to enable parallel parsing.
class asciireader {

  public:
  /**
  * @brief Read a whole file into memory
  *
  * @param file  The file to read.
  */
  void read_file(std::string file, bool exit_on_error = true)
  {
    FILE* fd = fopen(file.c_str(), MT_FOPEN_READ);
    if(!fd) {
      this->treat_file_error(file, errno, exit_on_error);
      return;
    }

    size_t filesize = get_file_size(fd);
    bytebuffer.resize(filesize);

    // read the entire textfile into the byte buffer
    fread(bytebuffer.data(), filesize, 1, fd);
    fclose(fd);

    // we do not know how many lines the file holds, we approximate
    // 50 bytes per line
    buff_dsp.reserve(filesize / 50);

    int line_cnt = 0;
    buff_dsp.push_back(0);
    size_t idx = 0, dsp_idx= 0;

    while(idx < filesize) {
      line_cnt++;

      if(bytebuffer[idx] == '\n') {
        dsp_idx++;
        buff_dsp.push_back(buff_dsp[dsp_idx - 1] + line_cnt);
        line_cnt = 0;
      }

      idx++;
    }
  }


  /**
  * @brief Return the number of read lines.
  *
  * @return Number of read lines.
  */
  size_t num_lines()
  {
    if(buff_dsp.size() == 0) return 0;
    else                     return buff_dsp.size() - 1;
  }

  /**
  * @brief Copy a line into a given buffer.
  *
  * @param line_nbr  The line number to copy.
  * @param buff      The provided buffer.
  * @param buffsize  The size of the provided buffer.
  *
  * @return
  */
  bool get_line(size_t line_nbr, char* buff, size_t buffsize)
  {
    if(bytebuffer.size() == 0)
      return false;

    if(line_nbr > (bytebuffer.size() - 1))
      return false;

    size_t linesize = buff_dsp[line_nbr + 1] - buff_dsp[line_nbr];

    if(buffsize < (linesize + 1))
      return false;

    for(size_t i=0; i<linesize; i++)
      buff[i] = bytebuffer[buff_dsp[line_nbr]+i];

    buff[linesize] = '\0';

    return true;
  }

  private:
  std::vector<unsigned char> bytebuffer;   ///< The buffer holding the read file.
  std::vector<size_t>        buff_dsp;     ///< The displacement between lines.

  /// Treat a file read error by parsing the errno error number.
  void treat_file_error(const std::string & file, int error_number, bool do_exit = true)
  {
    fprintf(stderr, "Error: could not open %s for reading! Reason:\n"
                    "%s\n", file.c_str(), strerror(error_number));
    if(do_exit)
      exit(1);
  }

  /// Get the size of the file behind the file-descriptor fd
  size_t get_file_size(FILE* fd)
  {
    size_t oldpos = ftell(fd);
    fseek(fd, 0L, SEEK_END);
    size_t sz = ftell(fd);
    fseek(fd, oldpos, SEEK_SET);

    return sz;
  }
};


