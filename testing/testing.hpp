
#ifndef _TESTING
#define _TESTING

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

#include <string>
#include <vector>
#include <list>


namespace tst
{

void read_text_line(FILE* fd,
                    char* buff,
                    const size_t buffsize,
                    std::string & line)
{
  char* ptr = fgets(buff, buffsize, fd);
  if(ptr) line = ptr;
  else    line = "";
}

int read_float(FILE* fd, double & val)
{
  int ret = fscanf(fd, "%lf ", &val);
  return ret;
}

size_t file_size(FILE* fd) {

  size_t oldpos = ftell(fd);
  fseek(fd, 0L, SEEK_END);
  size_t sz = ftell(fd);
  fseek(fd, oldpos, SEEK_SET);

  return sz;
}

void remove(std::string & str, const char c)
{
  size_t widx=0, ridx=0;
  for(ridx=0; ridx < str.size(); ridx++)
    if(str[ridx] != c) str[widx++] = str[ridx];

  str.resize(widx);
}

void remove_multiple(std::string & str, const char c)
{
  size_t widx=0, ridx=0;
  for(ridx=0; ridx < str.size(); ridx++) {
    if(str[ridx] != c)
      str[widx++] = str[ridx];
    else {
      // copy first one
      str[widx++] = str[ridx++];
      // jump over all remainig c's
      while(str[ridx] == c) ridx++;
    }
  }

  str.resize(widx);
}

void remove_indent(std::string & str, const char c)
{
  size_t widx=0, ridx=0;
  while(str[ridx] == c) ridx++;

  for(; ridx < str.size(); ridx++)
    str[widx++] = str[ridx];

  str.resize(widx);
}




class base_check
{
  protected:
  std::string _name;

  public:
  base_check(const char* name) : _name(name)
  {}

  const char* name() const
  {
    return _name.c_str();
  }

  virtual void evaluate(bool & fail, std::string & reason) const = 0;
};

class bin_file_equal : public base_check
{
  private:
  std::string _lhs_file;
  std::string _rhs_file;

  public:
  bin_file_equal(const char* name,
                 const char* lhs_file,
                 const char* rhs_file) :
                 base_check(name),
                 _lhs_file(lhs_file),
                 _rhs_file(rhs_file)
  {}

  void evaluate(bool & fail, std::string & reason) const
  {
    const size_t buffsize = 2048;
    char buffer[buffsize];

    fail = false;
    reason = "";

    FILE* lhs_fd = fopen(_lhs_file.c_str(), "rb");
    if(lhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _lhs_file.c_str(), strerror(errno));
      reason = buffer;
      return;
    }

    FILE* rhs_fd = fopen(_rhs_file.c_str(), "rb");
    if(rhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _rhs_file.c_str(), strerror(errno));
      reason = buffer;
      fclose(lhs_fd);
      return;
    }

    size_t lhs_size = tst::file_size(lhs_fd);
    size_t rhs_size = tst::file_size(rhs_fd);
    if(lhs_size != rhs_size) {
      fail = true;
      reason = "Files not of same size.";
      fclose(lhs_fd); fclose(rhs_fd);
      return;
    }

    std::vector<char> lhs_buff(lhs_size);
    std::vector<char> rhs_buff(lhs_size);
    fread(lhs_buff.data(), 1, lhs_size, lhs_fd);
    fread(rhs_buff.data(), 1, lhs_size, rhs_fd);
    fclose(lhs_fd); fclose(rhs_fd);

    for(size_t i=0; i<lhs_size; i++)
      if(lhs_buff[i] != rhs_buff[i]) {
        fail = true;
        reason = "binary missmatch.";
        return;
      }
  }
};

class text_file_equal : public base_check
{
  private:
  std::string _lhs_file;
  std::string _rhs_file;

  public:

  text_file_equal(const char* name,
                  const char* lhs_file,
                  const char* rhs_file) :
                  base_check(name),
                  _lhs_file(lhs_file),
                  _rhs_file(rhs_file)
  {}

  void evaluate(bool & fail, std::string & reason) const
  {
    const size_t buffsize = 2048;
    char buffer[buffsize];

    fail = false;
    reason = "";

    FILE* lhs_fd = fopen(_lhs_file.c_str(), "rb");
    if(lhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _lhs_file.c_str(), strerror(errno));
      reason = buffer;
      return;
    }

    FILE* rhs_fd = fopen(_rhs_file.c_str(), "rb");
    if(rhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _rhs_file.c_str(), strerror(errno));
      reason = buffer;
      fclose(lhs_fd);
      return;
    }

    std::string lhs_line, rhs_line;
    size_t line = 0;

    do {
      // read lines
      tst::read_text_line(lhs_fd, buffer, buffsize, lhs_line);
      tst::read_text_line(rhs_fd, buffer, buffsize, rhs_line);

      // remove indentation from lines
      tst::remove_indent(lhs_line, ' ');
      tst::remove_indent(rhs_line, ' ');
      // remove multiple spaces from lines
      tst::remove_multiple(lhs_line, ' ');
      tst::remove_multiple(rhs_line, ' ');

      // compare lines
      if(lhs_line.compare(rhs_line) != 0)
      {
        sprintf(buffer, "Missmatch at line %ld.", line);
        fail = true;
        reason = buffer;
        break;
      }
      line++;
    }
    while(lhs_line.size());

    fclose(lhs_fd);
    fclose(rhs_fd);
  }
};


class float_text_file_equal : public base_check
{
  private:
  std::string _lhs_file;
  std::string _rhs_file;
  double      _prec;

  public:

  float_text_file_equal(const char* name,
                  const char* lhs_file,
                  const char* rhs_file,
                  const double prec) :
                  base_check(name),
                  _lhs_file(lhs_file),
                  _rhs_file(rhs_file),
                  _prec(prec)
  {}

  void evaluate(bool & fail, std::string & reason) const
  {
    const size_t buffsize = 2048;
    char buffer[buffsize];

    fail = false;
    reason = "";

    FILE* lhs_fd = fopen(_lhs_file.c_str(), "rb");
    if(lhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _lhs_file.c_str(), strerror(errno));
      reason = buffer;
      return;
    }

    FILE* rhs_fd = fopen(_rhs_file.c_str(), "rb");
    if(rhs_fd == NULL) {
      fail = true;
      sprintf(buffer, "%s cannot be opened: %s", _rhs_file.c_str(), strerror(errno));
      reason = buffer;
      fclose(lhs_fd);
      return;
    }

    double lhs_val, rhs_val;
    int lr, rr;

    do {
      // read lines
      lr = read_float(lhs_fd, lhs_val);
      rr = read_float(rhs_fd, rhs_val);

      // compare lines
      if(lr != rr) {
        sprintf(buffer, "Number of values dont match");
        fail = true;
        reason = buffer;
        break;
      }
      if(lhs_val - rhs_val > _prec) {
        sprintf(buffer, "Values dont match: %lf != %lf", lhs_val, rhs_val);
        fail = true;
        reason = buffer;
        break;
      }
    }
    while(lr > 0);

    fclose(lhs_fd);
    fclose(rhs_fd);
  }
};



class execution_test
{
  private:
  std::list<base_check*> checks;
  std::string exec_command;

  public:

  execution_test(const char* command) : exec_command(command)
  {
    #ifdef TST_NO_EXEC_OUT
    exec_command += " 1> /dev/null 2> /dev/null";
    #endif
  }

  ~execution_test() {
    for(base_check* c : checks)
      delete c;
  }

  bool run() {

    bool had_errors = false;
    printf("\n\n#### Executing command \n\n"
           "%s\n\n", exec_command.c_str());

    int ret = system(exec_command.c_str());
    if(ret) {
      printf("\n## Run failed with exit status %d.\n", ret);
      had_errors = true;
      return had_errors;
    }
    else printf("\n## Run succeeded.\n");

    for(const base_check* c : checks) {
      bool failed = false;
      std::string fail_reason;

      printf("\n## Evaluating check \"%s\" : ", c->name());
      c->evaluate(failed, fail_reason);

      if(failed) {
        printf("failed. Reason: %s\n", fail_reason.c_str());
        had_errors = true;
      }
      else
        printf("succeeded.\n");
    }

    return had_errors;
  }

  void add_check(base_check* check) {
    checks.push_back(check);
  }
};



} // end of tst namespace
#endif


