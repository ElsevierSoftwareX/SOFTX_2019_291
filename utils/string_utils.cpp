#include "mt_utils.h"

#ifndef WINBUILD
#include <glob.h>      // for glob()
#endif

#include <libgen.h>    // for basename()

bool parse_param (const std::string & param, const std::string & flag,
                  std::string & res, const mt_check_type check)
{
  bool matched = false;
  if(param.compare(0, flag.size(), flag) == 0)
  {
    res.assign(param.begin()+flag.size(), param.end());
    matched = true;
  }

  if(matched) {
    switch(check) {
      case chk_nonzero:
        check_nonzero(res.size(), "parse_param"); break;
      case chk_fexists:
        check_file_exists(res, "parse_param"); break;
      case dont_check: break;
    }
  }

  return matched;
}

void fixBasename(std::string & base)
{
  size_t len = base.size();
  if( (len > 0) && (base[len-1] == '.') ) base.resize(len-1);
}

#ifndef WINBUILD
void mt_glob(const std::string& pattern, mt_vector<std::string> & ret)
{
  glob_t glob_result;
  //glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
  glob(pattern.c_str(), 0, NULL, &glob_result);

  // Initialise return vector
  ret.resize(glob_result.gl_pathc);

  // Put results in return vector
  for(unsigned int i=0; i<glob_result.gl_pathc; i++)
    ret[i] = std::string(glob_result.gl_pathv[i]);
  globfree(&glob_result);
}
#endif

bool endswith(const std::string value, const std::string ending) {
  if (ending.size() > value.size())
    return false;
  return std::equal(ending.rbegin(), ending.rend(), value.rbegin());
}

std::string mt_basename(const std::string path)
{
  char * cstr = new char[path.length()+1];
  std::strcpy(cstr, path.c_str());
  std::string ret = std::string(basename(cstr));
  delete [] cstr;

  return ret;
}

char* mt_basename(const char* path, mt_vector<char> & strbuff)
{
  if(strbuff.size() < strlen(path))
    strbuff.resize(strlen(path)+1);

  std::strcpy(strbuff.data(), path);
  char* ret = basename(strbuff.data());
  return ret;
}

void split_string(const std::string & input, char s, mt_vector<std::string> & list)
{
  std::size_t fnd = 0, ofnd = 0;
  int num_parts = 1, idx = 0;

  fnd = input.find(s);
  while (fnd != std::string::npos)
  {
    num_parts++;
    ofnd = fnd+1;
    fnd = input.find(s, ofnd);
  }

  list.resize(num_parts);

  fnd = input.find(s); ofnd = 0;
  while (fnd != std::string::npos)
  {
    list[idx++].assign(input.begin() + ofnd, input.begin() + fnd);
    ofnd = fnd+1;
    fnd = input.find(s, ofnd);
  }
  if(ofnd < input.size())
    list[idx].assign(input.begin() + ofnd, input.end());
}

bool find_extension(std::string fname, std::string ext, size_t & pos)
{
  pos = fname.rfind(ext);
  return (pos != std::string::npos) && (fname.size() == (pos + ext.size()));
}

bool remove_extension(std::string & fname)
{
  size_t pos = fname.rfind(".");
  if(pos != std::string::npos && pos != 0) {
    fname.resize(pos);
    return true;
  }

  return false;
}

int swap_char(std::string & str, const char from, const char to)
{
  int cnt = 0;

  for(size_t i=0; i<str.size(); i++)
    if(str[i] == from) {
      str[i] = to;
      cnt++;
    }

  return cnt;
}


