#ifndef STRING_UTILS_H
#define STRING_UTILS_H

#include <string>

/**
* @brief Parse a meshtool-style parameter into a result string.
*
* @param param   The parameter.
* @param flag    The flag we want to match.
* @param res     The result string containing the value of the parameter.
* @param check   The optional check to perform on the parameter.
*
* @return        True if the parameter contained the flag.
*/
bool parse_param (const std::string & param, const std::string & flag,
                  std::string & res, const mt_check_type check = dont_check);

/// globing
void mt_glob(const std::string& pattern, mt_vector<std::string> & ret);

/**
* @brief Check if a string ends with a certain substring.
*
* @param value    String we check.
* @param ending   Substring we check with.
*
* @return True if there is a match.
*/
bool endswith(const std::string value, const std::string ending);

/**
* @brief Wrapper for glibc basename.
*
* @param path  Input path.
*
* @return      Output basename.
*/
std::string mt_basename(const std::string path);
char*       mt_basename(const char* path, mt_vector<char> & strbuff);

/**
* @brief Add a '.' to the basename of mesh if it is not present.
*
* @param base The mesh basename.
*/
void fixBasename(std::string & base);

/**
* @brief Split a string in a list of substrings.
*
* Note that the delimiter is not stored.
*
* @param input  The input string.
* @param s      The substring delimiter.
* @param list   The list of substrings.
*/
void split_string(const std::string & input, char s, mt_vector<std::string> & list);

/**
* @brief Check if filename contains the given extension
*
* @param fname The filename
* @param ext   The extension
* @param pos   Position where the extension was found
*
* @return Whether the extension was found
*/
bool find_extension(std::string fname, std::string ext, size_t & pos);
bool remove_extension(std::string & fname);

/// swap a char to another one
int swap_char(std::string & str, const char from, const char to);
#endif
