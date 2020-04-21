/**
* @file progress.hpp
* @brief CLI progress class.
* @author Aurel Neic
* @version
* @date 2017-02-11
*/

#ifndef _MT_PROGRESS_H
#define _MT_PROGRESS_H


#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#ifndef WINBUILD
#include <sys/ioctl.h>
#endif
#include <unistd.h>

#include <string>

/**
* @brief Base class for tracking progress.
*
* @tparam V Integer type
*/
template<class V>
class progress
{
  protected:
  V _num_next;       ///< number of times next() will be triggered
  V _counter;        ///< number of times next() has been triggered
  std::string _msg;  ///< message that will be displayed before progress information.

  /// evaluate progress and display it if necessary, this needs to be implemented
  /// by each non-abstract class individually
  virtual void evaluate() = 0;

  public:
  progress(V num_next, const char* msg) :
           _num_next(num_next),
           _msg(msg)
  {
    _counter = 0;
  }

  /// increment progress counter by one
  void next()
  {
    ++_counter;
    this->evaluate();
  }

  /// increment progress counter by a custom value
  void increment(V inc)
  {
    _counter += inc;
    this->evaluate();
  }

  /// finish up progress display
  void finish() {
    _counter = _num_next;
    this->evaluate();
    printf("\n");
  }
};

/**
* @brief Do not display progress
*
* @tparam V Integer type.
*/
template<class V>
class progress_silent : public progress<V>
{
  private:
  /// evaluate progress, this needs to be implemented
  /// by each non-abstract class individually
  void evaluate() {}

  public:
  progress_silent(V num_next, const char* msg) :
           progress<V>(num_next, msg)
  {
    printf("%s ...\n", msg);
    fflush(stdout);
  }
};

/**
* @brief Display progress as a percentage
*
* @tparam V Integer type.
*/
template<class V>
class progress_percent : public progress<V>
{
  private:
  V _threshold;        ///< counter threshold for displaying
  std::string _blank;  ///< blank space before the progress is displayed

  /**
  * @brief Display the progress message
  *
  * @param percentage The displayed percentage
  */
  void display(V percentage)
  {
    const std::string & msg = progress<V>::_msg;

    printf("\r%s%s%d %c", msg.c_str(), _blank.c_str(), int(percentage), '%');
    fflush(stdout);
  }

  /**
  * @brief Evaluate if a new display is needed
  */
  void evaluate()
  {
    const V & num_next = progress<V>::_num_next;
    const V & counter  = progress<V>::_counter;
    V prog = V(float(counter) / float(num_next) * 100.0f);

    if(prog > _threshold)
    {
      _threshold = prog;
      this->display(prog);
    }
  }

  public:
  progress_percent(V num_next, const char* msg, V start = 60) :
                      progress<V>(num_next, msg),
                      _threshold(0)
  {
    int size = start - progress<V>::_msg.size();
    if(size < 0) size = 0;

    _blank.assign(size, ' ');
    this->display(0);
  }
};


/**
* @brief Display progress as a bar
*
* @tparam V Integer type
*/
template<class V>
class progress_bar : public progress<V>
{
  private:
  V _threshold;   ///< threshold for displaying

  const V _bar_len;       ///< total amount of characters representing progress
  const char _bar_char;   ///< character used to fill the bar

  std::string _bar;       ///< bar
  std::string _blank_a;   ///< blank space before bar starts
  std::string _blank_b;   ///< blank space not yet filled by bar

  /**
  * @brief Display progress
  *
  * @param nbars Number of characters representing the current bar.
  */
  void display(V nbars)
  {
    const std::string & msg = progress<V>::_msg;
    if(nbars > _bar_len) nbars = _bar_len;

    _bar.assign(nbars, _bar_char);
    _blank_b.assign(_bar_len - nbars, ' ');

    printf("\r%s%s[%s%s]", msg.c_str(), _blank_a.c_str(), _bar.c_str(), _blank_b.c_str());
    fflush(stdout);
  }

  /**
  * @brief Evaluate if a new display is needed
  */
  void evaluate()
  {
    const V & num_next = progress<V>::_num_next;
    const V & counter  = progress<V>::_counter;
    V prog = V(float(counter) / float(num_next) * float(_bar_len));

    if(prog > _threshold)
    {
      _threshold = prog;
      this->display(prog);
    }
  }

  public:
  /**
  * @brief Constructor
  *
  * @param num_next     Estimate of how many times next() will be issued.
  * @param msg          Message to display together with progress bar.
  * @param bar_start    Number of characters after which the bar should start.
  * @param bar_len      Number of characters making up the bar.
  * @param bar_char     Character used to represent the bar.
  */
  progress_bar(V num_next, const char* msg, V bar_start = 60, V bar_len = 20, char bar_char = '=') :
                  progress<V>(num_next, msg),
                  _threshold(0),
                  _bar_len(bar_len),
                  _bar_char(bar_char)
  {
    int size = bar_start - progress<V>::_msg.size();
    if(size < 0) size = 0;

    _blank_a.assign(size, ' ');
    this->display(0);
  }
};

/**
* @brief Display progress as a bar with added ETA
*
* @tparam V Integer type
*/
template<class V>
class progress_bar_eta : public progress<V>
{
  private:
  V _threshold;   ///< threshold for displaying

  const V _bar_len;       ///< total amount of characters representing progress
  const char _bar_char;   ///< character used to fill the bar

  std::string _bar;       ///< bar
  std::string _blank_a;   ///< blank space before bar starts
  std::string _blank_b;   ///< blank space not yet filled by bar
  std::string _blank_c;   ///< blank space covering ETA in case we are finished
  float _perc;            ///< progress percentage, used for ETA estimation
  struct timeval _tstart; ///< starting time, used for ETA estimation
  struct timeval _tnow;   ///< current time, used for ETA estimation

  char* _dynbuff;

  /**
  * @brief Display progress
  *
  * @param nbars Number of characters representing the current bar.
  */
  void display(V nbars)
  {
    int columns = 0;

    #ifndef WINBUILD
    if (isatty(fileno(stdin))) {
      struct winsize w;
      ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
      columns = w.ws_col;
    }
    #endif

    const std::string & msg = progress<V>::_msg;
    if(nbars > _bar_len) nbars = _bar_len;

    _bar.assign(nbars, _bar_char);
    _blank_b.assign(_bar_len - nbars, ' ');

    gettimeofday(&_tnow, NULL);
    float needed_time, tdiff;
    int nsec = 0, nmin = 0, nhour = 0;

    if(nbars > 0 && nbars < _bar_len) {
      tdiff = _tnow.tv_sec - _tstart.tv_sec;
      needed_time = (tdiff / _perc) * (1.0f - _perc);
      nhour = int(needed_time) / 3600;
      nmin  = (int(needed_time) % 3600) / 60;
      nsec  = int(needed_time) % 60;

      sprintf(_dynbuff, "\r%s%s[%s%s] %02dh:%02dm:%02ds ETA ",
             msg.c_str(), _blank_a.c_str(), _bar.c_str(), _blank_b.c_str(), nhour, nmin, nsec);
    }
    else {
      sprintf(_dynbuff, "\r%s%s[%s%s] ",
      msg.c_str(), _blank_a.c_str(), _bar.c_str(), _blank_b.c_str());
    }

    // we check if the message fits the terminal
    int len = strlen(_dynbuff);
    if(columns && len > columns) {
      int diff = len - columns;
      int oldsize = _blank_a.size();
      int newsize = diff < oldsize ? oldsize - diff : 0;

      // first we reduce the blanks
      _blank_a.assign(newsize, ' ');

      if(nbars > 0 && nbars < _bar_len)
        sprintf(_dynbuff, "\r%s%s[%s%s] %02dh:%02dm:%02ds ETA ",
            msg.c_str(), _blank_a.c_str(), _bar.c_str(), _blank_b.c_str(), nhour, nmin, nsec);
      else
        sprintf(_dynbuff, "\r%s%s[%s%s] ",
            msg.c_str(), _blank_a.c_str(), _bar.c_str(), _blank_b.c_str());

      _blank_a.assign(oldsize, ' ');
      len = strlen(_dynbuff);
      diff = len - columns;

      // if the message still doesnt fit, we trim if from the left
      if(diff > 0) {
        _dynbuff[0] = '\r';
        int start = diff+1, ridx, widx;
        for(ridx=start, widx=1; ridx < len; ridx++, widx++)
          _dynbuff[widx] = _dynbuff[ridx];

        _dynbuff[widx] = '\0';
      }
    }

    // output final message
    if(columns && (nbars == 0 || nbars == _bar_len)) {
      _blank_c.assign(columns-1, ' ');
      printf("\r%s", _blank_c.c_str());
    }
    printf("%s", _dynbuff);
    fflush(stdout);
  }

  /**
  * @brief Evaluate if a new display is needed
  */
  void evaluate()
  {
    const V & num_next = progress<V>::_num_next;
    const V & counter  = progress<V>::_counter;
    V prog = V(float(counter) / float(num_next) * float(_bar_len));
    _perc  = float(counter) / float(num_next);

    if(prog > _threshold)
    {
      _threshold = prog;
      this->display(prog);
    }
  }

  public:
  /**
  * @brief Constructor
  *
  * @param num_next     Estimate of how many times next() will be issued.
  * @param msg          Message to display together with progress bar.
  * @param bar_start    Number of characters after which the bar should start.
  * @param bar_len      Number of characters making up the bar.
  * @param bar_char     Character used to represent the bar.
  */
  progress_bar_eta(V num_next, const char* msg, V bar_start = 55, V bar_len = 12, char bar_char = '=') :
                  progress<V>(num_next, msg),
                  _threshold(0),
                  _bar_len(bar_len),
                  _bar_char(bar_char)
  {
    int size = bar_start - progress<V>::_msg.size();
    _dynbuff = new char[1024];

    if(size < 0) size = 0;
    gettimeofday(&_tstart, NULL);

    _blank_a.assign(size, ' ');
    this->display(0);
  }

  virtual ~progress_bar_eta() {
    delete [] _dynbuff;
  }
};




#endif

