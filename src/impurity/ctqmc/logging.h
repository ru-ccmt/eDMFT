// @Copyright 2007 Kristjan Haule
// 
#ifndef __LOGGING__
#define __LOGGING__

#include <iostream>

using namespace std;


enum { TRACE, STATS, TIMING };

#ifdef _LOGGING

map<int,ostream*> slog;
#define LOGV(ilog,name) do { (*slog[ilog])<<__FILE__<<" ["<<__LINE__<<"] : "<<#name<<" = "<<(name)<< std::endl;} while(false)
#define LOG(ilog,name) do { (*slog[ilog])<<__FILE__<<" ["<<__LINE__<<"] : "<<name<< std::endl;} while(false)
#define LOGN(ilog,name) do { (*slog[ilog])<<__FILE__<<" ["<<__LINE__<<"] : "<<name;} while(false)
#define LOGN0(ilog,name) do { (*slog[ilog])<<" "<<name;} while(false)

#else

#define LOGV(ilog,name) do{}while(false)
#define LOG(ilog,name) do{}while(false)
#define LOGN(ilog,name) do{}while(false)
#define LOGN0(ilog,name) do{}while(false)

#endif

/*
struct nullstream : ostream {
    nullstream() : ios(0), ostream(0) {}
};

class SLog {
    map<int, ostream*> _streams;
    map<int, bool> _enabled;
    nullstream nout;

public:
    void add_log(int ilog, ostream& stream, bool enabled = true) {
	_streams[ilog] = &stream;
	_enabled[ilog] = enabled;
    }

    void enable(int ilog) { _enabled[ilog] = true; }
    void disable(int ilog) { _enabled[ilog] = false; }

    const bool enabled(int ilog) { return _enabled[ilog]; }

    ostream& operator() (int ilog) { return _enabled[ilog] ? *_streams[ilog] : nout; }
};

#ifdef _LOG
class DoLog {
public:
  DoLog();
  ~DoLog();
  void AddStream( std::ostream& );
  void AddStream( const std::string& );
  void Write( const std::string& );
};
typedef DoLog Log;
#else
class NoLog {
public:
  inline NoLog() {}
  inline ~NoLog() {}
  inline void AddStream( std::ostream& ) {}
  inline void AddStream( const std::string& ) {}
  inline void Write( const std::string& ) {}
};
typedef NoLog Log;
#endif
*/

#endif /* __LOGGING__ */
