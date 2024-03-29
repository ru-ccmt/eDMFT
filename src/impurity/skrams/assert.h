// @Copyright 2007 Kristjan Haule
// 
#ifndef ASSERT_
#define ASSERT_

#ifdef NO_ARG_CHECK
#define Assert(condition, message)
#else /* NO_ARG_CHECK */
#include <iostream>
#define Assert(condition, message)\
{\
  if(!(condition)) std::cerr << (message) << std::endl;\
}
#endif /* NO_ARG_CHECK */

#ifdef CTMA_DEBUG
int debug_flag1, debug_flag2;
#define CTMA_LOG(x) x
#else
#define CTMA_LOG(x)
#endif

#endif /* ASSERT_ */
