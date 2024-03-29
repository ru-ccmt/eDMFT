// @Copyright 2007 Kristjan Haule
// 
#ifndef ASSERT_
#define ASSERT_

// #ifdef NO_ARG_CHECK
// #define Assert(condition, message)
// #else /* NO_ARG_CHECK */
// #include <iostream>
// #define Assert(condition, message)\
// {\
//   if(!(condition)) std::cerr << (message) << std::endl;\
// }
// #endif /* NO_ARG_CHECK */

#ifdef _DEBUG
#define LOG(x) x
#else
#define LOG(x)
#endif

#ifdef _DEBUG
#define CHECK(x) x
#else
#define CHECK(x)
#endif

#endif /* ASSERT_ */
