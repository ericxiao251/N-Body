// Wrap printf with MACRO LOG
#ifndef __LOG_H__
#define __LOG_H__

#if defined ENABLE_LOG && ENABLE_LOG > 0
  #define LOG(a) printf a
#else
  #define LOG(a) (void)0
#endif

#endif