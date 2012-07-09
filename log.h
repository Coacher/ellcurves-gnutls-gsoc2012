/* logging facility */
#ifndef _LOG_H
#define _LOG_H

/* global debug level */
extern int debug;

/* timing stuff */
#ifdef ENABLE_DEBUG
#include <time.h>

clock_t past, present, future;
clock_t secs[16];
clock_t next_secs[16];

#define MKTIMESTAMP(ARG) ARG=clock()

#define TIMEBETWEEN(PREV, NEXT) ((double) (NEXT - PREV ) / CLOCKS_PER_SEC)

#else
#define MKTIMESTAMP(ARG)
#define TIMEBETWEEN(PREV, NEXT) 0
#endif

/* output and other stuff */
#ifdef ENABLE_DEBUG
#define DEBUG_CALL(EXP) \
do {            \
    if (debug) { \
        EXP;    \
    }           \
} while(0)
#else 
#define DEBUG_CALL(EXP)
#endif

#ifdef ENABLE_DEBUG
#define DEBUG_CALL_L(LEVEL, EXP) \
do {            \
    if (debug >= LEVEL) { \
        EXP;    \
    }           \
} while(0)
#else 
#define DEBUG_CALL_L(LEVEL, EXP)
#endif

/* prints to stderr if debug is set to non-zero */
void dbg_msg(const char *fmt, ...);
/* prints to stderr if debug is gt level */
void dbg_msg_l(int level, const char *fmt, ...);

/* gmp analogs */
/* prints to stderr if debug is set to non-zero */
void gmp_dbg_msg(const char *fmt, ...);
/* prints to stderr if debug is gt level */
void gmp_dbg_msg_l(int level, const char *fmt, ...);
#endif
