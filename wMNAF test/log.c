/* logging facility */
#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <gmp.h>

#include "log.h"

/* write message to stderr if debug >= 0 */
void dbg_msg(const char *fmt, ...)
{
#ifdef ENABLE_DEBUG
	if (debug) {
        va_list args;
	
        va_start(args, fmt);
    	vfprintf(stderr, fmt, args);
    	va_end(args);
    }
#endif
}

/* write message to stderr if debug >= level */
void dbg_msg_l(int level, const char *fmt, ...)
{
#ifdef ENABLE_DEBUG    
	if (debug >= level) {
        va_list args;
	
        va_start(args, fmt);
    	vfprintf(stderr, fmt, args);
    	va_end(args);
    }
#endif
}

/* gmp analogs */
/* write message to stderr if debug >= 0 */
void gmp_dbg_msg(const char *fmt, ...)
{
#ifdef ENABLE_DEBUG
	if (debug) {
        va_list args;
	
        va_start(args, fmt);
    	gmp_vfprintf(stderr, fmt, args);
    	va_end(args);
    }
#endif
}

/* write message to stderr if debug >= level */
void gmp_dbg_msg_l(int level, const char *fmt, ...)
{
#ifdef ENABLE_DEBUG    
	if (debug >= level) {
        va_list args;
	
        va_start(args, fmt);
    	gmp_vfprintf(stderr, fmt, args);
    	va_end(args);
    }
#endif
}
