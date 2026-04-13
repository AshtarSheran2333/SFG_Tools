#if !defined(__UTILS_ERROR_MACROS__)
#define __UTILS_ERROR_MACROS__

#define error_io_check(ierr, msg) call utils_error_io_check(ierr, msg, __FILE__, __LINE__)
#define error_stop(msg) call utils_error_stop(msg, __FILE__, __LINE__)

#endif
