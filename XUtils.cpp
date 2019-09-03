#include "XUtils.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <csignal>

void string_copy(char *dest, const char* src, size_t n) {
#ifdef _WIN32
	strcpy_s(dest, n, src);
#else
	strncpy(dest, src, n);
#endif

}

void string_concat(char *dest, const char* src, size_t n) {
#ifdef _WIN32
	strcat_s(dest, n, src);
#else
	strncat(dest, src, n);
#endif

}

void string_print(char *str, size_t size, const char *format, ...) {
	va_list args;
	va_start(args, format);
#ifdef _WIN32
	vsprintf_s(str, size, format, args);
#else
	vsnprintf(str, size, format, args);
#endif
	va_end(args);
}

FILE* file_open(const char *path, const char *mode) {
#ifdef _WIN32
	FILE* f;
	errno_t fileOpenResult = fopen_s(&f, path, mode);
	return f;
#else
	return fopen(path, mode);
#endif
}