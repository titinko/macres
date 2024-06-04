# MacRes - A cross-platform resampler
**[Download this resampler](https://github.com/titinko/macres/releases)**

(Fork of https://github.com/ohac/tn_fnds)

MacRes is short for Mac Resampler

This is an executable designed to be called by https://github.com/titinko/utsu (UTSU).
Despite the name, it can be compiled for Windows, Mac, or Linux using the Makefile.

The Windows/Linux versions of this resampler can also be used with the UTAU program:
https://en.wikipedia.org/wiki/Utau.

## Building from source

This program uses two external libraries: [libpyin](https://github.com/Sleepwalking/libpyin) and [libgvps](https://github.com/Sleepwalking/libgvps).

You will need to git clone both repositories, install their libraries using make install, and update the LIB_PREFIX field in macres's Makefile to match the directory where the .a and .h files were generated. You will also need to manually edit pyin.h to wrap its function definitions with "extern C" as described in [this bug](https://github.com/Sleepwalking/libpyin/issues/10).
