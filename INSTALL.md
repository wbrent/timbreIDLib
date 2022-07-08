# Using timbreIDLib

timbreIDLib is a library that contains many different Pd objects. If you already have the compiled timbreIDLib binary, you can use it in Pd by either:

  <ol>
    <li> using the [declare] object to load timbreIDLib as needed per patch</li>
    <li> telling Pd to load the timbreIDLib binary file globally at startup via Pd's "Preferences/Startup" dialog</li>
  </ol>

For method 1, assuming timbreIDLib is installed to the default Externals Install Directory, you can simply use the [declare] object in your patch like this:

  > [declare -lib timbreIDLib]

For method 2, make a new startup path in the Startup dialog, and provide the path to the timbreIDLib binary. Note that you must not list the extension of the timbreIDLib file (i.e., pd_linux, pd_darwin, or dll). For example:

> /Users/yourname/Documents/Pd/externals/timbreIDLib

When using Pd's Startup dialog on Windows, note that you can specify the path with forward slashes, and spaces in the path are ok. There is still no need to append the library file extension (dll). For instance:

> C:/Users/Your Name/Documents/Pd/externals/timbreIDLib

Once you have specified the path to timbreIDLib in the Startup dialog, you must quit and restart Pd.

If timbreIDLib is loaded successfully (using either method), you will see a message in Pd's post window stating the timbreIDLib version number.




# Notes on FFTW

As of version 0.7, timbreIDLib uses the FFTW library available at http://www.fftw.org.

FFTW is included as a pre-compiled shared library with timbreIDLib's Windows package available through deken. Simply leave libfftw3f-3.dll in the timbreIDLib directory and everything should work fine. For the Linux and macOS packages on deken, FFTW is statically linked with timbreIDLib.

If you are compiling timbreIDLib from source (see below), you must link it with a single precision compilation of FFTW. To build FFTW in single precision under Linux, configure it like this:

> ./configure CFLAGS="-fPIC" --enable-float

and like this on macOS:

> ./configure CFLAGS="-arch arm64 -arch x86_64" --enable-float

and like this on a Raspberry Pi:

> ./configure CFLAGS="-fPIC" --enable-float

Then run:

> make
>
> sudo make install

On Linux and macOS, the FFTW library files should be installed to /usr/local/lib by default.




# Building timbreIDLib From Source

As of timbreIDLib 0.7.8, pd-lib-builder is used for building:

> https://github.com/pure-data/pd-lib-builder

A (likely outdated) version of pd-lib-builder is already included with the timbreIDLib source, but you can update it via the link above if needed. If you have already compiled and installed FFTW, you can make timbreIDLib from the command line by typing "make" in this directory. On Linux and macOS, timbreIDLib will statically link to the FFTW library.

On Windows, you can use MSYS2 (https://msys2.org) to compile timbreIDLib from the command line. You may need to edit the location and name of the FFTW library in timbreIDLib/Makefile because the libfftw3f-3.dll binary and FFTW header file may not be installed to the expected default locations on your system. Further, Windows expects the FFTW binary name to be "fftw3f-3" rather than the Linux/macOS expectation of "fftw3f." Simply change the ldlibs and cflags lines of timbreIDLib/Makefile to:

> ldlibs = -LC:/my/fftw/library/directory -lfftw3f-3

> cflags = -Iinclude -IC:/my/fftw/library/directory

Finally, after successful compilation, you will either have to set up an environment variable to point to the location of libfftw3f-3.dll, or simply put libfftw3f-3.dll directly in the timbreIDLib directory so that the timbreIDLib.dll binary can access the shared FFTW library.

The FFTW library for Windows is available precompiled at:

> http://www.fftw.org/install/windows.html

You will need the single precision version specifically. Their provided zip file contains several compiled versions of FFTW, but only libfftw3f-3.dll is required for timbreIDLib.
