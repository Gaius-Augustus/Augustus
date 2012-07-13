Prerequisite: GCC Version
~~~~~~~~~~~~~~~~~~~~~~~~~
    If your MinGW version isn't using at least GCC 3.4.5, it needs
    to be updated.  Older versions are known to not work with MySQL++.


Prerequisite: MySQL C Development Files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MySQL++ is built atop MySQL's C API library, so you need to have
    MySQL installed on your development system to get the current C API
    development files.

    If you do a default installation of MySQL, the development files
    probably won't be installed. As of this writing you have to do
    either a Complete or Custom install to get these files.  (They keep
    changing the way the Windows installer works, so this may not be
    true any more by the time you read this.)

    The MySQL++ Makefile assumes that you installed MySQL in
    
        C:\Program Files\MySQL\MySQL Server 5.0\
    
    If not, you have two options.

    The simplest is to edit Makefile.mingw.  This is a generated file,
    but if that's the only change to MySQL++ you need, it works fine.

    If you're doing deeper work on MySQL++, you should change the
    variable MYSQL_WIN_DIR at the top of mysql++.bkl instead.  Then to
    generate Makefile.mingw from that file, you will need the Win32
    port of Bakefile from http://bakefile.org/  The command to do
    that is:

        bakefile_gen -f mingw


Prerequisite: MySQL C API DLL Import Library
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Before you can build MySQL++ with MinGW, you will need to create
    a MinGW-compatible import library for MySQL's C API library.
    Using the current default install path for MySQL and assuming
    MySQL++ is in c:\mysql++, the commands to do this are:

        cd C:\Program Files\MySQL\MySQL Server 5.0\lib\opt
        dlltool -k -d c:\mysql++\libmysqlclient.def -l libmysqlclient.a


Building the Library and Example Programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    With the prerequisites above taken care of, you can build MySQL++
    with this command:

        mingw32-make -f Makefile.mingw

    Notice that we're using the MinGW-specific version of GNU make, not
    the Cygwin or MSYS versions.  Many things will break otherwise: path
    separator handling, shell commands used by the Makefile, etc.

    Speaking of Cygwin and MSYS, if you have either these or any other
    Unix emulation environment installed, be sure their executables
    aren't in the PATH when building MySQL++.  MinGW's version of GNU
    make does some funny things if it thinks it's running in the
    presence of Unixy tools, which will break the MySQL++ build.

    Once the library is built, you should run the examples.  At minimum,
    run resetdb and simple1.

    Once you're satisfied that the library is working correctly, you can
    run install.hta to automatically install the library files and
    headers in subdirectories under c:\mysql++.


Cygwin and MinGW Coexistence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    It's possible to have both Cygwin and MinGW installed and build
    with the MinGW tools without interference from the Cygwin bits.
    The main thing you have to take care of is that MinGW's bin
    directory must precede the Cygwin bin directory in the PATH,
    so that its tools are found first.  If you use Cygwin's bash
    as a command shell in preference to the DOS-like cmd.exe, you
    can use this shell script to temporarily set the environment to
    "MinGW mode" and make it easy to get back to "Cygwin mode":

        #!/bin/sh
        PATH=/c/mingw/bin:/c/windows:/c/windows/system32:/c/cygwin/bin
        echo "Say 'exit' to leave MinGW shell and restore Cygwin environment."
        /usr/bin/bash --rcfile ~/.mingwrc

    I recommend having at least this in the ~/.mingwrc file:

        alias make=mingw32-make
        PS1='MinGW: \W \$ '

    The prompt change reminds you that you are in a sub-shell set up for
    MinGW.  The alias for 'make' ensures you don't accidentally run
    Cygwin's make, which won't work with Makefile.mingw.  We could just
    leave /c/cygwin/bin out of the environment, but there are Cygwin
    tools we want access to, like vim.  As long as all the MinGW ones
    override those Cygwin also provides, we don't need to worry about
    having both in the PATH.  Besides, having the alias is nice for
    those who have 'make' committed to muscle memory.


Building on Linux
~~~~~~~~~~~~~~~~~
    You might wish to build MySQL++ with MinGW because you're
    not actually running Windows, but need Windows executables.
    The thought being that this lets you use GCC, the same compiler
    you're probably using to make native executables.  There are
    indeed ways to make this work.

    The most "native" way to do this is to run MinGW under Wine.
    Leonti Bielski provided these instructions:

        1. Install MinGW through Wine:

           $ wine MinGW-5.1.6.exe

        2. Add the MinGW directory to Wine's PATH with Wine regedit:

           http://winehq.org/site/docs/wineusr-guide/environment-variables

        3. Install MySQL under Wine, or at least unpack the Windows
           ZIP file version of MySQL in a place where Wine can find it.
           You don't need to run a Windows MySQL server under Wine.
           We're only doing this to get the MySQL C API library and
           its headers, which MySQL++ builds against.  The resulting
           MinGW build of MySQL++ can talk to a native MySQL server
           out in Wine's host environment or on some other machine.

        4. Modify Makefile.mingw to match the install location for
           the MySQL C API files.

        5. Create libmysqlclient.a as described above, except with
           minor differences for running under Wine:

           $ wine mingw32-dlltool -k -d /native/path/libmysqlclient.def...

        6. Build MySQL++ with:
        
           $ wine mingw32-make -f Makefile.mingw

    Another way is to build a Windows virtual machine, such as with
    VMware or VirtualBox.  In that case, you'd use the regular build
    instructions at the top of this document.

    You might think to avoid the need for Wine or Windows by use of a
    MinGW cross-compiler:

        $ ./configure --target=mingw32
        $ make

    Unfortunately, that currently doesn't work.

    The reason is that our autoconf build system assumes a
    typical POSIX type target, which MinGW is not.  We made this
    assumption because we have a perfectly good MinGW build option,
    Makefile.mingw.  But, that also won't work on a POSIX system
    because that Makefile assumes external commands run under cmd.exe,
    not some Unixy shell.  Thus the advice to build with Makefile.mingw
    under Windows or something sufficiently close to it.

    If you really wanted to, you could extend the autoconf build system
    to make it realize when it's being used to cross-compile for MinGW.
    Patches thoughtfully considered; see HACKERS.txt.
