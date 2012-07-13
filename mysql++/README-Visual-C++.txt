Visual C++ Version Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    MySQL++ is fully functional with Visual C++ 2005 and 2008.

    MySQL++ also works with VC++ 2003 (a.k.a. VC++ 7.1), with the
    exception of the SSQLS feature.  There was partial support for
    SSQLS with VC++ 2003 in MySQL++ v2, but a feature we added in
    MySQL++ v3.0 crashes the VC++ 2003 compiler when you try to use
    even simple SSQLS, so we had to remove support for this entirely
    for that platform.  (See the v3.0 section in the Breakages chapter
    of the user manual for details.)

    Older versions of Visual C++ are basically hopeless when it
    comes to building current versions of MySQL++.  They have too
    many weaknesses in their Standard C++ implementation to build a
    modern library like MySQL++.  If you cannot upgrade your compiler,
    my advice is that you're best off programming straight to the
    MySQL C API rather than try to make MySQL++ build.


Where Are the Project Files, and Why So Many Versions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    There are three sets of .sln and .vcproj files shipped with
    MySQL++, in the vc2003, vc2005 and vc2008 subdirectories.
    Other than the SSQLS issue brought up above, there no functional
    difference between these versions.  We ship separate project
    files for each version of Visual Studio partly to save you from
    having to walk through the project conversion wizard, and partly
    so you can build the library with multiple versions of Visual C++
    without conflicts among the object files.


Prerequisites
~~~~~~~~~~~~~
    You need to have the Windows version of the MySQL server installed
    on your development system, even if you always access a MySQL
    server on a different machine.  This is because in addition to
    installing the server itself, the official MySQL Windows binaries
    also install the client-side development files that MySQL++
    needs in order to communicate with a MySQL server.  We call this
    the MySQL C API; MySQL++ is a C++ wrapper for this C API.

    You have to do a Custom install to enable installation of these
    development files.  If you get an error about mysql-version.h or
    mysql.h when building MySQL++, go back and reinstall the MySQL
    server, paying careful attention to the options.

    If you've installed the development files and are still getting
    build errors, read on.


Using Nonstandard MySQL Installations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The Visual Studio project files that come with MySQL++ have
    everything set up correctly for the common case.  The biggest
    assumption in the settings is that you're using MySQL 5.0
    (the newest version considered suitable for production use
    at the time of this writing) and that you installed it in
    the default location:

        C:\Program Files\MySQL\MySQL Server 5.0\

    If you installed a different version, or it's in a different
    directory, you need to change the project file settings to
    reference the C API development files in that other location.
    There are two ways to do this.

    The hard way is to make 4 different changes to 39 separate
    project files.  If you're a talented Visual Studio driver, you
    can do this in as little as about 5 or 6 steps.  You might even
    get it right the first time.

    The easy way is to install Bakefile (http://bakefile.org/),
    change the value of the MYSQL_WIN_DIR variable near the top of
    mysql++.bkl in the top level of the MySQL++ source tree, and run
    rebake.bat.  This will rebuild all of the project files for you,
    using the new MySQL path in all the many places it's needed.


Building MySQL++ for 64-Bit
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The MySQL++ Visual Studio project files ship with the assumption
    that you're building for 32-bit Windows.  While the utility
    of running the MySQL server on a 64-bit system is clear, that
    doesn't mandate running its client programs in 64-bit mode, too.
    As a result, we haven't yet bothered to come up with an easy way
    to change this.  Here's the hard way:

    - Follow the steps above to change the MySQL install location,
      if necessary.

    - Open the solution file, then say Build > Configuration Manager

    - Pull down the "Active solution platform" list box, say New...,
      and add "x64", copying settings from the Win32 build, and
      allowing it to create new project platforms.

    - Pull the same list box down, say Edit..., and remove the Win32
      build, unless you actually need both versions.

    It should then build.  If you get PRJ2009 errors, complaining
    about BuildLog.htm, this is a bug in Visual Studio 2008's
    new parallel build feature, which seems to affect x64 builds
    more than x86.  Sometimes you can get around it by just saying
    "build" repeatedly, each time getting a few more modules built.
    Or, you can disable the feature by going to Tools > Options... >
    Projects and Solutions > Build and Run > # of parallel processes
    and setting the value to 1.

    If you don't follow the instructions above carefully, you can
    end up with a completely unbuildable solution.  (It happened
    to me once when preparing these instructions!)  If you simply
    cannot seem to make it work, it's often simplest to just blow
    away the MySQL++ source tree, unpack a fresh copy and try again,
    paying more attention to the details.


Building the Library and Example Programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    You must build both the Debug and Release versions of the library,
    because a release build of your program won't work with a Debug
    version of the MySQL++ DLL.  These DLLs get different names, so
    you can install them in the same directory if needed: mysqlpp_d.dll
    for the Debug version, and mysqlpp.dll for the Release version.

    With the library built, run at least the resetdb and simple1
    examples to ensure that the library is working correctly.
    In addition to the other generic examples, there are a few
    Visual C++ specific examples that you might want to look at in
    examples\vstudio.  See README-examples.txt for further details.

    Once you're sure the library is working correctly, you can run
    the install.hta file at the project root to install the library
    files and headers in a directory of your choosing.
    
    (Aside: You may not have come across the .hta extension before.
    It's for a rarely-used feature of Microsoft's Internet Explorer,
    called HTML Applications.  Know what Adobe AIR is?  Kinda like
    that, only without the compilation into a single binary blob which
    you must install before you can run it.  Just open install.hta
    in a text editor to see how it works.)


Using MySQL++ in Your Own Projects
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    This is covered in the user manual, chapter 9.


Working With Bakefile
~~~~~~~~~~~~~~~~~~~~~
    MySQL++'s top-level Visual Studio project files aren't
    maintained directly.  Instead, we use a tool called Bakefile
    (http://bakefile.org/) to generate them from mysql++.bkl. Since
    there are so many project files in MySQL++, it's often simpler to
    edit this source file and "re-bake" the project files from it than
    to make your changes in Visual Studio.

    To do this, download the native Windows version of Bakefile from the
    web site given above.  Install it, and then put the installation
    directory in your Windows PATH.  Then, open up a command window, cd
    into the MySQL++ directory, and type "rebake".  This will run
    rebake.bat, which rebuilds the Visual Studio project files from
    mysql++.bkl.

    There's more information about using Bakefile in HACKERS.txt.


If You Run Into Problems...
~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Especially if you have linking problems, make sure your project
    settings match the above.  Visual C++ is very picky about things
    like run time library settings.  When in doubt, try running one
    of the example programs.  If it works, the problem is likely in
    your project settings, not in MySQL++.

