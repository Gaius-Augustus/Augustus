Prerequisite: Build MySQL
~~~~~~~~~~~~~~~~~~~~~~~~~
    Before you can build MySQL++, you need to build the MySQL client
    library from source.  You're welcome to try linking MySQL++ to the
	native MySQL client libraries instead, but if you have problems with
	it, don't expect any support on the mailing list.
	
	You do not need to build a Cygwin version of the MySQL server.  It's
	still best to use the native Windows version of that.

    Older versions of MySQL do not build cleanly out of the box on
    Cygwin.  This has been fixed recently.  I tested these instructions
    on 5.0.67, with a contemporaneous version of Cygwin 1.5.x.

    You can build it with the standard configure && make && make install
    sequence, but a better configure command for this purpose is:

    $ ./configure --prefix=/usr --sysconfdir=/etc --localstatedir=/var \
        --infodir=/usr/share/info --mandir=/usr/share/man \
        --disable-shared --without-{debug,readline,libedit,server}


Building the Library and Example Programs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    With that done, you can build MySQL++ just as you would any other
    typical Unix program.  See README-Unix.txt for details.
