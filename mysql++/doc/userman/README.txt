The user manual is written in XML DocBook format, version 4.2.
It uses the official DocBook XSL stylesheets, and will build with
versions as old as 1.61.2.  (Why these versions?  They're what comes
with CentOS 3, the oldest system I still use.)

To make the HTML version of the user manual, just type 'make' in this
directory.  (Or 'make html' if you want to be pedantic.)  To make the
PDF version of the manual, say 'make pdf'.  To make both versions,
say 'make all'.



The most likely cause of build failures is not having the right
processing programs installed.  The DocBook processing model looks
like this:

    DocBook --[XSLT proc.]-+-> HTML
                   ^       |
                   |       +-> XSL-FO --[XSL-FO proc.]--> PDF
                   |
           {XSL stylesheets}

"DocBook" above is a file in DocBook XML format.  In this case,
it's in userman.dbx.

There are many alternatives for the tools in the square brackets:

    The first replaceable piece is an XSLT processor, which
    translates XML into other text formats, such as HTML or other
    varieties of XML.  We use xsltproc from the Gnome project.
    (http://xmlsoft.org/XSLT/)  There are numerous alternatives,
    but each supports different extensions to the standard, so it's
    simplest if everyone uses the same processor for a given document.

    We use the XSLT processor to do two transforms.  One is directly to
    HTML.  The other is to XSL-FO, an XML-based page layout language.
    We do this intermediary transform because XSLT is not good at
    creating binary formats like PDF.

    The second replaceable piece in the diagram above is an XSL-FO
    processor, which converts XSL-FO to a more directly useful page
    layout format, like PDF.  The user manual's build system supports
    several alternatives.

    The build system relies on a simple script in this directory --
    fo2pdf -- to find an XSL-FO formatter and run it.  It looks first
    for RenderX's XEP, which comes in a free-as-in-beer version for
    personal use.  (http://renderx.com/download/personal.html)  If
    you're in a commercial environment, RenderX wants you to use their
    commercial trial version which will format this manual without
    complaint, but it puts watermarks and blank pages into the output.
    They want $300 for the single-user to get clean output.  It's the
    same as the free personal version, just with a different license.
    You don't need the higher-end versions of XEP; they don't do
    anything we need here.

    If fo2pdf can't find XEP, it then looks for Antenna House's XSL
    Formatter (http://antennahouse.com/).  It's pretty much the same
    deal as XEP: crippled demo version for testing, and a single-user
    version for $300.  There is no free version for personal use,
    however.

    Failing all that, fo2pdf falls back to the only free-as-in-speech
    XSL-FO formmatter, Apache FOP (http://xmlgraphics.apache.org/fop/).
    FOP isn't yet available through most Unixy package download
    systems, so you'll have to download it directly from the source.
    The Debian repositories do have it, so on any Debian based distro
    (e.g. Ubuntu) you can just say "sudo apt-get install fop".

    You might be wondering why fo2pdf looks for FOP last, given that
    MySQL++ is itself free software and relies on a lot of other
    free software.  It's just that it's a good bet that if there's
    a commercial processor on the system, it was put there quite
    purposefully by someone who paid (or caused money to be paid) for
    it, and so wants it to be used.  The commercial vendors can still
    get money for their products because FOP hasn't caught up with
    them in several important areas.  That said, don't feel that you
    need to go and buy an XSL-FO processor just to build the manuals.
    We try to always keep the manual in a state where FOP can generate
    adequate output.


The third replaceable piece above is the DocBook XSL stylesheet set.
The stylesheets are the XSLT processor's rules, controlling how
the input XML gets transformed to the output format.  The standard
DocBook stylesheet set (link below) includes stylesheets for HTML and
XSL-FO output.  By default, xsltproc looks for these first on your
local system, and failing to find them, tries to download them on
the fly from the Internet.  Because this slows processing quite a bit
even if you have a fast Internet connection, and it obviously doesn't
work when your net connection is down, we've disabled this option.
Therefore, you must have the DocBook XSL stylesheets installed to
build the user manual.

Most Unixy type systems have pre-built DocBook XSL stylesheet packages
available:

    Red Hat/Fedora: docbook-style-xsl RPM package
    Mac OS X:       docbook-xsl Fink package (http://fink.sf.net)
    Cygwin:         docbook-xml?? package (?? = DocBook version)
    Ubuntu/Debian:  docbook-xsl, from the standard APT repository

(Please send the name of the package for your system to the mailing
list if it isn't listed above, and I'll add it to the list.)

If you can't find a package for your system, you can get the DocBook
stylesheets from the source: http://docbook.sourceforge.net/  They're
a bit tricky to set up correctly, so it's better to use a pre-built
package if you can.

If you are still having problems, post the details about it to the
MySQL++ mailing list, and I'll try to help you debug the problem.
You might also find the FOP and/or DocBook mailing lists helpful.




If you're looking to hack on the manual, here are some helpful
resources for getting up to speed on DocBook:

    Mills' "Installing And Using An XML/SGML DocBook Editing Suite"
    article:

        http://tinyurl.com/8alb2

        This is the best tutorial I've found.
        

    Walsh and Muellner's _DocBook: The Definitive Guide_ book, second
    edition, online version:

        http://www.docbook.org/tdg/en/html/docbook.html

        This is the official DocBook referece.


    DocBook FAQ:

        http://www.dpawson.co.uk/docbook/

        Go here when you have a question that the tutorials and
        references do not answer.


    The official DocBook site:

        http://docbook.org/

