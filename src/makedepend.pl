#!/usr/bin/perl
# @(#) makedepend  List make target dependencies
# Author: Dale Mensch, July 1993
# Modified by: Becca Thomas, January 1994

# CONFIGURATION PARAMETERS:

# Using GNU make?
$have_gmake = 0;                        # 0 = No, 1 = Yes

# Choose the source-file suffix:
$source_suffix = '.cc';                  # Possibilities: .C, .cpp, .cc

# Choose the object-file suffix:
$object_suffix = '.o';                  # Another possibility: .xo

# List directories containing header files that don't change:
$ignore_include_dirs = '/usr';  # Don't need to list these

# Location of standard header files:
$add_include_dirs = '-I../include';   # If not included automatically

# Command that runs preprocessor, includes headers, excludes comments:
$cpp_cmd = 'cc -E';

# MAIN PROGRAM:

# Check correct command line usage:
die "Usage: $0 target [target...]\n" unless (@ARGV);# At least one

# Open pipe to read make (or gmake) output:
if ($have_gmake) {
    open(MAKE, "gmake -n -W $ARGV[0] |");       # -W outdates target
} else {
    # Rename target because if non-existent it's treated as out-of-date:
    ($sample_target = $ARGV[0]) =~ s/$source_suffix$/$object_suffix/;
    if (-e $sample_target) {    # If the sample target exists
        $saved_target = "$sample_target"."~";   # Backup copy name
        rename($sample_target, $saved_target);  # Squirrel away
    }
    open(MAKE, "make -n $sample_target |");     # Make treats as dated
}

# Determine necessary include directives:
while (<MAKE>) {    # For each make-command output line
    chop;           # Remove trailing newline
    @cmd_line = split(/[ \t]+/, $_);            # Tokenize at whitespace

    while (@cmd_line) {
        $make_token = shift (@cmd_line);        # For each token
        if ($make_token =~ /^-I|^-D/) {         # Include/define?
            $directive_token = $make_token;     # Yes, save it
            if (! $is_directive{$directive_token}) {    # Not seen
                $is_directive{$directive_token} .= $directive_token;
                # Add to catenated string of include directories:
                $cpp_directives .= " ";         # (Catenation saves
                $cpp_directives .= $directive_token;# original order)
            }
        }
    }   # End of while ($make_token = ...
}       # End of while (<MAKE>) {

# Add any necessary include directories:
$cpp_directives .= " $add_include_dirs ";

# For each source file named on command line
while (@ARGV) {                         # While there are arguments
    undef %includes;                    # Reset array before next file
    $source_file = shift;               # Save name for later use
    if (! -r $source_file) {            # Can file be read?
        print STDERR "Source file \"$source_file\" not readable.\n";
        next;   # Try another source file, if available
    }
    # Preprocessor names all include directives:
    open(CPP, "$cpp_cmd $cpp_directives $source_file |");

    # Create list of all included files:
    while (<CPP>) {
        s#//.*$##;                      # Delete C++ "//" comment
        if (/^(#|#line)\s*\d+\s+["<](\S+)[">]\s+.*$/ &&
          ! m.$ignore_include_dirs. && ! m.built-in. ) { # Grab relevant includes
                $includes{$2} .= $2 unless $includes{$2};   # Save name
        }   # End of if (/^ ...
    }       # End of while (<CPP>) {

    # Report target name:
    $_ = $source_file;
    s/$source_suffix$/$object_suffix/;  # Create target name
    print "$_ :";                       # Dependency or rules line

    # List the prerequisites for the target:
    foreach $dependency (sort keys(%includes)) {
        print " \\\n";  # continuation line
        print "\t$includes{$dependency}";# Important: tab first
    }
    print "\n\n";                       # Empty lines between targets
}
rename($saved_target, $sample_target) if $saved_target;
