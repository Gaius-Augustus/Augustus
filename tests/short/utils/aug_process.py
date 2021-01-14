#!/usr/bin/env python3

import os
import subprocess
import unittest


def execute(testcase, cmd, out=subprocess.PIPE):
    """ Execute the specified cmd and writes the stdout into a file or return
        as string.
        Throws an assertion exception if the return code of the executed cmd was
        not zero or if stderr is not empty.
        Throws an assertion exception if stdout should be written to a file
        but was not.
    :param testcase: the calling unittest.TestCase or None
    :param cmd: cmd string executed in a shell
    :param out: a file name or None
    :return: an empty string if stdout was writte to a file or the stdout as string
    """
    isFile = isinstance(out, str)
    output = open(out, 'w') if isFile else out

    p = subprocess.Popen(cmd,
                         stdout=output,
                         stderr=subprocess.PIPE,
                         shell=True,
                         universal_newlines=True)
    rc = p.wait()
    error = p.stderr.read()
    p.stderr.close()

    testcase = unittest.TestCase() if testcase is None else testcase
    if error:
        print("error " + error)
    # testcase.assertEqual(error, '', error)
    testcase.assertEqual(rc, 0, f'Return code not 0! Error: {error}' )

    if isFile:
        testcase.assertTrue(os.path.isfile(out),
                            'Output file was not created as expected!')
    else:
        stdout = p.stdout.read()
        p.stdout.close()
        return stdout

    return ''
