# AUGUSTUS Test Cases

The test cases located in the folder `tests/short` are used to check whether AUGUSTUS is executable and optionally whether AUGUSTUS returns exactly the expected results.

## Requirements
To execute the test cases, the following requirements must be fulfilled.

1. The execution script must be called from `tests/short` as working directory.
2. The AUGUSTUS binaries must be accessible in `../../bin` and the directory structure should be the same as in the GitHub repository. The last point is important to open required test data.
3. Python version 3.6 or higher must be installed.
4. If MySQL should also be tested, then the corresponding Python module is needed. It can be installed as follows.
        
        pip3 install mysql-connector-python

## Execute Test Cases
If all requirements are fulfilled, the test cases can be executed from the directory `tests/short` as follows.

    ./execute_test.py [options] <testcase>

### Required Argument for the Test Case Module to be Executed
Semantically related test cases are grouped in modules. The argument `<testcase>` specifies the module whose tests are to be executed. There are currently the following different test modules available.

- **examples**: This module contains test cases that test the main program AUGUSTUS. The test cases are based on sample data located in the examples folder. The test cases can be found in the module `examples/test_examples.py`.
- **bam2wig**: This module contains test cases that test the program `bam2wig` and can be found in `auxprogs/bam2wig/test_bam2wig.py`.
- **homGeneMapping**: Test cases for the program `homGeneMapping`. The module can be found in `auxprogs/homgenemapping/test_homgenemapping.py`. To execute this tests further requirements have to be fulfilled as described [here](../../auxprogs/homGeneMapping/README.TXT).

### Optional Arguments
The following optional arguments can be used.

- **--compare**: After the test run, the generated results are compared with reference data. If both are not identical, an error is thrown and the diff is output.
The reference data is stored for each module in a folder named `<path/to/module/>expected_results`.
- **--html**: Should be used in combination with the `--compare` option. If a test fails, the diff is additionally stored in the folder `html_output` as HTML file.
- **--mysql**: CGP test cases of the `examples` module and so far also only these are also executed using a MySQL database. If this option is used, a MySQL database must be set up and the configuration file `examples/testconfig.json` must be adjusted accordingly.
- **--clean**: Remove all files created during a test run. If this option is set, no tests are executed.

### GitHub Actions Workflow (Test and Build)
The test cases mentioned here are also used for CI tests in the AUGUSTUS GitHub repository. GH-Actions are used for this and currently the test modules `examples` and `bam2wig` are called from the workflow *Test and Build*. For example, the call for the `examples` module looks like this.

    ./execute_test.py --compare --html --mysql examples

This call means that the usage of a MySQL database should also be tested, that the results should be compared with reference data, and that in case of an failed test run the diffs should also be saved as HTML.

## Create a New Test Case
There are two ways to add a new test case. Firstly, you can simply add it to an existing module. Secondly, you can create a new test module and add the test case there.

### Add a New Test Case to an Existing Module
If the new test case to be created semantically matches an existing module, it can be added there. After that, the expected reference results still have to be added to the `expected_results` folder of the module.

### Add a New Module
A new module should be in its own package. To integrate it into the existing execution structure, the module should implement the methods `execute(compare, html, mysql)` and `clean()`. Afterwards, it can be added to the list of possible test case modules (`testcases`) in `execute_test.py`.