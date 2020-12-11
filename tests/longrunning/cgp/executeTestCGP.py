#!/usr/bin/env python3

import subprocess
import os
import shutil
import argparse
import sys
import wget
import tarfile
import datetime
from memory_profiler import memory_usage
from concurrent.futures import ThreadPoolExecutor

# import util script from parent directory
sys.path.append('..')
import lr_util as util

# Minimal execution script for longrunning cgp test cases
# based on scripts/executeTestCGP.py by Giovanna Migliorelli.
# The script is executed by GH Actions.

parser = argparse.ArgumentParser(
    description='Execute Augustus long running test cases.')
parser.add_argument('-v', '--eval', action='store_true',
                    help='to evaluate accuracy (respect to the last prediction obtained by launching the script with --run option).')
parser.add_argument('-c', '--chunks',
                    # required=True,
                    nargs='+',
                    help='a list of one or more positive integers indicating the chunk/s to be processed (refer to documentation for a list of chunks over hg38.chr1).')
parser.add_argument('-g', '--augustusDir',
                    help='path to comparative Augustus executable.')
parser.add_argument('-l', '--evalDir',
                    help='path to Eval script.')
parser.add_argument('-d', '--dataDir',
                    help='path to the folder where the test data should be extracted.')
parser.add_argument('-j', '--jobs',
                    help='to set the maximum number of jobs executed in parallel. (default value 2)')
parser.add_argument('-r', '--pathToGitRepo',
                    help='path to the Augustus Git repository.')
args = parser.parse_args()

# used jobs for parallel execution of lr test
jobs = 2

# if not already existing, create dir to collect results for the current chunk
def make_dirs(paths_shared, paths, chunks):
    for chunk in chunks:
        if os.path.exists(paths[chunk]['result_dir']) == False:
            os.makedirs(paths[chunk]['result_dir'])
    if os.path.exists(paths_shared['accuracy']) == False:
        os.makedirs(paths_shared['accuracy'])
    if os.path.exists(paths_shared['joingenes_out_dir']) == False:
        os.makedirs(paths_shared['joingenes_out_dir'])


# create new genomes.tbl for this chunk
def make_genometbl_chunk(paths, chunk):
    cleanup_tbl_chunk(paths, chunk)
    print('Creating new genome tbl for current chunk...',
          paths[chunk]['tbl_test_file'])

    with open(paths[chunk]['tbl_test_file'], 'w') as f:
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            print(species, paths[chunk]['sqlitedb_dir'] +
                  species + '.MINIMAL.fasta', sep='\t', file=f)


def cleanup_db(paths, chunk, removeFASTA=True):
    print('Cleaning up', paths[chunk]['sqlitedb_dir'])

    if removeFASTA:
        if os.path.exists(paths[chunk]['sqlitedb_dir']):
            shutil.rmtree(paths[chunk]['sqlitedb_dir'])
        os.makedirs(paths[chunk]['sqlitedb_dir'])

    print('Cleaning up', paths[chunk]['sqlitedb_test_file'])
    if os.path.exists(paths[chunk]['sqlitedb_test_file']):
        os.remove(paths[chunk]['sqlitedb_test_file'])


def cleanup_tbl_chunk(paths, chunk):
    print('Cleaning up', paths[chunk]['tbl_test_file'])
    if os.path.exists(paths[chunk]['tbl_test_file']):
        os.remove(paths[chunk]['tbl_test_file'])


def cleanup_tbl(paths_shared):
    print('Cleaning up', paths_shared['tbl_file'])
    if os.path.exists(paths_shared['tbl_file']):
        os.remove(paths_shared['tbl_file'])


def port_test(paths_shared, paths, chunks):
    for chunk in chunks:
        print('Porting test for chunk', chunk, '...')
        make_genometbl_chunk(paths, chunk)


# parallel execution : acknowldgement Daniel Honsel (revisited code from test_case.py)
def run_test_parallel(paths_shared, paths, chunks):
    # create a command for each chunk
    args = []
    for chunk in chunks:
        cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'],
               '--alnfile=' + 'nofilegiven',
               '--speciesfilenames=' +
               paths[chunk]['tbl_test_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' +
               paths[chunk]['sqlitedb_test_file'],
               '--optCfgFile=../../../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
               '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' +
               paths_shared['working_dir'] + 'out' +
               str(chunk) + 'run', '--/Testing/testMode=run',
               '--/Testing/workingDir=' + paths_shared['working_dir']]

        output = paths[chunk]['result_dir'] + 'out.runTest'

        args.append([cmd, output, chunk])

    with ThreadPoolExecutor(max_workers=int(jobs)) as executor:
        print('Executing ' + jobs + ' jobs in parallel.')
        for cmd, output, chunk in args:
            print('Adding thread for chunk: ' + str(chunk) + '...')
            executor.submit(execute_test, cmd, output, chunk)


def execute(cmd, output, mode='w'):
    with open(output, mode) as file:
        p = subprocess.Popen(cmd,
                             stdout=file,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)

    p.wait()
    error = p.stderr.read()
    p.stderr.close()
    if error:
        print(error)


def execute_test(cmd, output, chunk, mode='w'):
    print('Runnning prediction on chunk', chunk,
          'using the minimal data set...')
    with open(output, mode) as file:
        p = subprocess.Popen(cmd,
                             stdout=file,
                             stderr=subprocess.PIPE,
                             universal_newlines=True)

    p.wait()
    error = p.stderr.read()
    p.stderr.close()
    if error:
        print(error)


# currently not used (returns accuracy for single chunks) - not parallelized
def run_evaluate(paths, chunk):
    print('Runnning evaluation on chunk', chunk, '...')
    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'],
           paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
    execute(cmd, paths[chunk]['result_dir'] + 'out.eval')


# currently not used (returns accuracy for single chunks) - parallelized
def run_evaluate_parallel(paths, chunks):
    proc_list = []

    for chunk in chunks:
        print('Runnning evaluation on chunk', chunk, '...')
        cmd = [paths_shared['eval_bin'], paths_shared['anno_file'],
               paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
        execute(cmd, paths[chunk]['result_dir'] + 'out.eval')

        filename = paths[chunk]['result_dir'] + 'out.eval'

        with open(filename, 'w') as file:
            proc_list.append(
                subprocess.Popen(cmd,
                                 stdout=file,
                                 stderr=subprocess.PIPE,
                                 universal_newlines=True))
    for p in proc_list:
        p.wait()

    for p in proc_list:
        error = p.stderr.read()
        p.stderr.close()
        if error:
            print(error)


# currently in use (returns accuracy after merging the contributes from all chunks) - parallelized
def run_evaluate_global(paths_shared, paths, chunks):
    print('Runnning evaluation on chunks', chunks, '...')

    with open(paths_shared['joingenes_out_dir'] + 'jgGFFs', 'w') as f:
        for chunk in chunks:
            print(paths_shared['working_dir'] + 'out' +
                  str(chunk) + 'run/hg38.cgp.gff', '1', sep='\t', file=f)

    cmd = [paths_shared['joingenes_bin'], '-f' + paths_shared['joingenes_out_dir'] +
           'jgGFFs', '-o' + paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['joingenes_out_dir'] + 'out.joingenes')

    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'],
           paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['accuracy'] + 'out.eval')


def init_paths_shared(augustusDir, workingDir, evalDir):
    paths_shared = {
        'eval_dir': evalDir,  # '/data/eval',    # path to eval
        # '../',                                    # path to augustus binaries
        'augustus_dir': augustusDir,
        # '../examples/cgp12way/',                        # path to working directory (it contains tree, genome tbl, SQLite db and there results will be written)
        'working_dir': workingDir,
        # path to annotation for hg38.chr1
        'anno_file': workingDir + 'ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.FILTERED.gtf',
    }

    # todo : move here all outXresult
    paths_shared.update({'log': paths_shared['working_dir'] + 'LOG/'})
    # results returned by eval
    paths_shared.update(
        {'accuracy': paths_shared['working_dir'] + 'ACCURACY/'})
    # wrong paths_shared.update({'joingenes_bin' : paths_shared['augustus_dir'] + '/auxprogs/joingenes'})   # path to joingenes exec dir
    # output from joingenes
    paths_shared.update(
        {'joingenes_out_dir': paths_shared['working_dir'] + 'JOINGENES/'})
    # path to SQlite for full length genomes (only if a new data set is to be built)
    paths_shared.update(
        {'sqlitedb_file': paths_shared['working_dir'] + 'SQLITE/12way.db'})
    # path to genomes.tbl for length genomes (only if a new data set is to be built)
    paths_shared.update(
        {'tbl_file': paths_shared['working_dir'] + 'GENOMETBL/genomes.tbl'})
    # path to tree
    paths_shared.update(
        {'tree_file': paths_shared['working_dir'] + 'TREE/ucsc12way.nwk'})
    # path to augustus exec
    paths_shared.update(
        {'augustus_bin': paths_shared['augustus_dir'] + 'bin/augustus'})
    # path to joingenes exec
    paths_shared.update(
        {'joingenes_bin': paths_shared['augustus_dir'] + 'auxprogs/joingenes/joingenes'})
    # path to eval
    paths_shared.update(
        {'eval_bin': paths_shared['eval_dir'] + 'evaluate_gtf.pl'})

    return paths_shared


def init_paths(chunks):
    paths = {}
    for chunk in chunks:
        dictionary = {}

        # the following paths automatically reflect changes in paths_shared
        # will contain minimal fasta after their extraction
        dictionary.update(
            {'sqlitedb_dir': paths_shared['working_dir'] + 'minimalFasta' + str(chunk) + '/'})
        dictionary.update(
            {'sqlitedb_test_file': paths_shared['working_dir'] + 'SQLITE/12wayTEST_' + str(chunk) + '.db'})
        dictionary.update(
            {'tbl_test_file': paths_shared['working_dir'] + 'GENOMETBL/genomesTEST_' + str(chunk) + '.tbl'})
        dictionary.update(
            {'result_dir': paths_shared['working_dir'] + 'out' + str(chunk) + 'result/'})

        paths[chunk] = dictionary

    return paths


def expand_dir(path):
    tmp = path
    if len(path) > 0 and path[len(path)-1] != '/':
        tmp += '/'
    return tmp


# Get prepared data for the test case.
def get_test_data(dataDir):
    if os.path.exists(dataDir):
        shutil.rmtree(dataDir)
    os.mkdir(dataDir)

    url = 'http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/cgp12way.tgz'
    filename = os.path.join(dataDir, 'cgp12way.tar.gz')
    wget.download(url, out=filename)
    if (os.path.isfile(filename)):
        tar = tarfile.open(filename)
        tar.extractall(dataDir)
        tar.close()
        os.remove(filename)


def execute_lr_test(paths_shared, paths, chunks):
    start = datetime.datetime.now()
    mem_usage = memory_usage(
        (run_test_parallel, (paths_shared, paths, chunks)), interval=60, include_children=True)
    max_mem_usage = max(mem_usage) / 1000.0
    end = datetime.datetime.now()
    execution_time = (end - start).total_seconds() / 60.0
    return {'default' : {'execution_time' : execution_time, 'used_memory' : max_mem_usage}}


if __name__ == '__main__':
    if args.chunks is None:
        print('No chunks specified, please make use of --chunks to pass a non empty list of positive integers...')
        sys.exit()

    chunks = [int(x) for x in list(dict.fromkeys(args.chunks))]
    # range valid for chr1 chunk size 2.5 Mb, chunk overlap 0.5 Mb
    chunks = [x for x in chunks if x > 0 and x < 126]

    if len(chunks) == 0:
        print('No valid chunks specified...')
        sys.exit()

    if args.augustusDir is None:
        print('Path to comparative augustus executable required, please make use of --augustusDir to pass the path...')
        sys.exit()
    augustusDir = str(expand_dir(args.augustusDir))

    if args.dataDir is None:
        print('Path to test data folder required, please make use of --dataDir to pass the path...')
        sys.exit()
    dataDir = str(expand_dir(args.dataDir))

    # set working directory according to test data stored on the webserver
    workingDir = dataDir + 'cgp12way/'

    evalDir = ''
    if args.eval:
        if args.eval is None:
            print(
                'Path to Eval script required, please make use of --evalDir to pass the path...')
            sys.exit()
        evalDir = str(expand_dir(args.evalDir))

    get_test_data(dataDir)

    paths_shared = init_paths_shared(augustusDir, workingDir, evalDir)
    paths = init_paths(chunks)

    make_dirs(paths_shared, paths, chunks)

    port_test(paths_shared, paths, chunks)

    if args.jobs:
        jobs = args.jobs

    util.check_memory()

    res = execute_lr_test(paths_shared, paths, chunks)

    if args.eval:
        run_evaluate_global(paths_shared, paths, chunks)

    # collect commit information for database storage
    if args.pathToGitRepo is None:
        info = 'NoInformation', 'NoInformation'
    else:
        info = util.commit_info(args.pathToGitRepo)

    util.store_additional_data(
        info[1], info[0], res, workingDir + 'additional_information.json')
