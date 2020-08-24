#!/usr/bin/env python3
import subprocess
import os
import shutil
import argparse
import sys
import wget
import tarfile
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


# Minimal execution script for longrunning cgp test cases
# based on scripts/executeTestCGP.py by Giovanna Migliorelli.
# The script is executed by GH Actions. 

parser = argparse.ArgumentParser(description='Execute Augustus long running test cases.')
parser.add_argument('-v', '--eval', action='store_true',
                    help='to evaluate accuracy (respect to the last prediction obtained by launching the script with --run option).')
parser.add_argument('-c', '--chunks',
                    #required=True,
                    nargs='+',
                    help='a list of one or more positive integers indicating the chunk/s to be processed (refer to documentation for a list of chunks over hg38.chr1).') 
parser.add_argument('-g', '--augustusDir',
                    help='path to comparative Augustus executable.')  
parser.add_argument('-l', '--evalDir',
                    help='path to Eval script.')  
parser.add_argument('-w', '--workingDir',
                    help='path to data set used in testing (link).')
parser.add_argument('-j', '--jobs',
                    help='to set the maximum number of jobs executed in parallel. (default value 2)')                    
args = parser.parse_args()


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
    print('Creating new genome tbl for current chunk...', paths[chunk]['tbl_test_file'])

    with open(paths[chunk]['tbl_test_file'], 'w') as f:
        for species in ['hg38', 'rheMac3', 'mm10', 'rn6', 'oryCun2', 'bosTau8', 'canFam3', 'loxAfr3', 'echTel2', 'dasNov3', 'monDom5', 'galGal4']:
            print(species, paths[chunk]['sqlitedb_dir'] + species + '.MINIMAL.fasta', sep='\t', file=f)

    
def cleanup_db(paths, chunk, removeFASTA = True):
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
def run_test_parallel(paths_shared, paths, chunks, jobs=2):

    proc_list = []

    # create a command for each chunk
    args = []
    for chunk in chunks:
        cmd = [paths_shared['augustus_bin'], '--species=human', '--treefile=' + paths_shared['tree_file'], 
        '--alnfile=' + 'nofilegiven',
        '--speciesfilenames=' + paths[chunk]['tbl_test_file'], '--softmasking=1', '--alternatives-from-evidence=0', '--dbaccess=' + paths[chunk]['sqlitedb_test_file'],
        '--optCfgFile=../../../config/cgp/cgp_param_21features_accuracy_largerGrid.convDivCorrect.cfg',
        '--stopCodonExcludedFromCDS=true', '--/CompPred/logreg=on', '--/CompPred/outdir=' + paths_shared['working_dir'] + 'out' + str(chunk) + 'run', '--/Testing/testMode=run',
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
    print('Runnning prediction on chunk', chunk, 'using the minimal data set...')    
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
    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
    execute(cmd, paths[chunk]['result_dir'] + 'out.eval')


# currently not used (returns accuracy for single chunks) - parallelized
def run_evaluate_parallel(paths, chunks):
    proc_list = []

    for chunk in chunks:
        print('Runnning evaluation on chunk', chunk, '...')
        cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff']
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


# currently in use (returns accuracy after merging the contributes from all chunks) - parallelized
def run_evaluate_global(paths_shared, paths, chunks):
    print('Runnning evaluation on chunks', chunks, '...')
    
    with open(paths_shared['joingenes_out_dir'] + 'jgGFFs', 'w') as f:
        for chunk in chunks:
            print(paths_shared['working_dir'] + 'out' + str(chunk) + 'run/hg38.cgp.gff', '1', sep='\t', file=f)  
        
    cmd = [paths_shared['joingenes_bin'], '-f' + paths_shared['joingenes_out_dir'] + 'jgGFFs', '-o' + paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['joingenes_out_dir'] + 'out.joingenes')

    cmd = [paths_shared['eval_bin'], paths_shared['anno_file'], paths_shared['joingenes_out_dir'] + 'joingenes.gff']
    execute(cmd, paths_shared['accuracy'] + 'out.eval')


def init_paths_shared(augustusDir, workingDir, evalDir):
    paths_shared = {
    'eval_dir' : evalDir,           #'/data/eval',    # path to eval
    'augustus_dir' : augustusDir,   # '../',                                    # path to augustus binaries
    'working_dir' : workingDir,     # '../examples/cgp12way/',                        # path to working directory (it contains tree, genome tbl, SQLite db and there results will be written)
    'anno_file' : workingDir + 'ENSEMBL/ensembl.ensembl_and_ensembl_havana.chr1.CDS.gtf.dupClean.FILTERED.gtf',  # path to annotation for hg38.chr1
    }

    paths_shared.update({'log' : paths_shared['working_dir'] + 'LOG/'})                             # todo : move here all outXresult
    paths_shared.update({'accuracy' : paths_shared['working_dir'] + 'ACCURACY/'})                   # results returned by eval
    # wrong paths_shared.update({'joingenes_bin' : paths_shared['augustus_dir'] + '/auxprogs/joingenes'})   # path to joingenes exec dir
    paths_shared.update({'joingenes_out_dir' : paths_shared['working_dir'] + 'JOINGENES/'})         # output from joingenes
    paths_shared.update({'sqlitedb_file' : paths_shared['working_dir'] + 'SQLITE/12way.db'})        # path to SQlite for full length genomes (only if a new data set is to be built)
    paths_shared.update({'tbl_file' : paths_shared['working_dir'] + 'GENOMETBL/genomes.tbl'})       # path to genomes.tbl for length genomes (only if a new data set is to be built)
    paths_shared.update({'tree_file' : paths_shared['working_dir'] + 'TREE/ucsc12way.nwk'})         # path to tree
    paths_shared.update({'augustus_bin' : paths_shared['augustus_dir'] + 'bin/augustus'})           # path to augustus exec
    paths_shared.update({'joingenes_bin' : paths_shared['augustus_dir'] + 'auxprogs/joingenes/joingenes'})  # path to joingenes exec
    paths_shared.update({'eval_bin' : paths_shared['eval_dir'] + 'evaluate_gtf.pl'})                # path to eval

    return paths_shared  


def init_paths(chunks):
    paths = {}
    for chunk in chunks:
        dictionary = {}

        # the following paths automatically reflect changes in paths_shared
        dictionary.update({'sqlitedb_dir' : paths_shared['working_dir'] + 'minimalFasta' + str(chunk) + '/'})               # will contain minimal fasta after their extraction
        dictionary.update({'sqlitedb_test_file' : paths_shared['working_dir'] + 'SQLITE/12wayTEST_' + str(chunk) + '.db'})
        dictionary.update({'tbl_test_file' : paths_shared['working_dir'] + 'GENOMETBL/genomesTEST_' + str(chunk) + '.tbl'})
        dictionary.update({'result_dir' : paths_shared['working_dir'] + 'out' + str(chunk) + 'result/'})

        paths[chunk] = dictionary

    return paths


def expand_dir(path):
    tmp = path
    if len(path)>0 and path[len(path)-1] != '/':
        tmp += '/'
    return tmp


# Get prepared data for the test case. 
def get_test_data(dataDir):
    if os.path.exists(dataDir):
        shutil.rmtree(dataDir)

    url = 'http://bioinf.uni-greifswald.de/bioinf/downloads/data/aug-test/cgp12way.tgz'
    filename = os.path.join(Path(dataDir).parent, 'cgp12way.tar.gz')
    wget.download(url, out=filename)
    if (os.path.isfile(filename)):
        tar = tarfile.open(filename)
        tar.extractall(Path(dataDir).parent)
        tar.close()
        os.remove(filename)


if __name__ == '__main__':
    if args.chunks is None:
        print('No chunks specified, please make use of --chunks to pass a non empty list of positive integers...')
        sys.exit()

    chunks = [int(x) for x in list(dict.fromkeys(args.chunks))]
    chunks = [x for x in chunks if x>0 and x<126]          # range valid for chr1 chunk size 2.5 Mb, chunk overlap 0.5 Mb

    if len(chunks) == 0:
        print('No valid chunks specified...')
        sys.exit()

    if args.augustusDir is None:
        print('Path to comparative augustus executable required, please make use of --augustusDir to pass the path...')
        sys.exit()
    augustusDir = str(expand_dir(args.augustusDir))  

    if args.workingDir is None:
        print('Path to data set used in testing required, please make use of --workingDir to pass the path...')
        sys.exit()
    workingDir = str(expand_dir(args.workingDir))

    evalDir = ''
    if args.eval:
        if args.eval is None:
            print('Path to Eval script required, please make use of --evalDir to pass the path...')
            sys.exit()
        evalDir = str(expand_dir(args.evalDir))

    get_test_data(workingDir)

    paths_shared = init_paths_shared(augustusDir, workingDir, evalDir)
    paths = init_paths(chunks)

    make_dirs(paths_shared, paths, chunks)

    port_test(paths_shared, paths, chunks)
    if args.jobs:
        run_test_parallel(paths_shared, paths, chunks, jobs=args.jobs)
    else:
        run_test_parallel(paths_shared, paths, chunks)

    if args.eval:
        run_evaluate_global(paths_shared, paths, chunks)
