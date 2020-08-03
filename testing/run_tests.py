#!/usr/bin/env python
import argparse, os, colorama, filecmp
from colorama import Fore
colorama.init()
def compare_with_reference(test, calc):
    """
        Compare files in two dirrectories

        Parameters:
            test (str):  directory with calculated files
            calc (str):  directory with reference (expected) files
        Returns:
            not fail: fail/pass of the comparision
            missing: list of  missing files
            different: list of files with differences
    """
    listOfFiles = list()
    for (dirpath, dirnames, filenames) in os.walk(test):
            listOfFiles += [os.path.join(dirpath, file) for file in filenames]
    fail = False
    missing = []
    different = []
    for fr in listOfFiles:
        ref_f = fr
        cal_f = fr.replace(test,calc)
        #print(cal_f)
        fail_ok = True
        if os.path.exists(cal_f):
            comp = filecmp.cmp(ref_f, cal_f)
            if comp:
                print(f"checking against {ref_f}....."+Fore.GREEN+"OK"+Fore.RESET)
            else:
                print(f"checking against {ref_f}....."+Fore.RED+"FAIL"+Fore.RESET)
                different.append(cal_f)
                os.system(f"diff {cal_f} {ref_f}  > {cal_f}.diff")
                fail=True
        else:
            fail = True
            print(Fore.RED+f"Missing expected file: {cal_f}"+Fore.RESET )
            missing.append(cal_f)

    return (not fail,missing,different)

parser = argparse.ArgumentParser(description='Run integration tests')
parser.add_argument('--only-compare',
                    help='Do not clean tetsing run tmp and out directories.\
                    Only compare against a reference files', action='store_true')
parser.add_argument('--threads', help='Number of threads. Currently testing is done only on one node', type=int, default=4)
parser.add_argument('--singularity-mount-point', help='Please indicate singularity mount points separated via comma withoutspaces', type=str, default="")
parser.add_argument('--testing-result', help='file for testing summary', type=str, default="testing_res.txt")
args = parser.parse_args()
#check singularity mount point
sing_path = ""
sing_path = args.singularity_mount_point
if args.singularity_mount_point  == "":
    if not os.environ["SINGULARITY_BINDPATH"] == "":
        sing_path = os.environ["SINGULARITY_BINDPATH"]
        print("Will use Singularity binding path:", os.environ["SINGULARITY_BINDPATH"])
    else:
        raise Exception("Any directory fo singularity binding is not defined.\n Pease define an enviromental variable SINGULARITY_BINDPATH or provide it via variable --singularity-mount-point  ")
else:
    os.environ["SINGULARITY_BINDPATH"]=args.singularity_mount_point
    print(f'Seting the enviromental variable SINGULARITY_BINDPATH to {os.environ["SINGULARITY_BINDPATH"]}')
#Extract the testing files
test1_r1=os.path.abspath("../datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R1_001.fastq.gz")
test1_r2=os.path.abspath("../datasets/testingdata/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R2_001.fastq.gz")
if not os.path.exists(test1_r1) or not os.path.exists(test1_r2):
    print("Extractiong reference files...")
    os.system(f"cd ../ ; snakemake --use-conda --conda-frontend mamba --configfile testing/testing.yaml -j {args.threads}  extract_testing_file ")
print(f"For the test set #1 these input files will be used:\n\t{test1_r1}\n\t{test1_r2}")
test1_refs = "testingdata/reference"
if not os.path.exists(test1_refs):
    print("Extractiong reference set")
    os.system(f"tar xvzf {test1_refs}.tar.gz")
else:
    print(f"For testing the files from the folder {test1_refs} will be used")


#first test seet of files more ... be added later
calc1="testingdata/calculated"
run_calcualtions=True
if args.only_compare:
    if os.path.exists(calc1):
        print(f"Analysing dir {calc1} wo calculatioins")
        run_calcualtions=False
    else:
        print(f"Directory {calc1} doesn't exist. Option --only-comapre is ignored")
if run_calcualtions:
    print(f"Cleaning the directory {calc1}")
    os.system(f"rm -r {calc1} ")
    print("Starting testing run")
    run_log = calc1+"._run_log.txt"
    print(f"\t...log could be checked here: {run_log}")
    os.system(f" cd ../ ;  snakemake --nolock  --use-conda --conda-frontend mamba --configfile testing/testing.yaml -j {args.threads} >testing/{run_log} 2>&1  ")

print("TESTING GENERATED FILES WITH EXPECTED. SET#1")
ok,missing,different=compare_with_reference(test1_refs,calc1)
print(".............................analysis completed..")
f_res = open(args.testing_result,"w")
if ok:
    print(Fore.GREEN + "ALL TESTS PASSED" + Fore.RESET )
    print(Fore.GREEN + "ALL TESTS PASSED" + Fore.RESET,file=f_res )
else:
    print(Fore.RED + "SOME TESTS FAILED" + Fore.RESET )
    print(Fore.RED + "SOME TESTS FAILED" + Fore.RESET,file=f_res )
    if len(missing) > 0:
        print("Missing:")
        print("Missing:", file=f_res)
        for g in missing:
            print(f"\t{g}")
            print(f"\t{g}",file=f_res)
    if len(different):
        print("Different:")
        print("Different:",file=f_res)
        for g in different:
            print(f"\t check {g}.diff")
            print(f"\t check {g}.diff",file=f_res)
print(f"Summary written to: {args.testing_result}")
f_res.close()



