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
        cluster_data = False # specific traeatment
        comp = True
        if os.path.exists(cal_f):
            # some files randomly changes and some modifications must be done
            # before comparing
            ref_f_proc = ref_f + ".proc"
            cal_f_proc = cal_f + ".proc"
            os.system(f" cp {ref_f} {ref_f_proc} ; cp {cal_f} {cal_f_proc} ")
            if ref_f.find("contigs.csv") > -1:
                #remove first column - whole output can different by random
                os.system(f'''cut -d"," -f2- {ref_f} | sort > {ref_f_proc}''')
                os.system(f'''cut -d"," -f2- {cal_f} | sort > {cal_f_proc}''')
            elif ref_f.find("chosen_clusters.csv") > -1:
                #choose a particular column - whole output can be different by random
                os.system(f'''cut -d"," -f5 {ref_f}  > {ref_f_proc}''')
                os.system(f'''cut -d"," -f5 {cal_f}  > {cal_f_proc}''')
                cluster_data = True
                f1 = open(ref_f_proc,"r")
                ls1 = list(map(int,f1.readlines()[1:]))
                f2 = open(cal_f_proc,"r")
                ls2 = list(map(int,f2.readlines()[1:]))
                ls1.sort()
                ls2.sort()
                for i in range(len(ls1)):
                    if abs(ls1[i]-ls2[i]) > 2 :
                        comp = False
            elif ref_f.find("all_clusters.csv") > -1:
                #choose a particular column - whole output can be different by random
                os.system(f'''cut -d"," -f4 {ref_f}  > {ref_f_proc}''')
                os.system(f'''cut -d"," -f4 {cal_f}  > {cal_f_proc}''')
                cluster_data = True
                f1 = open(ref_f_proc,"r")
                ls1 = list(map(int,f1.readlines()[1:]))
                f2 = open(cal_f_proc,"r")
                ls2 = list(map(int,f2.readlines()[1:]))
                ls1.sort()
                ls2.sort()
                for i in range(len(ls1)):
                    if abs(ls1[i]-ls2[i]) > 2 :
                        comp = False
            elif ref_f.find(".tsv"):
                #sort first column - different by random
                os.system(f'''cut -d"," -f2- {ref_f} | sort > {ref_f_proc}''')
                os.system(f'''cut -d"," -f2- {cal_f} | sort  > {cal_f_proc}''')
            if comp:
                comp = filecmp.cmp(ref_f_proc, cal_f_proc)
            if comp:
                print(f"checking against {ref_f}....."+Fore.GREEN+"OK"+Fore.RESET)
            else:
                print(f"checking against {ref_f}....."+Fore.RED+"FAIL"+Fore.RESET)
                different.append(cal_f)
                os.system(f''' echo Reference content: > {cal_f}.diff ''')
                os.system(f''' cat {ref_f_proc} >> {cal_f}.diff ''')
                os.system(f''' echo Calculated content: >> {cal_f}.diff ''')
                os.system(f''' cat {cal_f_proc} >> {cal_f}.diff ''')
                os.system(f''' echo Diff: >> {cal_f}.diff ''')
                os.system(f"diff {cal_f_proc} {ref_f_proc}  >> {cal_f}.diff")
                fail=True
            os.system(f"rm {ref_f_proc} {cal_f_proc} ")
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
parser.add_argument('--testing-result', help='file for testing summary', type=str, default="testing_res.txt")
parser.add_argument('--jenkins-cov', help='Snakemake cov calculations', action='store_true')
args = parser.parse_args()
#Extract the testing files
test1_r1=os.path.abspath("../datasets/testingdata/remote/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R1_001.fastq.gz")
test1_r2=os.path.abspath("../datasets/testingdata/remote/zymo_simulated_real_insert_distribution/Zymo10.1Sim_S1_L001_R2_001.fastq.gz")

if not os.path.exists(test1_r1) or not os.path.exists(test1_r2):
    print("Downloading reference files...")
    os.system(f"cd ../ ; snakemake --use-conda --conda-frontend mamba --configfile testing/testing.yaml -j {args.threads} extract_testing_file")
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

if args.jenkins_cov:
    run_calcualtions=False
    awk = "awk '{print $2}'"
    cmd = f"cd ../ && snakemake extract_testing_file --quiet --use-conda -F -n --conda-frontend mamba --configfile testing/testing.yaml -j 1 | head -n -1 | tail -n +3 | {awk} | sort | uniq >> testing/cov_extract_testing_file.log"
    print(cmd)
    status = os.system(cmd)
    
    if status != 0:
        sys.exit(5)
    cmd = f"cd ../ && snakemake main --quiet --use-conda -F -n --conda-frontend mamba --configfile testing/testing.yaml -j 1 | head -n -1 | tail -n +3 | {awk} | sort | uniq >> testing/cov_main.log"
    print(cmd)
    status = os.system(cmd)
    if status != 0:
        sys.exit(5)

    os.system("cat cov_*.log | sort | uniq > USED_RULES")
    os.system("cd .. && grep -Poh 'rule .+:' $(ls snakefiles/*.smk| grep -vP '/0_') | sed 's/rule //' | sed 's/ //g' | sed 's/://' | sort | uniq > testing/ALL_RULES")
    os.system("grep -xvf USED_RULES ALL_RULES > LEFT_RULES")


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
    os.system(f" cd ../ ;snakemake main --nolock --use-conda --conda-frontend mamba --configfile testing/testing.yaml -j {args.threads} > testing/{run_log} 2>&1")

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
