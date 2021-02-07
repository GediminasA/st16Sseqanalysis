//
def tests_started = false
def COUNT_PASS = ''
def COUNT_FAIL = ''
def COUNT_TOTAL = ''
def COUNT_MAIN_FAIL = ''
def DATE_TAG = java.time.LocalDateTime.now()
def MAIN_DIR = "/mnt/jenkins/st16SseqAnalysis"
echo "RUNNING WITH PARAMS: ${params.WITH_PARAMS}"
def BRANCH = ""

if(params.WITH_PARAMS){
    BRANCH = params.BRANCH
}
else{
    BRANCH = env.GIT_BRANCH
}
BRANCH = BRANCH.replaceAll(/[^\w]/, "_")

pipeline {
    agent {
        node {
            label 'master'
            customWorkspace "${MAIN_DIR}/st16SseqAnalysis_$BRANCH"
        }
    }
    environment {
        MINI_PATH="${MAIN_DIR}/MINIs/$BRANCH"
        JULIA_DIR="${MAIN_DIR}/JULIA_st16SseqAnalysis"
        PATH_JULIA="$JULIA_DIR/julia-1.5.3/bin"
        ENV_PATH="${MAIN_DIR}/ENVs/TEST_ENV_$BRANCH"
        ENV_PATH_J="${MAIN_DIR}/ENVs/TEST_JAVA11_$BRANCH"
        PKGS_DIR="/mnt/jenkins/CONDA_PKGS"
        PATH="$MINI_PATH/bin:$PATH_JULIA:$PATH"
        JULIA_PKGDIR="/mnt/jenkins/JULIA_PKGS"
        JULIA_DEPOT_PATH="$JULIA_DIR/.julia"
    }
    stages {
        stage('Install CONDA') {
            steps {
                sh """#!/usr/bin/env bash
                    wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
                    bash Miniconda3-latest-Linux-x86_64.sh -b -p $MINI_PATH
                    conda config --add channels conda-forge
                    conda config --add channels bioconda
                    conda config --add channels anaconda
                    conda config --add channels r
                    conda config --add pkgs_dirs $PKGS_DIR
                """
            }
            post {
                success {
                    sh 'rm Miniconda3-latest-Linux-x86_64.sh'
                }
            }
        }

        stage('Create Envs') {
            steps {
                parallel (
                    "Install Julia":{
                        sh """#!/usr/bin/env bash
                            mkdir $JULIA_DIR
                            wget https://julialang-s3.julialang.org/bin/linux/x64/1.5/julia-1.5.3-linux-x86_64.tar.gz
                            tar -xvzf julia-1.5.3-linux-x86_64.tar.gz -C $JULIA_DIR
                            rm julia-1.5.3-linux-x86_64.tar.gz
                            cd scripts/julia_modules/st16SseqJuliaTools
                            ./setup.sh
                        """
                    },
                    "Create Java ENV": {
                        sh '''#!/usr/bin/env bash
                            source deactivate
                            conda create -p $ENV_PATH_J --yes
                            source activate $ENV_PATH_J
                            conda install -c conda-forge openjdk=11.0.1
                            source deactivate
                        '''
                    },
                    "Create Workflow ENV": {
                        sh '''#!/usr/bin/env bash
                            conda install -c conda-forge mamba
                            mamba env create -p $ENV_PATH -f envs/main.yaml
                        '''
                    }
                )
            }
        }

        stage('Coverage') {
            steps {
                parallel (
                    "Snakemake Coverage": {
                        script{
                            tests_started = true
                            sh """#!/usr/bin/env bash
                                source activate $ENV_PATH
                                cd testing
                                ./run_tests.py --jenkins-cov
                            """
                            sh """#!/usr/bin/env bash
                                source activate $ENV_PATH
                                python testing/generate_cov.py --inp \$(ls snakefiles/*.smk| grep -vP "/0_" | tr "\n" " ") --out smk-cov.xml --uncov testing/LEFT_RULES
                            """
                        }
                    }
                )
            }
        }

        stage('Workflow') {
            steps {
                parallel (
                    "Workflow tests" : {
                        script{
                            tests_started = true
                            sh """#!/usr/bin/env bash
                                source activate $ENV_PATH
                                cd testing
                                ./run_tests.py --threads 12
                            """
                        }
                    },
                    "SonarQube": {
                        withSonarQubeEnv('SonarQube') {
                            script {
                                def sonarScanner = tool name: 'SQscanner', type: 'hudson.plugins.sonar.SonarRunnerInstallation'
                                sh """#!/usr/bin/env bash
                                    source activate $ENV_PATH_J
                                    ${sonarScanner}/bin/sonar-scanner -DskipTests -Dcobertura.skip
                                    source deactivate
                                """
                            }
                        }
                    }
                )
            }
        }
    }
    post {
        cleanup{
            sh """#!/usr/bin/env bash
                rm $ENV_PATH -r
                rm $ENV_PATH_J -r
                rm $JULIA_DIR -r
                rm $MINI_PATH -r
                rm testing/testingdata/reference -r
                rm testing/testingdata/calculated -r
                rm testing/testing_res.txt
                rm testing/cov_*.log
                rm testing/LEFT_RULES
                rm testing/ALL_RULES
                rm testing/USED_RULES
                rm smk-cov.xml
            """
        }
    }
}