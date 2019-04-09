"""
Snakefile for  genotyping

Author: Marc Jan Bonder
Affiliation: EMBL-EBI
Date: Tuesday 10 February 2019
#Run: snakemake --snakefile ./snakemake --jobs 50 --latency-wait 30 --cluster-config ../cluster.json --cluster 'bsub -q {cluster.queue} -n {cluster.n} -R "rusage[mem={cluster.memory}]" -M {cluster.memory} -o ./cis_eqlts.20181120.o -e ./cis_eqlts.20181120.e'
"""

import glob									# Unix style pathname pattern expansion. 
import os									# Allows interaction with UNIX filesystem
from subprocess import run					# Allows you to run subprocesses
import pandas as pd							# data manipulation and analysis
import re									# regex
from os.path import join					# Join file paths intelligently

shell.prefix("set -euo pipefail;") 

def _multi_arg_start(flag, files):
    flag += " "
    return " ".join(flag + f for f in files)

def _multi_arg_end(flag, files):
    flag = " "+flag
    return " ".join(f + flag for f in files)

def _multi_arg_both_ends(flag1, flag2, files):
    flag1 += " "
    flag2 = " "+flag2
    return " ".join(flag1 + f + flag2 for f in files)

def flatenChunk(chunk):
    return chunk.replace(":", "_").replace("-", "_")

def extendChunk(chunk):
    relChunk = chunk.pop()
    chunkSplitted = relChunk.split("_")
    return chunkSplitted[0]+":"+chunkSplitted[1]+"-"+chunkSplitted[2]

#Variables
chunkFile = './chr.txt'
genotypeFile = './hipsci_pipeline/geuvadis_CEU_test_data/Genotypes/Geuvadis' #Can be .bed or .bgen, no extenion needed 
annotationFile = './hipsci_pipeline/geuvadis_CEU_test_data/Expression/Geuvadis_CEU_Annot.txt'
phenotypeFile = './hipsci_pipeline/geuvadis_CEU_test_data/Expression/phenos.txt'
covariateFile = './hipsci_pipeline/geuvadis_CEU_test_data/Expression/cov.txt'
kinshipFile = './kinship.rel'
sampleMappingFile = './mappingFile.txt' #Not needed
numberOfPermutations = '10000'
minorAlleleFrequency = '0.000000001'
hwe = '0.00001'
callRate = '0.95'
windowSize = '250000'
blockSize = '1500'
outputFolder = './Rel/'

finalQtlRun = './Rel/qtl_results_all.txt'

with open(chunkFile,'r') as f:
    chunks = [x.strip() for x in f.readlines()]

qtlOutput = []
for chunk in chunks:
    #print(chunk)
    processedChunk = flatenChunk(chunk)
    #print(processedChunk)
    processedChunk=expand(outputFolder+'{chunk}.finished',chunk=processedChunk )
    qtlOutput.append(processedChunk)


## flatten these lists
qtlOutput = [filename for elem in qtlOutput for filename in elem]
#finalQtlRun = [filename for elem in finalQtlRun for filename in elem]

rule all:
    input:
        qtlOutput, finalQtlRun

rule run_qtl_mapping:
    input:
        af = annotationFile,
        pf = phenotypeFile,
        cf = covariateFile,
        kf = kinshipFile,
        smf = sampleMappingFile
    output:
        './Rel/{chunk}.finished'
    params:
        gen=genotypeFile,
        od = outputFolder,
        np = numberOfPermutations,
        maf = minorAlleleFrequency,
        hwe = hwe,
        cr = callRate,
        w = windowSize,
        bs = blockSize
    run:
        chunkFull = extendChunk({wildcards.chunk})
        shell(
            "python ./hipsci_pipeline/limix_QTL_pipeline/run_QTL_analysis.py "
            "--plink {params.gen} "
            " -af {input.af} "
            " -pf {input.pf} "
            " -cf {input.cf} "
            " -od {params.od} "
            " --kinship_file {input.kf} "
            " --sample_mapping_file {input.smf} "
            " -gr {chunkFull} "
            " -np {params.np} "
            " -maf {params.maf} "
            " -hwe {params.hwe} "
            " -cr {params.cr} "
            " -c -gm standardize "
            " -w {params.w} "
            " --block_size {params.bs} ")
        shell("touch {output}")


rule aggregate_qtl_results:
    input:
        IF = outputFolder,
        OF = outputFolder,
        finalFiles = qtlOutput
    output:
        './Rel/qtl_results_all.txt'
    run:
        shell(
            "python ./hipsci_pipeline/post-processing_QTL/minimal_postprocess.py "
            "-id {input.IF} "
            " -od {input.OF} "
            " -sfo ")
