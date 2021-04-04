from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####
# Config file
configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

# Samples: List of samples and conditions
samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

# List of sample+unit information (e.g. paths, builds, etc.)
units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### test space #####
#u = units.loc[ ('net-037', '1'), ["fq1", "fq2"] ].dropna()
#print(u)
#print([ f"{u.fq1}", f"{u.fq2}" ])
#print("|".join(samples.index))
#print("|".join(units["unit"]))

##### wildcard constraints #####
wildcard_constraints:
    sample = "|".join(samples.index),
    unit = "|".join(units["unit"])

####### helpers ###########
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    #u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
    u = units.loc[ (wildcards.sample, '1'), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq1}", f"{u.fq2}" ]

def get_snp_paths(wildcards):
    ''' Assembles the paths for snp/indels for indel realignment'''
    var_build = "snp_" + config['common']['build']
    print(var_build)
    known = list(config['params']['gatk'][var_build].values())
    return known

def get_indel_paths(wildcards):
    ''' Assembles the paths for snp/indels for indel realignment'''
    var_build = "indel_" + config['common']['build']
    print(var_build)
    known = list(config['params']['gatk'][var_build].values())
    return known

def get_samples():
    print(samples.index.tolist())
    return samples.index.tolist()

def combine_args(input_args):
    format_args = " ".join(input_args)
    return format_args

def get_rgid(wildcards):
    """ Files in a raw @RG header for bwa mem alignment """
    dat = units.loc[ (wildcards.sample, '1'), ['platform', 'library'] ].dropna()
    rg=("@RG" +
        "\\tID:" + wildcards.sample +
        "\\tSM:" + wildcards.sample +
        "\\tPL:" + f"{dat.platform}" +
        "\\tPU:L001" +
        "\\tLB:" + f"{dat.library}")
    rg = "-R '" + rg + "'"
    return f"{rg}"
