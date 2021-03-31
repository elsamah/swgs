from snakemake.utils import validate
import pandas as pd

# this container defines the underlying OS for each job when using the workflow
# with --use-conda --use-singularity
container: "docker://continuumio/miniconda3"

##### load config and sample sheets #####

configfile: "config/config.yaml"
validate(config, schema="../schemas/config.schema.yaml")

samples = pd.read_csv(config["samples"], sep="\t").set_index("sample", drop=False)
samples.index.names = ["sample_id"]
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_csv(
    config["units"], dtype=str, sep="\t").set_index(["sample", "unit"], drop=False)
units.index.names = ["sample_id", "unit_id"]
units.index = units.index.set_levels(
    [i.astype(str) for i in units.index.levels])  # enforce str in index
validate(units, schema="../schemas/units.schema.yaml")

report: "../report/workflow.rst"

##### wildcard constraints #####

#wildcard_constraints:
#    sample = "|".join(samples.index),
#    unit = "|".join(units["unit"])

####### helpers ###########
def get_fastqs(wildcards):
    """Get raw FASTQ files from unit sheet."""
    u = units.loc[ (wildcards.sample, wildcards.unit), ["fq1", "fq2"] ].dropna()
    return [ f"{u.fq1}", f"{u.fq2}" ]
