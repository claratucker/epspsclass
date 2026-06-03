# workflow/rules/epspsclass.smk
# ==============================
# Snakemake rule for EPSPS classification with EPSPSClass.
#
# Usage in your Snakefile:
#   include: "path/to/epspsclass.smk"
#
# Required config keys (set in config.yaml or --config):
#   input_fasta:  path to protein FASTA (can be wildcarded)
#   results_dir:  directory for output TSV files
#
# Optional:
#   epspsclass_threshold: identity threshold (default 40.0)
#
# Example config.yaml:
#   input_fasta: "data/epsps_sequences.fasta"
#   results_dir: "results/classification"
#   epspsclass_threshold: 40.0
#
# Example wildcard usage:
#   expand("{results_dir}/{sample}_epspsclass.tsv",
#          results_dir=config["results_dir"],
#          sample=SAMPLES)

rule epspsclass_classify:
    """Classify EPSPS sequences by glyphosate sensitivity class."""
    input:
        fasta = "{prefix}.fasta",
    output:
        tsv   = "{prefix}_epspsclass.tsv",
    log:
        "{prefix}_epspsclass.log",
    params:
        threshold = config.get("epspsclass_threshold", 40.0),
    conda:
        "../envs/epspsclass.yaml"
    shell:
        """
        epspsclass classify \
            --input  {input.fasta} \
            --output {output.tsv} \
            --threshold {params.threshold} \
            --summary \
            2> {log}
        """


rule epspsclass_classify_by_exposure_tier:
    """
    Classify sequences split by exposure tier (Tier 1/2/3).
    Expects input FASTAs named: {tier}_sequences.fasta
    Produces per-tier classification TSVs for downstream dN/dS stratification.
    """
    input:
        fasta = "data/{tier}_sequences.fasta",
    output:
        tsv   = "results/classification/{tier}_epspsclass.tsv",
    log:
        "logs/epspsclass_{tier}.log",
    params:
        threshold = config.get("epspsclass_threshold", 40.0),
    conda:
        "../envs/epspsclass.yaml"
    shell:
        """
        epspsclass classify \
            --input  {input.fasta} \
            --output {output.tsv} \
            --threshold {params.threshold} \
            --summary \
            2> {log}
        """


rule epspsclass_validate:
    """Print marker table to log — run once to confirm installation."""
    output:
        touch("results/.epspsclass_validated"),
    log:
        "logs/epspsclass_validate.log",
    conda:
        "../envs/epspsclass.yaml"
    shell:
        "epspsclass validate-markers > {log} 2>&1"
