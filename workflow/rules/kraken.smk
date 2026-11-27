# rules for kraken-based analysis

rule downloadtaxonomy_kraken:
    """
    download taxonomy for kraken to create a custom database
    """
    output:
        flag = touch("CustomFishDB/taxonomy_step.done")
    params:
        db = "CustomFishDB"
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2-build --download-taxonomy --db {params.db} --use-ftp"

rule add_sequence_to_db:
    """
    add fasta file of BOLD to library for after compilation of the database
    """
    input:
        flag = "CustomFishDB/taxonomy_step.done",
        bold = config["BOLD"]
    output:
        flag = touch("CustomFishDB/library_step.done")
    params:
        db = "CustomFishDB"
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2-build --add-to-library {input.bold} --db {params.db}"

rule build_db:
    """
    compilation of database for kraken
    """
    input:
        flag = "CustomFishDB/library_step.done"
    output:
        flag = touch("CustomFishDB/build_step.done")
    params:
        db = "CustomFishDB"
    threads: 8
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2-build --build --threads {threads} --db {params.db}"



# if you have already created the database, the pipeline will skip the steps above and go directly to the kraken analysis. If you want to re-create the database, 
# delete the CustomFishDB folder and re-run the pipeline.

rule run_kraken:
    """
    run Kraken on custom database with reads from PoolSeq data
    """
    input:
        unpack(get_fastqs),
        db_flag = "CustomFishDB/build_step.done",
        checksums = INITIAL_CHECKSUMS
    output:
        report = OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_kraken.report",
        out = temp(OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_kraken.out")
    params:
        db = "CustomFishDB",
        conf = config["confidence"]
    threads: 10
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --confidence {params.conf} --gzip-compressed --paired {input.r1} {input.r2} --report {output.report} --output {output.out}"