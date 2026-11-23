rule downloadtaxonomy_kraken:
    """
    download taxonomy for kraken
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
    add fasta file of BOLD to library
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
    compilation of DB
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

rule run_kraken:
    """
    run Kraken
    """
    input:
        unpack(get_fastqs),
        db_flag = "CustomFishDB/build_step.done"
    output:
        report = WDIR + "/{species}/PoolSeq/MGX/check_contamination/{species}_{sex}_kraken.report",
        out = temp(WDIR + "/{species}/PoolSeq/MGX/check_contamination/{species}_{sex}_kraken.out")
    params:
        db = "CustomFishDB",
        conf = config["confidence"]
    threads: 10
    conda:
        "../envs/kraken.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --confidence {params.conf} --gzip-compressed --paired {input.r1} {input.r2} --report {output.report} --output {output.out}"