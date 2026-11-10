configfile: "config/config.yaml"

SAMPLES = config["samples"]
DB_DIR = config["db_dir"]

rule all:
    input:
        expand("results/kraken2/{sample}_report.txt", sample=SAMPLES)
        expand("results/summary/{sample}_contamination.txt", sample=SAMPLES)

rule download_taxonomy:
    output:
        directory(f"{DB_DIR}/taxonomy")
    log:
        "logs/download_taxonomy.log"
    shell:
        """
        kraken2-build --download-taxonomy --db {DB_DIR} > {log} 2>&1
        """
    
rule add_bold_to_library:
    input:
        fasta= config["bold_fasta"]
        taxonomy= f"{DB_DIR}/taxonomy"
    output:
        flag = touch(f"{DB_DIR}/library_added.flag")
    log:
        "logs/add_library.log"
    shell:
        """
        kraken2-build --add-to-library {input.fasta} --db {DB_DIR} > {log} 2>&1
        """

rule build_kraken_db:
    input:
        flag = f"{DB_DIR}/library_added.flag"
    output:
        db_file = f"{DB_DIR}/hash.k2d"
    threads: 8
    log:
        "logs/build_db.log"
    shell:
        """
        kraken2-build --build --db {DB_DIR} --threads {threads} > {log} 2>&1
        """

rule kraken2_classify:
    input:
        read = "reads/{sample}.fastq",
        db = f"{DB_DIR}/hash.k2d"
    output:
        report = "results/kraken2/{sample}_report.txt",
        output = "results/kraken2/{sample}_kraken.txt"
    threads: 8
    params:
        db = DB_DIR
    log:
        "logs/kraken2/{sample}.log"
    shell:
        """
        kraken2 --db {params.db} \
                --threads {threads} \
                {input.read} \
                --report {output.report} \
                --output {output.output} \
                {input.read} > {log} 2>&1
        """


rule parse_contamination_report:
    input:
        report = "results/kraken2/{sample}_report.txt"
    output:
        summary = "results/summary/{sample}_contamination.txt"
    shell:
        """
        echo "=== Bilan de classification : {wildcards.sample} ===" > {output.summary}
        echo -e "% Reads\tNbre Reads\tNiveau\tTaxon" >> {output.summary}
        echo "----------------------------------------------------" >> {output.summary}
        awk -F'\\t' '$4 ~ /U|S/ {{printf "%.2f%%\t%s\t%s\t%s\\n", $1, $2, $4, $6}}' {input.report} >> {output.summary}
        """

rule format_bold_fasta:
    """
    Formatte le fichier FASTA BOLD d'origine en ajoutant les TaxIDs NCBI aux entêtes
    """
    input:
        raw_fasta = "resources/bold_raw.fasta"
    output:
        formatted_fasta = "resources/bold_references_formatted.fasta"
    log:
        "logs/format_bold.log"
    shell:
        """
        python3 format_bold_kraken.py {input.raw_fasta} {output.formatted_fasta} > {log} 2>&1
        """