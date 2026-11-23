KMER_DIR = WDIR + "/{species}/Poolseq/MGX/check_contamination/kmer/{sex}"

rule make_blast_db:
    "creation of blast database"
    input:
        bold = config["BOLD"]
    output:
        flag = touch(CustomBlastBoldDB/build_db.done)
    params:
        db = "CustomBlastBoldDB"
    conda:
        "../envs/kmer.yaml"
    shell:
        "Makeblastdb -in {input.bold} -dbtype nucl -out {params.db}"

rule kmtricks_pipeline:
    input:
        unpack(get_fastqs)
        fof = "fof.txt"
    output:
        matrix = KMER_DIR + "/resultats_kmtricks/matrices/matrix_0.count"
    params:
        run_dir = KMER_DIR + "/resultats_kmtricks"
        size_kmer = config["kmer_size"]
    threads: 8
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        rm -rf {params.run_dir}
        kmtricks pipeline --file {input.fof} \
            --run-dir {params.run_dir} \
            --kmer-size {params.size_kmer}
            --hard-min 1 \
            --nb_partitions 1 \
            --threads {threads} \
            --mode kmer:count:bin
        """

rule kmtricks_dumb:
    "extract kmer to fasta file"
    input:
        matrix = KMER_DIR + "/resultats_kmtricks/matrices/matrix_0.count"
    output:
        fasta = KMER_DIR + "/kmers_clean.fasta"
    params:
        run_dir = KMER_DIR + "/resultats_kmtricks"
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        kmtricks dump --run-dir {params.run_dir} --input {input.matrix} \
            | awk '{{print ">"NR"\\n"$1}}' > {output.fasta}
        """

def write_read_lists(r1_in, r2_in, in_txt, out_txt, r1_out, r2_out):
    with open(in_txt, "w") as f_in:
        f_in.write(f"{r1_in}\n{r2_in}\n")
    with open(out_txt, "w") as f_out:
        f_out.write(f"{r1_out}\n{r2_out}\n")

rule prepare_back2seq_lists:
    input:
        unpack(get_fastqs)
    output:
        in_list = temp(KMER_DIR + "/in_reads.txt"),
        out_list = temp(KMER_DIR + "/out_reads.txt")
    params:
        r1_out = KMER_DIR + "/reads_extract_1.fastq",
        r2_out = KMER_DIR + "/reads_extract_2.fastq"
    run:
        write_read_lists(
            r1_in=input.r1,
            r2_in=input.r2,
            in_txt=output.in_list,
            out_txt=output.out_list,
            r1_out=params.r1_out,
            r2_out=params.r2_out
        )

rule back_to_sequences:
    """
    Récupère les reads contenant les k-mers extraits
    """
    input:
        unpack(get_fastqs),
        kmers = KMER_DIR + "/kmers_propres.fasta"
        in_list = KMER_DIR + "/in_reads.txt",
        out_list = KMER_DIR + "/out_reads.txt"
    output:
        r1 = KMER_DIR + "/reads_extract_1.fastq",
        r2 = KMER_DIR + "/reads_extract_2.fastq"
    params:
        kmer_size = config["kmer_size"]
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        back_to_sequences --in-kmers {input.kmers} \
            --in-filelist {output.in_list} \
            --out-filelist {output.out_list} \
            -k {params.kmer_size}
        """


rule spades_assembly:
    """
    Assemble les reads extraits avec SPAdes
    """
    input:
        r1 = KMER_DIR + "/reads_extract_1.fastq",
        r2 = KMER_DIR + "/reads_extract_2.fastq"
    output:
        contigs = KMER_DIR + "/assemblage_final/contigs.fasta"
    params:
        outdir = KMER_DIR + "/assemblage_final"
    threads: 10
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        rm -rf {params.outdir}
        spades.py --only-assembler --careful \
            --threads {threads} \
            --s1 {input.r1} --s1 {input.r2} \
            -o {params.outdir}
        """


rule blast_contigs:
    """
    Aligne les contigs assemblés contre la base BOLD
    """
    input:
        contigs = KMER_DIR + "/assemblage_final/contigs.fasta",
        db = BLAST_DB_FILES
    output:
        blast = WDIR + "/{species}/Poolseq/MGX/results/kmers/{species}_{sex}_blast_output.tsv"
    params:
        db_prefix = "db/CustomBold1DB"
    threads: 10
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        blastn -query {input.contigs} \
            -db {params.db_prefix} \
            -num_threads {threads} \
            -outfmt 6 > {output.blast}