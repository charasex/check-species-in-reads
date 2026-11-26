# rules for kmer-based analysis

# directory for kmers results
KMER_DIR = OUTDIR + "/{species}/PoolSeq_check/kmers/{sex}"

rule make_blast_db:
    "creation of custom blast database"
    input:
        bold = config["BOLD_blast"]
    output:
        flag = touch("dbBlast/build_db.done")
    params:
        db = "dbBlast/CustomBlastBoldDB"
    conda:
        "../envs/kmer.yaml"
    shell:
        "makeblastdb -in {input.bold} -dbtype nucl -out {params.db}"

rule kmtricks_pipeline:
    "run kmtricks pipeline to extract kmers from BOLD database"
    input:
        unpack(get_fastqs),
        fof = "fof.txt"
    output:
        matrix = temp(KMER_DIR + "/resultats_kmtricks/matrices/matrix_0.count")
    params:
        run_dir = KMER_DIR + "/resultats_kmtricks",
        size_kmer = config["kmer_size"]
    threads: 8
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        rm -rf {params.run_dir}
        kmtricks pipeline --file {input.fof} \
            --run-dir {params.run_dir} \
            --kmer-size {params.size_kmer} \
            --hard-min 1 \
            --nb-partitions 1 \
            --threads {threads} \
            --mode kmer:count:bin
        """

rule kmtricks_dumb:
    "extract kmer to fasta file"
    input:
        matrix = KMER_DIR + "/resultats_kmtricks/matrices/matrix_0.count"
    output:
        fasta = temp(KMER_DIR + "/kmers_clean.fasta")
    params:
        run_dir = KMER_DIR + "/resultats_kmtricks"
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        kmtricks dump --run-dir {params.run_dir} --input {input.matrix} \
            | awk '{{print ">"NR"\\n"$1}}' > {output.fasta}
        """

# write two double file lists for back_to_sequences, one for input reads and one for output reads

def write_read_lists(r1_in, r2_in, in_txt, out_txt, r1_out, r2_out):
    with open(in_txt, "w") as f_in:
        f_in.write(f"{r1_in}\n{r2_in}\n")
    with open(out_txt, "w") as f_out:
        f_out.write(f"{r1_out}\n{r2_out}\n")

rule prepare_back2seq_lists:
    """
    prepare the input and output file lists for back_to_sequences
    """
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
    find the target reads containing into the kmers extracted from BOLD and write them to new fastq files
    """
    input:
        unpack(get_fastqs),
        kmers = KMER_DIR + "/kmers_clean.fasta",
        in_list = KMER_DIR + "/in_reads.txt",
        out_list = KMER_DIR + "/out_reads.txt"
    output:
        r1 = temp(KMER_DIR + "/reads_extract_1.fastq"),
        r2 = temp(KMER_DIR + "/reads_extract_2.fastq")
    params:
        kmer_size = config["kmer_size"]
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        back_to_sequences --in-kmers {input.kmers} \
            --in-filelist {input.in_list} \
            --out-filelist {input.out_list} \
            -k {params.kmer_size}
        """


rule spades_assembly:
    """
    assemble the reads extracted with back_to_sequences into contigs using SPAdes
    """
    input:
        r1 = KMER_DIR + "/reads_extract_1.fastq",
        r2 = KMER_DIR + "/reads_extract_2.fastq"
    output:
        contigs = temp(KMER_DIR + "/assemblage_final/contigs.fasta")
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
    align the assembled contigs against the BOLD database
    """
    input:
        contigs = KMER_DIR + "/assemblage_final/contigs.fasta",
        db_flag = "dbBlast/build_db.done"
    output:
        blast = OUTDIR + "/{species}/PoolSeq_check/{species}_{sex}_blast_output.tsv"
    params:
        db = "dbBlast/CustomBlastBoldDB"
    threads: 10
    conda:
        "../envs/kmer.yaml"
    shell:
        """
        blastn -query {input.contigs} \
            -db {params.db} \
            -num_threads {threads} \
            -outfmt 6 > {output.blast}
        
        # add header to the blast output file
        sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore' {output.blast}
        """