import hashlib
from pathlib import Path


def calculate_sha256(filepath):
    sha256 = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


# Initial checksum calculation
rule initial_checksums:
    input:
        unpack(get_fastqs)
    output:
        INITIAL_CHECKSUMS
    run:
        Path(output[0]).parent.mkdir(parents=True, exist_ok=True)
        with open(output[0], "w") as out:
            for filepath in [input.r1, input.r2]:
                checksum = calculate_sha256(filepath)
                out.write(f"{checksum}  {filepath}\n")


# Final integrity verification
rule verify_input_integrity:
    input:
        unpack(get_fastqs),
        initial = INITIAL_CHECKSUMS,
        results = TARGETS
    output:
        VERIFIED_CHECKSUMS
    run:
        initial_checksums = {}
        with open(input.initial) as f:
            for line in f:
                checksum, filepath = line.rstrip("\n").split("  ", 1)
                initial_checksums[filepath] = checksum

        modified_files = []
        for filepath in [input.r1, input.r2]:
            current_checksum = calculate_sha256(filepath)
            if filepath not in initial_checksums or current_checksum != initial_checksums[filepath]:
                modified_files.append(filepath)

        if modified_files:
            raise ValueError(
                "Input FASTQ files have been modified:\n"
                + "\n".join(modified_files)
            )

        with open(output[0], "w") as f:
            f.write("All input FASTQ files are unchanged.\n")

        print(f"All input FASTQ files passed SHA-256 verification for {wildcards.species}_{wildcards.sex}.")