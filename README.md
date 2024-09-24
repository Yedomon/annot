# Step 0: Repeats


```



#!/usr/bin/env python3

import os
import argparse
import subprocess

def build_database(db, input_file):
    print("Building database ...")
    subprocess.run(["BuildDatabase", "-name", db, input_file], check=True)

def run_repeat_modeler(db, threads):
    print("Running RepeatModeler ...")
    subprocess.run(["RepeatModeler", "-database", db, "-LTRStruct", "-threads", str(threads)], check=True)

def run_repeat_masker(input_file, pa):
    print("Running RepeatMasker ...")
    subprocess.run(["RepeatMasker", "-pa", str(pa), "-gff", "-lib", "consensi.fa", input_file], check=True)

def main():
    parser = argparse.ArgumentParser(description="Script to run RepeatModeler and RepeatMasker on a genome.")
    
    parser.add_argument("-i", "--input", required=True, help="Input genome FASTA file")
    parser.add_argument("-t", "--threads", type=int, default=16, help="Number of threads to use for RepeatModeler (default: 16)")
    parser.add_argument("-p", "--pa", type=int, default=16, help="Number of processors for RepeatMasker (default: 16)")
    parser.add_argument("-d", "--database", help="Name of the database (default: basename of input file without extension)")

    args = parser.parse_args()

    input_file = args.input
    db = args.database if args.database else os.path.splitext(os.path.basename(input_file))[0]
    threads = args.threads
    pa = args.pa

    build_database(db, input_file)
    run_repeat_modeler(db, threads)

    # Link the consensi.fa.classified file
    consensi_file = "consensi.fa"
    if not os.path.exists(consensi_file):
        consensi_dir = next((d for d in os.listdir(".") if d.startswith("RM_")), None)
        if consensi_dir:
            os.symlink(os.path.join(consensi_dir, "consensi.fa"), consensi_file)
        else:
            print("Error: Could not find consensi.fa file.")
            return

    run_repeat_masker(input_file, pa)

if __name__ == "__main__":
    main()



```


# Step 1: Data trimming

trim.sh

```
mkdir trimmed

for i in *_1.fastq.gz

do

base=$(basename $i _1.fastq.gz)

fastp --detect_adapter_for_pe \
        --overrepresentation_analysis \
        --correction --cut_right --thread 2 \
        --html trimmed/${base}.fastp.html --json trimmed/${base}.fastp.json \
        -i ${base}_1.fastq.gz -I ${base}_2.fastq.gz \
        -o trimmed/Trimmed_${base}_1.fq.gz -O trimmed/Trimmed_${base}_2.fq.gz

done

```

# Step 2: Mapping





Simple


```
#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
from pathlib import Path

def parse_arguments():
    """
    Parse the command line arguments passed to the script.
    """
    parser = argparse.ArgumentParser(
        description="A script to perform genome alignment using STAR."
    )

    parser.add_argument(
        "--genomeDir", 
        type=str, 
        required=True, 
        help="Path to the STAR genome directory."
    )
    
    parser.add_argument(
        "--fastaFile", 
        type=str, 
        required=True, 
        help="Path to the masked genome FASTA file."
    )
    
    parser.add_argument(
        "--workDir", 
        type=str, 
        required=True, 
        help="Path to the directory containing the paired-end read files."
    )

    parser.add_argument(
        "--threads", 
        type=int, 
        default=10, 
        help="Number of threads to use for STAR alignment (default: 10)."
    )

    return parser.parse_args()

def generate_star_index(genomeDir, fastaFile, threads):
    """
    Generate STAR genome index.
    """
    print(f"Generating STAR index at {genomeDir} using {fastaFile}")
    subprocess.run([
        "STAR", "--runMode", "genomeGenerate",
        "--genomeDir", genomeDir,
        "--genomeFastaFiles", fastaFile,
        "--runThreadN", str(threads),
        "--genomeSAindexNbases", "12"
    ], check=True)
    print("Genome indexing completed.")

def perform_star_alignment(read1, read2, genomeDir, workDir, sampleName, threads):
    """
    Perform the alignment of paired-end reads using STAR.
    """
    print(f"Aligning sample: {sampleName}")
    
    outPrefix = os.path.join(workDir, f"{sampleName}_")
    
    subprocess.run([
        "STAR", "--genomeDir", genomeDir,
        "--readFilesIn", read1, read2,
        "--readFilesCommand", "zcat",
        "--outSAMstrandField", "intronMotif",
        "--runThreadN", str(threads),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outReadsUnmapped", "Fastx",
        "--twopassMode", "Basic",
        "--limitBAMsortRAM", "150000000000",
        "--outFileNamePrefix", outPrefix,
        "--outFilterScoreMinOverLread", "0",
        "--outFilterMatchNmin", "0",
        "--outFilterMatchNminOverLread", "0"
    ], check=True)
    
    print(f"Alignment completed for {sampleName}")

def main():
    # Parse the command line arguments
    args = parse_arguments()

    # Ensure the working directory exists
    if not os.path.exists(args.workDir):
        raise FileNotFoundError(f"Working directory {args.workDir} does not exist.")

    # Step 1: Generate the STAR genome index
    generate_star_index(args.genomeDir, args.fastaFile, args.threads)

    # Step 2: Perform STAR alignment for each pair of FASTQ files
    input_pattern = os.path.join(args.workDir, '*_1.fq.gz')
    input_files_r1 = glob.glob(input_pattern)

    if len(input_files_r1) == 0:
        raise FileNotFoundError(f"No paired-end files matching pattern {input_pattern}")

    for read1 in input_files_r1:
        read2 = read1.replace('_1.fq.gz', '_2.fq.gz')

        if not os.path.exists(read2):
            print(f"Warning: Paired file for {read1} not found, skipping...")
            continue

        sample_name = Path(read1).stem.replace('_1', '')  # Extract the sample name
        perform_star_alignment(read1, read2, args.genomeDir, args.workDir, sample_name, args.threads)

if __name__ == "__main__":
    main()


```



Run



```
python mapping.py --genomeDir /path/to/genomeDir \
                  --fastaFile /path/to/masked.fasta \
                  --workDir /path/to/readsDir \
                  --threads 20


```


My Run



```

source /home/africabp/00.apps/miniconda3/bin/activate

conda activate star_env

cd /home/africabp/02.analyses/01.mapping

python3 mapping.py --genomeDir /home/africabp/02.analyses/01.mapping --fastaFile chr4.fasta --workDir /home/africabp/02.analyses/01.mapping/reads_trimmed/ --threads 10 &> log.txt &


```




Note: You can put just put the trimmed data in the work dir simply then

Singularity

```
#!/usr/bin/env python3

import os
import glob
import argparse
import subprocess
from pathlib import Path

def parse_arguments():
    """
    Parse the command line arguments passed to the script.
    """
    parser = argparse.ArgumentParser(
        description="A script to perform genome alignment using STAR through Singularity."
    )

    parser.add_argument(
        "--singularityImage", 
        type=str, 
        required=True, 
        help="Path to the STAR Singularity image (e.g., /path/to/star.sif)."
    )

    parser.add_argument(
        "--genomeDir", 
        type=str, 
        required=True, 
        help="Path to the STAR genome directory."
    )
    
    parser.add_argument(
        "--fastaFile", 
        type=str, 
        required=True, 
        help="Path to the masked genome FASTA file."
    )
    
    parser.add_argument(
        "--workDir", 
        type=str, 
        required=True, 
        help="Path to the directory containing the paired-end read files."
    )

    parser.add_argument(
        "--threads", 
        type=int, 
        default=10, 
        help="Number of threads to use for STAR alignment (default: 10)."
    )

    return parser.parse_args()

def run_singularity_star_command(singularity_image, star_args):
    """
    Run STAR command inside the Singularity container.
    """
    singularity_command = ["singularity", "exec", singularity_image, "STAR"] + star_args
    subprocess.run(singularity_command, check=True)

def generate_star_index(singularity_image, genomeDir, fastaFile, threads):
    """
    Generate STAR genome index using Singularity.
    """
    print(f"Generating STAR index at {genomeDir} using {fastaFile}")
    
    star_args = [
        "--runMode", "genomeGenerate",
        "--genomeDir", genomeDir,
        "--genomeFastaFiles", fastaFile,
        "--runThreadN", str(threads),
        "--genomeSAindexNbases", "12"
    ]
    
    run_singularity_star_command(singularity_image, star_args)
    
    print("Genome indexing completed.")

def perform_star_alignment(singularity_image, read1, read2, genomeDir, workDir, sampleName, threads):
    """
    Perform the alignment of paired-end reads using STAR in Singularity.
    """
    print(f"Aligning sample: {sampleName}")
    
    outPrefix = os.path.join(workDir, f"{sampleName}_")
    
    star_args = [
        "--genomeDir", genomeDir,
        "--readFilesIn", read1, read2,
        "--readFilesCommand", "zcat",
        "--outSAMstrandField", "intronMotif",
        "--runThreadN", str(threads),
        "--outSAMtype", "BAM", "SortedByCoordinate",
        "--outReadsUnmapped", "Fastx",
        "--twopassMode", "Basic",
        "--limitBAMsortRAM", "150000000000",
        "--outFileNamePrefix", outPrefix,
        "--outFilterScoreMinOverLread", "0",
        "--outFilterMatchNmin", "0",
        "--outFilterMatchNminOverLread", "0"
    ]
    
    run_singularity_star_command(singularity_image, star_args)
    
    print(f"Alignment completed for {sampleName}")

def main():
    # Parse the command line arguments
    args = parse_arguments()

    # Ensure the working directory exists
    if not os.path.exists(args.workDir):
        raise FileNotFoundError(f"Working directory {args.workDir} does not exist.")

    # Step 1: Generate the STAR genome index
    generate_star_index(args.singularityImage, args.genomeDir, args.fastaFile, args.threads)

    # Step 2: Perform STAR alignment for each pair of FASTQ files
    input_pattern = os.path.join(args.workDir, '*_1.fq.gz')
    input_files_r1 = glob.glob(input_pattern)

    if len(input_files_r1) == 0:
        raise FileNotFoundError(f"No paired-end files matching pattern {input_pattern}")

    for read1 in input_files_r1:
        read2 = read1.replace('_1.fq.gz', '_2.fq.gz')

        if not os.path.exists(read2):
            print(f"Warning: Paired file for {read1} not found, skipping...")
            continue

        sample_name = Path(read1).stem.replace('_1', '')  # Extract the sample name
        perform_star_alignment(args.singularityImage, read1, read2, args.genomeDir, args.workDir, sample_name, args.threads)

if __name__ == "__main__":
    main()
```




Run

```

python mapping.py --singularityImage /path/to/star.sif \
                  --genomeDir /path/to/genomeDir \
                  --fastaFile /path/to/masked.fasta \
                  --workDir /path/to/readsDir \
                  --threads 20



```

#Step 3: Braker annotation
