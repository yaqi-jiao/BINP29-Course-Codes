# Long read de novo assembly exercise

**Description**: This repository includes the workflow of long read data denovo assembly.

**The workflow includes:**

    1. quality control
    2. de novo assembly
    3. Quast quality assessment
    4. BUSCO assessment
   
## 1. Data structure

```text
.
├── Result
│   ├── 00_fastqc
│   ├── 01_denovo_assembly
│   └── 02_BUSCO
│ 
├── Script
│
└── README.md
```

## 2. Requirements

The following softwares and packages are required:

- FastQc v0.11.9
- Hifiasm-0.25.0-r726
- gfatools 0.5.5
- quast 5.3.0
- BUSCO 6.0.0
  
## 3. workflow

### Step1: Fastqc

```bash
mkdir 00_fastqc
cd 00_fastqc/
fastqc ../raw_data/SRR13577846.fastq.gz -o .
```
Output html file can be found in `Result/00_fastqc`

### Step2: Assembly

```bash
mkdir 01_denovo_assembly
cd 01_denovo_assembly
/usr/bin/time -v hifiasm -o asm -t 10 ../raw_data/SRR13577846.fastq.gz 2> hifiasm.time.txt
gfatools gfa2fa asm.bp.p_ctg.gfa > asm.bp.p_ctg.fa
```
Running time: System time (seconds) 90.85

Output file can be found in `Result/01_denovo_assembly`

### Step3: Summary Statistics using quast

```bash
usr/bin/time -v quast.py asm.bp.p_ctg.fa -o ../Result/01_denovo_assembly/quast_report
```

Running time: System time (seconds) 0.92

Output report file can be found in `Result/01_denovo_assebly/quast_report`


### Step4: BUSCO summary

```bash
/usr/bin/time -v busco -i ../01_denovo_assembly/asm.bp.p_ctg.fa -l saccharomycetes_odb10 -o ../Result/02_BUSCO/busco_asm -m genome -c 10
```

Running time: System time (seconds) 56.33

Output file can be fount in `Result/02_BUSCO/busco_asm`




