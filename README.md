# pub-gcr
This is the repository that we use to host some public docker images for utility use

<details>
  <summary>
    
  ## mpileup2matrix &#x1F4D9; 
  
  </summary>
  
  ### What does it do?
  mpileup2matrix is a docker image that takes a list of input fastq files from Nanopore sequencer and trims and aligns them against a reference sequence. It will then generate an mpileup file (*.mpileup) and two matrices: one is the coverage matrix and the other is the indel matrix, both are table delimited and on a per position basis.
  
  ### How to install it?
  
  <b>Step 1</b>
  Install docker (if you haven't done it) [link to installation page](https://docs.docker.com/engine/install/)
    
  <b>Step 2</b>
  Install git (if you haven't done it) [link to installation page](https://docs.github.com/en/desktop/installing-and-authenticating-to-github-desktop/installing-github-desktop)

  <b>Step 3</b>
  Run `git clone` of this repository:
       
  ```bash
    gh repo clone quantumsky-lab/pub-gcr
  ```

  or

  ```bash
    git clone https://github.com/quantumsky-lab/pub-gcr.git
  ```

  <b>Step 4</b>
  Use `cd` to nagivate to `pub-gcr/mpileup2matrix` and run:

  ```bash
    docker build -t mpileup2matrix .
  ```

  If you are using an Apple Silicon device (such as M1/2 chips), then you should run:

  ```bash
    docker buildx build --platform linux/amd64 -t mpileup2matrix .
  ```

  ### How to run it?

  You can get the helper information by running:

  ```bash
    docker run --rm mpileup2matrix -h
  ```

  You will get a print message that looks like this:

    usage: mpileup2matrix.py [-h] --infile-list INFILE_LIST --infile-vol INFILE_VOL --reference REFERENCE [--temp-dir TEMP_DIR] [--keep-temp] --prefix PREFIX [--blastn BLASTN]
                         [--makeblastdb MAKEBLASTDB] [--trimmomatic TRIMMOMATIC] [--homopolymer HOMOPOLYMER] [--min-map MIN_MAP]

    Run reads mapping with Jorna default settings for genome editing data.
    
    optional arguments:
      -h, --help            show this help message and exit
      --infile-list INFILE_LIST, -i INFILE_LIST
                            A list of input files in fastq format in text file; if paired-end, they should be in the same line, separated by comma. NOTE: no directory should be supplied.
      --infile-vol INFILE_VOL, -d INFILE_VOL
                            directory where the infiles are stored
      --reference REFERENCE, -r REFERENCE
                            Reference sequence in fasta format
      --temp-dir TEMP_DIR, -t TEMP_DIR
                            Where the intermediate files should live.
      --keep-temp, -k       Turns on temp dir keeping when specified.
      --prefix PREFIX, -o PREFIX
                            Output file prefix
      --blastn BLASTN       Path to blastn
      --makeblastdb MAKEBLASTDB
                            Path to makeblastdb
      --trimmomatic TRIMMOMATIC
                            Path to trimmomatic jar
      --homopolymer HOMOPOLYMER
                            Homopolymer length threshold
      --min-map MIN_MAP     Mininum match length threshold

  Please follow the example shown below to learn how to run the image.
  
  ### Example

  In the repository,  you will find a folder named `test`, which contains necessary files that you will use to do a test run. The files include:
    
    
      test
      ├── INIP.fa
      ├── INIP_samples.txt
      └── data
          ├── B1_01.fastq
          ├── B2_02.fastq
          ├── B3_03.fastq
          └── B4_04.fastq
      
      1 directory, 7 files
      
  `INIP.fa` is the reference sequence in fasta format. `INIP_samples.txt` is a plain text file that contains the `*.fastq` files, one file per line. 

  The `data` folder, is where all the raw `*.fastq` files are stored. 

  Note that these files are on your local drive. To run the docker image, we will use the `-v` option in docker to mount the local directory to become a virtual directory on the container. The way to do it is:

  ```bash
    docker run --rm -v $PWD:/root/data mpileup2matrix [options]
  ```
 
</details>
