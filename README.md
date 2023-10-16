# pub-gcr
This is the repository that we use to host some public docker images for utility use

<details>
  <summary><b>mpileup2matrix</b> &#x1F4D9; </summary>
  
  ### What does it do?
  mpileup2matrix is a docker image that takes a list of input fastq files from Nanopore sequencer and trims and aligns them against a reference sequence. It will then generate an mpileup file (*.mpileup) and two matrices: one is the coverage matrix and the other is the indel matrix, both are table delimited and on a per position basis.
  
  ### How to run it?
  
    1. Install docker (if you haven't done it) [link to installation page](https://docs.docker.com/engine/install/)
    2. Install git (if you haven't done it)
    3. Run `git clone` of this repository:
       ```bash
       git clone 
</details>
