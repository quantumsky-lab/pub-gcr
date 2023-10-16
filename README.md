# pub-gcr
This is the repository that we use to host some public docker images for utility use

<details>
  <summary>
    
  ## mpileup2matrix &#x1F4D9; 
  
  </summary>
  
  ### What does it do?
  mpileup2matrix is a docker image that takes a list of input fastq files from Nanopore sequencer and trims and aligns them against a reference sequence. It will then generate an mpileup file (*.mpileup) and two matrices: one is the coverage matrix and the other is the indel matrix, both are table delimited and on a per position basis.
  
  ### How to run it?
  
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
 
</details>
