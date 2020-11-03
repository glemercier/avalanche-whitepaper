# Avalanche Whitepaper
Avalanche Whitepaper in Tex format for translation to other languages

## Paper generation

### Dependencies

This project has been written in LaTex and requires the `pdflatex` command to build.

Install it on your machine as follows depending on your operating system

**macOS**

Install [homebrew](https://brew.sh/), then install the `mactex` tool:

```
brew install mactex
```

**Ubuntu/debian**

```
sudo apt-get install texlive-latex-base
```

### Run the generation

The project comes with the following make targets:
 - **paper.pdf**: English version of the paper
 - **paper_fr.pdf**: French version of the paper
 - **all**: builds all of the above targets
 
For instance, to generate the French version of the paper, type:
 
```
make paper_fr.pdf
```
 
