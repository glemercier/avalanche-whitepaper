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
sudo apt install texlive-latex-base texlive-science texlive-lang-french texlive-fonts-extra texlive-latex-extra make
```

### Run the generation

The project comes with the following make targets:
 - **paper.pdf**: English version of the paper
 - **paper_fr.pdf**: French version of the paper (default)
 - **all**: builds all of the above targets

For instance, to generate the French version of the paper, type:

```
make
```

## Submissions

Anyone can participate and help translate the whitepaper. Here are the steps to follow to make your own submissions:
 - Go to the [Issues](https://github.com/glemercier/avalanche-whitepaper/issues) and look for a task you are willing to work on
 - Check that it has not already been assigned to someone
 - Add a comment mentioning you volunteer for the task. One of the reviewer will assign the task to you
 - When you have done the translation of the part associated to the issue, submit your work as a [Pull Request](https://github.com/glemercier/avalanche-whitepaper/pulls)
 - The reviewers will check your work and may make some comments for you to modify
 - Once the reviewers are happy with you submission, it will be merged into the source. Then feel free to pick another task!

