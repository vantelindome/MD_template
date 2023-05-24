# Writing publishment section

## Requirements

- Graphs of output from 6analyse
- `pandoc >= 2.12`, `pandoc-crossref`(Use a version for `pandoc`), and `xelatex`.
  - `.github/workflows/analyse.yml` shows how to constract build environment on Ubuntu. The method might include installation of extra packages. If you want to try your machine, please be careful your diskspace.

## Usage

Run `make` on local.

Github Actions automatically build your writing by `.github/workflows/analyse.yml`. You can check outputs at `output` branch. Please make sure that build will occur if you push commits to `main` branch.

## Futher info

`makefile` and `.github/workflows/analyse.yml` are all. There are no magic. You can reproduce this work.
