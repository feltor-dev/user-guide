# The Feltor guide book project

The user guide for the feltor project.

## Install
In order to generate the static website we use
[jupyter-book](https://jupyterbook.org).

```bash
# install jupyter-book
conda install jupyter-book -c conda-forge
# clone this repository
git clone https://github.com/feltor-dev/user-guide
# build the book
jupyter-book build path/to/user-guide
```
In order to locally generate the simulation data you will need the
[Feltor](https://github.com/feltor-dev/feltor) code repository.  Follow the
quick-start guide to install.  It is recommended to keep Feltor and this
repository next to each other.  If you prefer not to, you need to set the
`FELTOR_PATH` in all jupyter notebooks manually.

## Usage
Build the book with
```bash
jupyter-book build path/to/user-guide
```
To publish changes after the book was built we use the
`ghp-import` python package (using the `-o` option
to only keep a single commit on the `gh-pages` branch)
```bash
pip install ghp-import
cd path/to/user-guide
ghp-import -n -f -p -o _build/html
```
## Last succesful build
```bash
Jupyter Book      : 1.0.0
External ToC      : 1.0.1
MyST-Parser       : 2.0.0
MyST-NB           : 1.0.0
Sphinx Book Theme : 1.1.0
Jupyter-Cache     : 0.6.1
NbClient          : 0.5.13
```
## Author
Matthias Wiesenberger
