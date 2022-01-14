# The Feltor guide book project

The user guid for the feltor project.

## Install
In order to generate the static website we use
[jupyter-book](https://jupyterbook.org).

```bash
# install jupyter-book
pip install -U jupyter-book
# clone this repository
git clone https://github.com/feltor-dev/user-guide
# build the book
jupyter-book build path/to/user-guide
# we use ghp-import to publish changes on the github pages
pip install ghp-import
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
To publish changes after the book was built
```bash
cd path/to/advection
ghp-import -n -f -p -o _build/html
```

## Author
Matthias Wiesenberger
