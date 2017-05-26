# Bash script used to update the docs page
# Commented section is meant for setup (run once only)

mkdir -p docs-sphinx
mkdir -p docs

## Set up sphinx metadata in docs-sphinx folder
#sphinx-apidoc -F -o docs-sphinx pyprot/


# Copy configuration script (tracked by git) to docs-sphinx
cp docs-config.py docs-sphinx/conf.py

# Make html website in docs-sphinx
cd docs-sphinx && make html

# Copy files in docs
cd ../ && cp -R docs-sphinx/_build/html/. docs/

# Create .nojekyll file in docs/ (allows the use of folders with _ prefix)
touch ./docs/.nojekyll