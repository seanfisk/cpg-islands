# to import this environment into your shell:
#
#     source environment.sh
#

export PYTHONPATH="$PWD"

# fun aliases

alias whatamidoinghere='git status'
# generate tags for Emacs
alias gentags='ctags -eR cpg_islands tests setup.py test.py'
alias t='python test.py'
