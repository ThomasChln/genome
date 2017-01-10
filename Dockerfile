from rocker/r-base

run apt-get update
run apt-get install -t unstable -y libcurl4-openssl-dev
run apt-get install -t unstable -y libssl-dev
run R -e "install.packages('devtools')"

run R -e "install.packages('magrittr')"

run R -e "install.packages('data.table')"

run R -e "source('https://bioconductor.org/biocLite.R')"
run R -e "devtools::install_github('genometra/SNPediaR/pkg')"

run apt-get install pandoc -y
run R -e "install.packages('rmarkdown')"
run R -e "install.packages('DT')"
