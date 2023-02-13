# FASE

Analysis of RNA-Sequencing data using FASE (Finding Alternative Splicing Events).

The pipeline is based on differential alternative splicing events and predicts the transcript structure and their concentration along with survival analysis. This is the first kind of pipeline that takes advantage of differential alternative splicing events for finding novel transcripts that are neglected due to low expression in transcript level statistical analysis. The pipeline generates reproducible results.

It has seven modules for downstream analysis of RNA-Sequencing data. These include: differential alternative splicing, transcript structure, transcript concentration, survival analysis, network analysis, differential gene expression, and differential junction expression.

## Installation:
### Using R
download.file(url = 'https://github.com/harshsharma-cb/FASE/archive/refs/heads/main.zip', destfile = 'FASE-main.zip')

install.packages('FASE-main.zip', type = 'source', repos = NULL)

### Using devtools in R
devtools::install_github('FASE')

## Citation:
Sharma, H., Pani, T., Dasgupta, U., Batra, J., and Sharma, R.D., Prediction of transcript structure and concentration using RNA-Seq data, Briefings in Bioinformatics, 2023;, bbad022, https://doi.org/10.1093/bib/bbad022
