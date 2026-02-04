
using CSV
using DataFrames
using LocalAncestry

chromosome = "22"
output = get_local_ancestries("test/data/current_1kgenomes_2000_1.vcf.gz", "test/data/target_1kgenomes_2000_1.vcf.gz", "test/data/ga_own_1.csv")
