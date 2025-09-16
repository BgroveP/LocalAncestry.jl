
function samples(path)
    # Open file and initialize string
    file = open_vcf(path)
    s::String = "1234"
    while ~eof(file)
        s = readline(file)
        if s[1:7] == "#CHROM\t"
            close(file)
            return string.(split(replace(s, '\n' => ""), '\t')[10:end])
        end
    end
    close(file)
    error("Reached end of file without finding samples")
end
