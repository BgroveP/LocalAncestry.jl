function open_vcf(path::AbstractString; opentype="r")
    if path[(end-6):end] == ".vcf.gz"
        file = GZip.open(path, opentype)
    elseif path[(end-3):end] == ".vcf"
        file = open(path, opentype)
    else
        error("Input file was not a (compressed) .vcf file")
    end
    return file
end
