function getVCFrows(file::String, c::Union{String,Int})
    # Initialize
    record = VCF.Record()
    loci::Int = 0
    foundFocalChromosome::Bool = false
    this_chromosome = [string(c), "chr" * string(c)]
    reader = VCF.Reader(open(file, "r"))

    while !eof(reader)
        try
            read!(reader, record)
        catch e
            break
        end
        recordChromosome = VCF.chrom(record)

        if any(recordChromosome .== this_chromosome)
            loci += 1
            if !foundFocalChromosome
                foundFocalChromosome = true
            end
        else
            if foundFocalChromosome
                break
            end
        end
    end

    close(reader)

    return loci
end

function getVCFrows2(file::String, c::Union{String,Int})
    # Initialize
    loci::Int = 0
    foundFocalChromosome::Bool = false
    this_chromosome = [string(c), "chr" * string(c)]
    infile = open(file)

    while !eof(infile)

        if any(readuntil(infile, '\t') .== this_chromosome)
            loci += 1
            if !foundFocalChromosome
                foundFocalChromosome = true
            end
        else
            if foundFocalChromosome
                break
            end
        end
        _ = readline(infile)
    end
    close(infile)

    return loci
end

function getVCFdata(f, c, l, n)

    # Initialize
    x = zeros(Int8, 2 * n, l)
    this_chromosome = [string(c), "chr" * string(c)]
    foundFocalChromosome = false
    record = VCF.Record()
    xRows = sort([collect(3:3:(3*n)); collect(1:3:(3*n))])
    i = 0
    dict_record = Dict{UInt8,UInt8}(0x30 => 0, 0x31 => 1, 0x7c => 0)
    reader = VCF.Reader(open(f, "r"))

    while !eof(reader)
        try
            read!(reader, record)
        catch e
            break
        end
        recordChromosome = VCF.chrom(record)

        if any(recordChromosome .== this_chromosome)
            i += 1
            x[:, i] = Int8.(
                get.([dict_record], record.data[vvToV(record.genotype)][xRows], 2)
            )
            if !foundFocalChromosome
                foundFocalChromosome = true
            end
        else
            if foundFocalChromosome
                break
            end
        end
    end
    close(reader)
    return x
end

function getVCFformat!(s)
    spos = position(s)
    for _ in 1:7
        _ = readuntil(s, '\t')
    end
    formatstring = split(readuntil(s, '\t'), ":")
    formatN = length(formatstring)
    formatGT = formatN > 1 ? findfirst(formatstring .== "GT") : 0
    seek(s, spos)

    return (formatGT, formatN)
end
function getVCFdata2(f::String, c::Union{Int,String}, l::Int, n::Int)
    
    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])

    # Iterables
    x = zeros(Int8, 2 * n, l) 

    # Dictionaries
    dict_record = Dict{UInt8,Int8}(0x30 => 0, 0x31 => 1)

    # Start file
    infile = open(f, "r")

    _ = readuntil(infile, '\t')
    _ = readuntil(infile, '\n')

    # Find chromosome
    while !any(readuntil(infile, '\t') .== this_chromosome)
        _ = readuntil(infile, '\n')
    end

    while j <= l
    
    end
end
function skipnchars(infile::IOStream, x::Int)

    for _ in 1:x
        skipchars(isntspace, infile)
        skipchars(isspace, infile)
    end
    return nothing
end

function isntspace(x::Char)
    return !isspace(x)
end

function getVCFdata2(f::String, c::Union{Int,String}, l::Int, n::Int; maxchunksize::Int=200000)

    # Constants
    thischunksize::Int = 0

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    nsave::Int = 0
    endoffile::Bool = false

    # Iterables
    x = zeros(Int8, 2 * n, l) .- 1
    saved::UInt8 = 0x30
    i::Int = 1
    j::Int = 1
    k1::Int = 1
    k2::Int = 1

    # Dictionaries
    const uint8bar = 0x7c
    const uint8tab = 0x09
    const uint8newline = 0x0a
    dict_record = Dict{UInt8,Int8}(0x30 => 0, 0x31 => 1)
    dict_lookup(dict) = key -> dict[key]
    lookup = dict_lookup(dict_record)

    # Start file
    infile = open(f, "r")
    thischunk = Vector{UInt8}(undef, maxchunksize)
    _ = readuntil(infile, '\t')
    _ = readuntil(infile, '\n')

    # Find chromosome
    while !any(readuntil(infile, '\t') .== this_chromosome)
        _ = readuntil(infile, '\n')
    end

    # Loop through the file
    while !endoffile

        thischunksize = readbytes!(infile, thischunk, maxchunksize)
        endoffile = thischunksize < (maxchunksize)
        thisloc = findall(thischunk .== uint8bar)
        thisnewline = findall(thischunk .== uint8newline)
        thistab = findall(thischunk .== uint8bar)


        # If no relevant data in chunk, skip loop or both skip loop and increment locus counter
        if (length(thisloc) == 0) && (length(thisnewline) == 0)
            continue
        elseif (length(thisloc) == 0) && (length(thisnewline) > 0)
            j = j + 1
            continue
        end

        k1 = 1
        k2 = 0
        k3 = 
        if length(thisnewline) > 0
            for nl in thisnewline
                if j > l
                    break
                end
                k1 = k2 + 1
                k2 = searchsortedfirst(thisloc[k1:end], nl) - 2 + k1
                # Before
                if thisloc[k1] != 1
                    x[i:2:(2*(k2-k1+1)+i-1), j] = lookup.(thischunk[thisloc[k1:k2].-1])
                else
                    nothing
                end

                # After
                if thisloc[k2] != maxchunksize
                    x[(i+1):2:(2*(k2-k1+1)+i), j] = lookup.(thischunk[thisloc[k1:k2].+1])
                end
                j = j + 1
            end

        end

        # After last newline
        k1 = k2 + 1
        x[i:2:(length(thisloc)-k1+i), j] = lookup.(thischunk[thisloc[k1:2:end].-1])
        x[(i+1):2:(length(thisloc)-k1+i), j] = lookup.(thischunk[thisloc[k1:2:end].+1])
    end
    close(infile)

    return x
end
using BenchmarkTools
@benchmark getVCFdata2(f, c, l, n)
@benchmark LocalAncestry.getVCFdata(f, c, l, n)

function getVCFindividuals(f)::Vector{String}
    r = VCF.Reader(open(f, "r"))
    x = r.header.sampleID
    close(r)
    return x
end

function readVCF(file::String, chromosome::Int64)

    # Checks
    assertFile(file, ".vcf")
    assertPositiveInteger(
        chromosome;
        errormessage="chromosome was $(chromosome), but should be an integer larger than zero",
    )

    # Get the number of loci from the focal chromosome
    individuals = getVCFindividuals(file)
    loci = getVCFrows(file, chromosome)
    x = getVCFdata(file, chromosome, loci, length(individuals))

    return x, individuals
end

function readVCF2(file::String, chromosome::Int64)

    # Checks
    assertFile(file, ".vcf")
    assertPositiveInteger(
        chromosome;
        errormessage="chromosome was $(chromosome), but should be an integer larger than zero",
    )

    # Get the number of loci from the focal chromosome
    individuals = getVCFindividuals(file)
    loci = getVCFrows2(file, chromosome)
    x = getVCFdata(file, chromosome, loci, length(individuals))

    return x, individuals
end

function readTrue(file::String, map::String, chromosome::Int64, individuals::Vector{String})
    genomic_map = CSV.read(map, DataFrame)
    columns_from_file = findall(genomic_map.chromosome .== chromosome)
    out = zeros(Int32, 2 * length(individuals), length(columns_from_file))
    outInd = string.(zeros(Int8, 2 * length(individuals)))
    inviduals_strings = string.(individuals)
    open(file, "r") do reader
        iterator = 1
        for line in eachline(reader)
            # Convert the line to an integer
            splitted_string = split(line, "\t")
            numbers = parse.(Int, splitted_string[2:end])
            if any(splitted_string[1] .== inviduals_strings)
                out[iterator, :] = numbers[columns_from_file]
                outInd[iterator] = splitted_string[1]
                iterator = iterator + 1
            end
        end
    end

    return out, outInd
end

function readRfmix(path::String, mappath::String)
    if !isfile(path)
        error("Assignment file not found!")
    end

    if !isfile(mappath)
        error("Map file not found!")
    end

    # Read the data
    data = CSV.read(path, DataFrame; header=2)
    genomic_map = CSV.read(mappath, DataFrame)

    individuals = unique(replace.(names(data)[5:end], r":::.*" => ""))
    populations = unique(replace.(names(data)[5:end], r".*:::hap\d+:::" => ""))
    haplotypes = unique(
        replace.(replace.(names(data)[5:end], r".*:::hap" => ""), r":::.*" => "")
    )
    chromosome = data.chromosome[1]
    nMarkers = sum(genomic_map.chromosome .== chromosome)
    nBlocks = size(data, 1)
    blocks = Vector{UnitRange{Int64}}(undef, size(data, 1))
    for (i, j) in enumerate(data.genetic_marker_index)
        blockstart = j + 1
        blockend = i == nBlocks ? nMarkers : data.genetic_marker_index[i+1]
        blocks[i] = blockstart:blockend
    end

    # Fill library
    library = OrderedDict{UnitRange{Int64},Vector{Matrix{Int8}}}()
    for i in 1:nBlocks
        library[blocks[i]] = [zeros(Int8, 1, 1)]
    end

    # Fill class
    class = OrderedDict{String,Vector{Union{Missing,String}}}()
    levels1 = [i * "_hap" * h for i in individuals for h in haplotypes]
    levels2 = [i * ":::hap" * h * ":::" for i in individuals for h in haplotypes]

    for (i, ind) in enumerate(levels1)
        probcolumns = levels2[i] .* populations
        class[ind] = [
            populations[argmax([data[j, c] for c in probcolumns])] for j in 1:nBlocks
        ]
    end

    # Fill probabilities
    probs = OrderedDict{String,OrderedDict{String,Vector{Union{Missing,Float64}}}}()
    for (i, ind) in enumerate(levels1)
        probs[ind] = OrderedDict{String,Vector{Union{Missing,Float64}}}()
        for p in populations
            probs[ind][p] = data[:, levels2[i].*p]
        end
    end

    # Return
    return probs, class, library
end

