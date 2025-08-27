# Read VCF prior
function getVCFrows(file::String, c::Union{String,Int})
    # Initialize
    loci::Int = 0
    foundFocalChromosome::Bool = false
    this_chromosome = [string(c), "chr" * string(c)]
    infile = vcf_open(file)

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
function getVCFindividuals(f)::Vector{String}
    r = VCF.Reader(open(f, "r"))
    x = r.header.sampleID
    close(r)
    return x
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



function readVCF(file::String, chromosome::Int64)

    # Checks
    #assertFile(file, ".vcf")
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

# Read VCF new
function skipline!(s::GZip.GZipStream, buf::Vector{UInt8})
    while !eof(s)
        # Read a chunk
        ptr = gzgets(s, buf)
        ptr == C_NULL && error("Error reading from GZipStream")
        # Find the first newline in the chunk
        idx = findfirst(==(UInt8('\n')), buf)
        if idx !== nothing
            # Found a newline; skip to the next line
            return true
        end
    end
    # Reached EOF without finding a newline
    return false
end

function readline!(s::GZip.GZipStream, buf::Vector{UInt8}, pos::Int; b::Int = 10000)
    # Grow buffer until we find an eol character
    while !eof(s)
        gzgets(s, pointer(buf) + pos - 1, length(buf) - pos + 1)
        pos = findnext(x -> x == UInt8('\0'), buf, pos)::Int
        if buf[pos-1] == UInt8('\n')
            return true
        else
            resize!(buf, length(buf) + b)
        end
    end
    # Reached EOF without finding a newline
    return false
end

function getVCFloci(path::AbstractString, c::Union{AbstractString,Int})

    # Constants 
    buffersize = 1000000

    # Initialize
    loci = Vector{String}()
    chromosome = ["", "chr"] .* string.([c, c])
    file = vcf_open(path)
    bufvec = Vector{UInt8}(undef, buffersize)

    # Get loci
    while !eof(file)
        try
            if any(readuntil(file, '\t') .== chromosome)

                push!(loci, readuntil(file, '\t'))
            end
            skipline!(file, bufvec)
        catch
            if !eof(file)
                error("Encountered unknown error while reading file")
            end
        end
    end

    return loci
end

function getVCFloci2(path::AbstractString, c::Union{AbstractString,Int})

    # Constants 
    buffersize = 10000

    # Initialize
    loci = Vector{String}()
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = string2UInt8.(this_chromosome)

    file = vcf_open(path)
    bufvec = Vector{UInt8}(undef, buffersize)
    pos1::Int = 1
    pos2::Int = 1
    
    # Scroll past headers
    while readline!(file,bufvec,pos1)
        if bufvec[2] != UInt8('#')
            break
        end
    end
    
    # Read data lines
    while readline!(file, bufvec, pos1)
        if vcf_row_chromosome_check(bufvec, chromosome)
            push!(loci, vcfrow_locus(bufvec, pos1, pos2))
        end
    end

    return loci
end

function vcf_row_chromosome_check(bufvec::Vector{UInt8}, chromosome::Vector{Vector{UInt8}})
    return any(all(bufvec[1:length(chromosome[1])] .== chromosome[1]) || all(bufvec[1:length(chromosome[2])] .== chromosome[2]))
end
function getVCFindividuals2(path::AbstractString)

    # Constants 
    initialbuffersize = 100000

    # Initialize
    file = vcf_open(path)
    bufvec = Vector{UInt8}(undef, initialbuffersize)

    # Get loci
    while !eof(file)
        skip(file, 1)
        if Char(peek(file)) != '#'
            return split(readline(file)[1:(end-1)], '\t')[10:end]
        end
        skipline!(file, bufvec)

    end
    return [""]
end

function vcfrow_locus(b::Vector{UInt8}, pos1::Int, pos2::Int)

    pos1 = 1
    pos2 = 1
    # Find first separator
    while b[pos1] != UInt8('\t')
        pos1 = pos1 + 1
    end

    # Find second separator
    pos1 = pos1 + 1
    pos2 = pos1
    while b[pos2+1] != UInt8('\t')
        pos2 = pos2 + 1
    end

    return join(Char.(b[pos1:pos2]), "")
end

function vcfrow_haplotype(b::Vector{UInt8}, locusout::Vector{Int8}, pos1::Int, pos2::Int, pos3::Int)

    pos1 = 0 # Used as '\t' counter
    pos2 = 1 # Used for marking position
    pos3 = 1
    # Find first separator
    while pos1 < 8
        if b[pos2] == UInt8('\t')
            pos1 = pos1 + 1
        end
        pos2 = pos2 + 1
    end

    # Find second separator
    pos1 = pos2 # Now also used for marking position
    while b[pos1+1] != UInt8('\t')
        pos1 = pos1 + 1
    end

    # Get haplotypes
    if (b[pos2] == UInt8('G')) & (b[pos1] == UInt8('T')) & ((pos1 - pos2) == 1)
        # If only field is "GT"
        while (pos3 < length(locusout)) & (pos1 < length(b))
            if b[pos1+1] == UInt8('|')
                locusout[pos3] = b[pos1] - 48
                locusout[pos3+1] = b[pos1+2] - 48
                pos3 = pos3 + 2
            end
            pos1 = pos1 + 1
        end
    else
        # If there are multiple fields
        error("Can't read vcfs with multiple fields yet")
    end

end

function vcf_open(path::AbstractString; opentype="r")
    if path[(end-6):end] == ".vcf.gz"
        file = GZip.open(path, opentype)
    elseif path[(end-3):end] == ".vcf"
        file = open(path, opentype)
    else
        error("Input file was not a (compressed) .vcf file")
    end

    return file
end

function readVCF2(path::AbstractString, c::Union{AbstractString,Int})

    # Constants 
    buffersize = 10000

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = string2UInt8.(this_chromosome)
    chromosome_char_lengths = length.(chromosome)

    loci = getVCFloci(path, c)
    individuals = getVCFindividuals2(path)
    x = zeros(Int8, 2 * length(individuals), length(loci))
    file = vcf_open(path)

    locusdict = Dict{String,Int}(loci .=> 1:length(loci))

    # Read
    # Initialize
    ## Containers
    buffer = zeros(UInt8, buffersize)
    tmplocus = zeros(Int8, 2 * length(individuals))
    chromosome_match = [true, true]

    # Integers for iteration
    pos1::Int = 1
    pos2::Int = 1
    pos3::Int = 1

    # Read file until the end of file
    while readline!(file, buffer, pos1)

        # Skip all headers
        if buffer[1] != UInt8('#')

            # Check if the line contains the focal chromosome
            chromosome_match .= true
            for k in eachindex(chromosome_match)
                pos = 1
                while pos <= chromosome_char_lengths[k]
                    if buffer[pos] != chromosome[k][pos]
                        chromosome_match[k] = false
                        break
                    else
                        pos = pos + 1
                    end
                end
            end

            # If line contains focal chromosome, read haplotypes
            if any(chromosome_match)
                vcfrow_haplotype(buffer, tmplocus, pos1, pos2, pos3)
                x[:, locusdict[(vcfrow_locus(buffer, pos1, pos2))]] = tmplocus
            end
        end
    end

    # Return
    return x, loci, individuals
end

function readVCF3(path::AbstractString, c::Union{AbstractString,Int})

    # Constants 
    buffersize = 10000

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = string2UInt8.(this_chromosome)
    chromosome_char_lengths = length.(chromosome)

    loci = getVCFloci(path, c)
    individuals = getVCFindividuals2(path)
    x = zeros(Int8, 2 * length(individuals), length(loci))
    file = vcf_open(path)

    locusdict = Dict{String,Int}(loci .=> 1:length(loci))

    # Read
    # Initialize
    ## Containers
    buffer = zeros(UInt8, buffersize)
    tmplocus = zeros(Int8, 2 * length(individuals))
    chromosome_match = [true, true]

    # Integers for iteration
    pos1::Int = 1
    pos2::Int = 1
    pos3::Int = 1

    # Read file until the end of file
    while readline!(file, buffer, pos1)

        # Skip all headers
        if buffer[1] != UInt8('#')

            # Check if the line contains the focal chromosome
            chromosome_match .= true
            for k in eachindex(chromosome_match)
                pos = 1
                while pos <= chromosome_char_lengths[k]
                    if buffer[pos] != chromosome[k][pos]
                        chromosome_match[k] = false
                        break
                    else
                        pos = pos + 1
                    end
                end
            end

            # If line contains focal chromosome, read haplotypes
            if any(chromosome_match)
                vcfrow_haplotype(buffer, tmplocus, pos1, pos2, pos3)
                x[:, locusdict[(vcfrow_locus(buffer, pos1, pos2))]] = tmplocus
            end
        end
    end

    # Return
    return x, loci, individuals
end
