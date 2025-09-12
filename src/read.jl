
function readVCF(path::AbstractString, c::Union{AbstractString,Int})

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = LocalAncestry.string2UInt8.(this_chromosome)
    chromosome_char_lengths = length.(chromosome)

    loci = LocalAncestry.getVCFloci(path, c)
    individuals = LocalAncestry.getVCFindividuals(path)
    x = zeros(Int8, 2 * length(individuals), length(loci))
    file = LocalAncestry.vcf_open(path)

    locusdict = Dict{String,Int}(loci .=> 1:length(loci))

    # Read
    # Initialize
    ## Containers
    buffer = zeros(UInt8, READLINE_BUFFER_SIZE)
    tmplocus = zeros(Int8, 2 * length(individuals))
    chromosome_match = [true, true]

    # Integers for iteration
    pos::Int = 1
    pos1::Int = 1
    pos2::Int = 1
    pos3::Int = 1

    # Read file until the end of file
    line = 0
    while LocalAncestry.readline!(file, buffer)
        line = line + 1
        if all(buffer .== UInt8('a'))
            break
        end
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
                LocalAncestry.vcfrow_haplotype(buffer, tmplocus, pos1, pos2, pos3)
                x[:, locusdict[(LocalAncestry.vcfrow_locus(buffer, pos1, pos2))]] = tmplocus
            end
        end
        buffer .= UInt8('a')
    end

    # Return
    return x, loci, string.(individuals)
end



function getVCFindividuals(path::AbstractString)

    # Initialize
    file = vcf_open(path)
    bufvec = Vector{UInt8}(undef, READLINE_BUFFER_SIZE)

    # Get loci
    while readline!(file, bufvec)
        if bufvec[2] != UInt8('#')
            return replace.(split(join(Char.(bufvec[1:(findfirst(bufvec .== UInt8('\n'))-1)]), ""), '\t')[10:end], ['\r'] .=> "")
        end
    end
    close(file)
    return [""]
end

function readline!(s::GZip.GZipStream, b::Vector{UInt8})

    # Initialize 
    p::Union{Int,Nothing} = 1

    # Read until we either hit the end of the file, or a newline
    while !eof(s)
        gzgets(s, pointer(b) + p - 1, length(b) - p + 1) # Reads into buffer
        p = findnext(x -> x == UInt8('\0'), b, p) # Finds the \0 symbol which marks the end of the input read above
        if ~isnothing(p) && (b[p-1] == UInt8('\n')) # Did we find a newline?
            return true # yes: stop here with success
        else
            resize!(b, length(b) + READLINE_BUFFER_SIZE) # no: grow the vector and go again
        end
    end
    # Reached EOF without finding a newline

    return false
end

function readline!(s::IOStream, b::Vector{UInt8})
    d::UInt8 = UInt8('\n')
    p::Int = 1

    while !eof(s)
        Base.@_lock_ios s n = Int(ccall(:jl_readuntil_buf, Csize_t, (Ptr{Cvoid}, UInt8, Ptr{UInt8}, Csize_t), s.ios, d, pointer(b, p), (length(b) - p + 1) % Csize_t))
        p += n
        if b[p-1] == d
            return true # yes: stop here with success
        else
            resize!(b, length(b) + READLINE_BUFFER_SIZE) # no: grow the vector and go again
        end

    end

    return false
end

function getVCFloci(path::AbstractString, c::Union{AbstractString,Int})

    # Initialize
    loci = Vector{String}()
    chromosome = ["", "chr"] .* string.([c, c])
    file = vcf_open(path)
    bufvec = Vector{UInt8}(undef, READLINE_BUFFER_SIZE)

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
    close(file)

    return loci
end
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

function skipline!(s::IOStream, b::Vector{UInt8})
    _ = readline(s)
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



function read_reference_ancestries(ancestries::Union{DataFrame,String})

    if isa(ancestries, String)
        ancestries = CSV.read(ancestries, DataFrame)
    end

    return standardize_ancestries(ancestries)
end

function standardize_ancestries(ancestries::DataFrame)

    # Check for mandatory headers
    mandatory_headers = ["individual", "population"]
    for m in mandatory_headers
        if ~any(m .== names(ancestries))
            error("The reference ancestries don't contain the mandatory header: " * m)
        end

        if typeof(ancestries[:, m]) != Vector{String}
            ancestries[:, m] = string.(ancestries[:, m])
        end
    end

    # Check or fill optional columns
    optional_columns = ["omit"]
    for o in optional_columns
        if ~any(o .== names(ancestries))
            # If not present, insert standard value
            ancestries[:, o] .= 0
        end
    end

    return ancestries
end


function readVCF_reference(path::AbstractString, c::Union{AbstractString,Int}, ancestries::DataFrame)

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = LocalAncestry.string2UInt8.(this_chromosome)
    chromosome_char_lengths = length.(chromosome)

    # Get focal individuals and the popdict
    individuals = LocalAncestry.getVCFindividuals(path)
    vcfcols, xrows, popDict = get_popdict(individuals, ancestries)
    loci = LocalAncestry.getVCFloci(path, c)


    x = zeros(Int8, 2 * length(individuals), length(loci))
    file = LocalAncestry.vcf_open(path)

    locusdict = Dict{String,Int}(loci .=> 1:length(loci))

    # Read
    # Initialize
    ## Containers
    buffer = zeros(UInt8, READLINE_BUFFER_SIZE)
    tmplocus = zeros(Int8, 2 * length(individuals))
    chromosome_match = [true, true]

    # Integers for iteration
    pos::Int = 1
    pos1::Int = 1
    pos2::Int = 1
    pos3::Int = 1

    # Read file until the end of file
    line = 0
    while LocalAncestry.readline!(file, buffer)
        line = line + 1
        if all(buffer .== UInt8('a'))
            break
        end
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
                LocalAncestry.vcfrow_haplotype(buffer, tmplocus, pos1, pos2, pos3)
                x[:, locusdict[(LocalAncestry.vcfrow_locus(buffer, pos1, pos2))]] = tmplocus
            end
        end
        buffer .= UInt8('a')
    end

    # Return
    return x, loci, string.(individuals)
end

function getVCFloci2(path::AbstractString, c::Union{AbstractString,Int}, popDict, vcfcols, xrows)

    # Initialize
    this_chromosome = ["", "chr"] .* string.([c, c])
    chromosome = LocalAncestry.string2UInt8.(this_chromosome)
    chromosome_char_lengths = length.(chromosome)
    file = LocalAncestry.vcf_open(path)

    # Read
    # Initialize
    ## Containers
    buffer = zeros(UInt8, READLINE_BUFFER_SIZE)
    tmplocus = zeros(Int8, 2 * length(individuals))
    chromosome_match = [true, true]
    MAF = Vector{Float64}()
    loci = Vector{String}()

    # Integers for iteration
    pos::Int = 1
    pos1::Int = 1
    pos2::Int = 1
    pos3::Int = 1

    # Read file until the end of file
    p = zeros(Float64, length(keys(popDict)))
    while LocalAncestry.readline!(file, buffer)
        if all(buffer .== UInt8('a'))
            break
        end
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
                LocalAncestry.vcfrow_haplotype(buffer, tmplocus, pos1, pos2, pos3)
                push!(loci, LocalAncestry.vcfrow_locus(buffer, pos1, pos2))
                x = @view tmplocus[vcfcols][xrows]
                push!(MAF, ponelocus(x, popDict, p))
            end
        end
        buffer .= UInt8('a')
    end
    close(file)
    return loci, MAF
end

function ponelocus(x, popDict, p)
    for (popi, pop) in enumerate(keys(popDict))
        p[popi] = LocalAncestry.mean(x[popDict[pop]]) + NEARZERO_FLOAT
    end
    p .= min.(p, 1 .- p)
    return maximum(p)
end
