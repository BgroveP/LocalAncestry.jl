function getHaploBlocks(initSize, stepSize, threshold, data, popDict, startFrom=1)
    nCol = size(data, 2)
    haploLib = OrderedDict()
    while in(startFrom, 1:nCol)
        if (startFrom == nCol) || (startFrom+stepSize+initSize >= nCol)
            thisBlock = startFrom:nCol
            ia, haplotypes = LocalAncestry.computeIA(data[:, thisBlock], popDict)
            haploLib[thisBlock] = haplotypes
            break
        else
            key, value = LocalAncestry.haploSearch(
                initSize, stepSize, threshold, data, startFrom, popDict
            )
            haploLib[key] = value
            startFrom = last(key)+1
        end
    end
    return haploLib, length(haploLib)
end

function getHaploBlocks2(initSize, stepSize, threshold, data, popDict, direction=1)

    # Initialize
    nCol = size(data, 2)
    haploLib = OrderedDict()
    notdone = true
    spots = 1:nCol
    loci = direction == 1 ? (1:nCol) : (nCol:-1:1)

    # Asserts
    assertInRange(
        startFrom,
        loci,
        "While building the haplotype library, the start position not within data.",
    )
    assertInRange(direction, (-1:2:1), "The direction must either be -1 or 1.")

    # Estimate
    while notdone
        if !in(startFrom + stepSize, loci)
            thisBlock = startFrom:nCol
            ia, haplotypes = computeIA(data[:, thisBlock], popDict)
            haploLib[thisBlock] = haplotypes
            notdone = false
            break
        else
            key, value = haploSearch(
                initSize, stepSize, threshold, data, startFrom, popDict
            )
            haploLib[key] = value
            startFrom = last(key)+1
        end
    end

    return haploLib, length(haploLib)
end

