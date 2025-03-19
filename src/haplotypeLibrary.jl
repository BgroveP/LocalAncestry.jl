function getHaploBlocks(initSize,stepSize,threshold, data, popDict,startFrom=1)
    nCol = size(data,2)
    #NOW IT IS AN ORDERED DICT#
    haploLib = OrderedDict()
    while in(startFrom, 1:nCol)
        if (startFrom == nCol) || (startFrom+stepSize >= nCol) ###use first pop size
            thisBlock = startFrom:nCol #set end of block as chr. end
            ia, haplotypes = computeIA(data[:, thisBlock],popDict) #only to get haplotypes
            haploLib[thisBlock] = haplotypes
            break
        else
            key,value = haploSearch(initSize,stepSize,threshold,data,startFrom, popDict)
            haploLib[key] = value
            startFrom = last(key)+1
        end  
    end
    return haploLib
end
