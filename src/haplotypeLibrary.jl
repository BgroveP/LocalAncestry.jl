function getHaploBlocks(initSize,stepSize,threshold, data, popDict,startFrom=1)
    nCol = size(data,2)
    haploLib = OrderedDict()
    while in(startFrom, 1:nCol)
        if (startFrom == nCol) || (startFrom+stepSize >= nCol) 
            thisBlock = startFrom:nCol 
            ia, haplotypes = ARV.computeIA(data[:, thisBlock],popDict) 
            haploLib[thisBlock] = haplotypes
            break
        else
            key,value = ARV.haploSearch(initSize,stepSize,threshold,data,startFrom, popDict)
            haploLib[key] = value
            startFrom = last(key)+1
        end  
    end
    return haploLib, length(haploLib)
end
