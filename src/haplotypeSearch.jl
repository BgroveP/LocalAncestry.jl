function haploSearch(initSize, stepSize, threshold, data, startPos, popDict)
    nCol = size(data, 2)
    lastSetChrEnd = false
    thisBlock = startPos:(startPos+initSize-1)
    currentIA, currentHaplo = computeIA(data[:, thisBlock], popDict)
    while true
        if last(thisBlock) == nCol
            break #ensure stopping of continuation if already came to an end
        end
        thisBlock = rangeChange(thisBlock; lastInc=stepSize)
        if last(thisBlock) > nCol ###use first pop size
            origLast = last(thisBlock)
            lastSetChrEnd = true
            thisBlock = first(thisBlock):nCol #set end of block as chr. end
        end
        #global is needed for printing in line 26. No other affect.
        tempIA, tempHaplo = computeIA(data[:, thisBlock], popDict)
        stat = (tempIA - currentIA) / currentIA
        if (stat <= threshold) #|| (currentIA == 0)
            #==if an increase in region size, exceeds chr end, above I set the end to be chr end.
            Here, if this  increase does not increase IA, I get back to original region size (which exceeds chr end),
            then reduce it back to what it was before the increase.
            ==#
            thisBlock = if lastSetChrEnd == true
                rangeChange(first(thisBlock):origLast; lastDec=stepSize)
            else
                rangeChange(thisBlock; lastDec=stepSize)
            end
            tempIA, tempHaplo = computeIA(data[:, thisBlock], popDict) #return to previous haplotype also!
            break
        end
        currentIA, currentHaplo = tempIA, tempHaplo
    end
    return thisBlock, currentHaplo
end

function haploSearch2(initSize, stepSize, threshold, data, startPos, popDict)
    nCol = size(data, 2)
    lastSetChrEnd = false
    thisBlock = startPos:(startPos+initSize-1)
      n = zeros(Int64, length(popDict))
    for (i, k) in enumerate(keys(popDict))
        n[i] = length(popDict[k])
    end
    currentIA, currentHaplo, currentHaplocount = computeIA3(data, thisBlock, popDict, n)
    while true
        if last(thisBlock) == nCol
            break #ensure stopping of continuation if already came to an end
        end
        thisBlock = LocalAncestry.rangeChange(thisBlock; lastInc=stepSize)
        if last(thisBlock) > nCol ###use first pop size
            origLast = last(thisBlock)
            lastSetChrEnd = true
            thisBlock = first(thisBlock):nCol #set end of block as chr. end
        end
        #global is needed for printing in line 26. No other affect.
        tempIA, tempHaplo, tempHaplocount = computeIA3(data, thisBlock, popDict, n)
        stat = (tempIA - currentIA) / currentIA
        if (stat <= threshold) #|| (currentIA == 0)
            #==if an increase in region size, exceeds chr end, above I set the end to be chr end.
            Here, if this  increase does not increase IA, I get back to original region size (which exceeds chr end),
            then reduce it back to what it was before the increase.
            ==#
            thisBlock = if lastSetChrEnd == true
                LocalAncestry.rangeChange(first(thisBlock):origLast; lastDec=stepSize)
            else
                LocalAncestry.rangeChange(thisBlock; lastDec=stepSize)
            end
            tempIA, tempHaplo = computeIA3(data, thisBlock, popDict, n) #return to previous haplotype also!
            break
        end
        currentIA, currentHaplo, currentHaplocount = tempIA, tempHaplo, tempHaplocount
    end
    return thisBlock, currentHaplo, currentHaplocount
end

function haploSearch4(initSize, stepSize, threshold, data, startPos, popDict)
    nCol = size(data, 2)
    lastSetChrEnd = false
    thisBlock = startPos:(startPos+initSize-1)
      n = zeros(Int64, length(popDict))
    for (i, k) in enumerate(keys(popDict))
        n[i] = length(popDict[k])
    end
    currentIA, currentHaplo, currentHaplocount = computeIA5(data, thisBlock, popDict, n)
    while true
        if last(thisBlock) == nCol
            break #ensure stopping of continuation if already came to an end
        end
        thisBlock = LocalAncestry.rangeChange(thisBlock; lastInc=stepSize)
        if last(thisBlock) > nCol ###use first pop size
            origLast = last(thisBlock)
            lastSetChrEnd = true
            thisBlock = first(thisBlock):nCol #set end of block as chr. end
        end
        #global is needed for printing in line 26. No other affect.
        tempIA, tempHaplo, tempHaplocount = computeIA3(data, thisBlock, popDict, n)
        stat = (tempIA - currentIA) / currentIA
        if (stat <= threshold) #|| (currentIA == 0)
            #==if an increase in region size, exceeds chr end, above I set the end to be chr end.
            Here, if this  increase does not increase IA, I get back to original region size (which exceeds chr end),
            then reduce it back to what it was before the increase.
            ==#
            thisBlock = if lastSetChrEnd == true
                LocalAncestry.rangeChange(first(thisBlock):origLast; lastDec=stepSize)
            else
                LocalAncestry.rangeChange(thisBlock; lastDec=stepSize)
            end
            tempIA, tempHaplo = computeIA3(data, thisBlock, popDict, n) #return to previous haplotype also!
            break
        end
        currentIA, currentHaplo, currentHaplocount = tempIA, tempHaplo, tempHaplocount
    end
    return thisBlock, currentHaplo, currentHaplocount
end
