function accumulateRD = chan_Accumulate(fft2dDataDB)
    [channelNum,rangeBin,dopplerBin] = size(fft2dDataDB);
    accumulateRD = zeros(rangeBin,dopplerBin);

    for channelId = 1:channelNum
        accumulateRD = accumulateRD + (abs(squeeze(fft2dDataDB(channelId,:,:))));
    end
end
