% function compCoffVec = generateCompCoff(parameter,speedBin)
% 
%     compCoffVec = zeros(parameter.virtualAntenna,1);
%     txNum = length(parameter.txAntenna);
%     rxNum = length(parameter.rxAntenna);
%     phi = 2 * pi * (speedBin - parameter.dopplerBin/2 - 1) / parameter.dopplerBin;
%     delta = phi / txNum;
%     for txId = 1:txNum
%         for rxId = 1:rxNum
%             compCoffVec((txId-1) * rxNum + rxId) = exp(-1i * (txId-1) * delta);
%         end
%     end
% end
function compCoffVec = generateCompCoff(parameter,speedBin)

    compCoffVec = zeros(parameter.virtualAntenna,1);
    txNum = length(parameter.txAntenna);
    rxNum = length(parameter.rxAntenna);
    phi = 2 * pi * (speedBin - parameter.dopplerBin/2 - 1) / parameter.dopplerBin;
    for txId = 1:txNum
        for rxId = 1:rxNum
            compCoffVec((txId-1) * rxNum + rxId) = exp(-1i * (txId-1) * phi);
        end
    end
end