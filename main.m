clc;
clear;

numOfBits = 1000;
numOfRuns = 20;

Fs = 10;

randomSeq = randi([0, 1], numOfBits, 1) * 2 - 1;
repeatedRandomSeq = repelem(randomSeq, Fs, 1);

L = 10;

hOrg = MultipathChannel(L, numOfRuns);

SNR = [-15:0];
energyPerBit = 1;
ber = cell(numOfRuns, 1);
repBer = cell(numOfRuns, 1);

for j = 1: numOfRuns
  h = tril(toeplitz(transpose(hOrg(:, j))));
  hInv = inv(h);
  for i = 1: length(SNR)
    noise = awgn(randomSeq, SNR(i), energyPerBit);
    repeatedNoise = awgn(repeatedRandomSeq, SNR(i), energyPerBit);

    noiseCells = mat2cell(reshape(noise, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));
    signal = mat2cell(reshape(randomSeq, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));
    repeatedSig = mat2cell(reshape(repeatedRandomSeq, 10, (numOfBits * Fs) / 10), [10], repelem([1], (numOfBits * Fs) / 10));
    repeatedNoiseCells = mat2cell(reshape(repeatedNoise, 10, (numOfBits * Fs) / 10), [10], repelem([1], (numOfBits * Fs) / 10));

    Y = cellfun(@(x, n) h * x + n, signal, noiseCells, "UniformOutput", false);
    repeatedY = cellfun(@(x, n) h * x + n, repeatedSig, repeatedNoiseCells, "UniformOutput", false);
    X = cellfun(@(y) hInv * y, Y, "UniformOutput", false);
    repeatedX = cellfun(@(y) hInv * y, repeatedY, "UniformOutput", false);
    
    sigMat = reshape(cell2mat(signal), numOfBits, 1);
    xMat = reshape(cell2mat(X), numOfBits, 1);
    unrepeatedX = transpose(passThroughMatchedFilter(reshape(cell2mat(repeatedX), numOfBits * Fs, 1), Fs));

    meanVal = (max(xMat) + min(xMat)) / 2;

    xMat = xMat >= meanVal;

    ber{j}(i) = sum(xor((sigMat + 1) / 2, xMat));
    repBer{j}(i) = sum(xor((sigMat + 1) / 2, unrepeatedX));
  end
end

avgBer = zeros(1, length(SNR));
avgRepBer = zeros(1, length(SNR));

for i = 1: numOfRuns
  avgBer = avgBer + ber{i};
  avgRepBer = avgRepBer + repBer{i};
end

avgBer = ((avgBer / numOfRuns) / numOfBits) * 100;
avgRepBer = ((avgRepBer / numOfRuns) / numOfBits) * 100;

figure;
plot(SNR, avgBer);
title("Zero Forcing Method")
xlabel("SNR");
ylabel("BER");
xlim([-15,0]);

figure;
plot(SNR, avgRepBer);
title("Matched Filter Method")
xlabel("SNR");
ylabel("BER");
xlim([-15,0]);