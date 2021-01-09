clc;
clear;

numOfBits = 1000;
numOfRuns = 100;

randomSeq = randi([0, 1], numOfBits, 1) * 2 - 1;

L = 10;

hOrg = MultipathChannel(L, numOfRuns);

SNR = [-15:0];
energyPerBit = 1;
ber = cell(numOfRuns, 1);

for j = 1: numOfRuns
  h = tril(toeplitz(transpose(hOrg(:, j))));
  hInv = inv(h);
  for i = 1: length(SNR)
    noise = awgn(randomSeq, SNR(i), energyPerBit);

    noiseCells = mat2cell(reshape(noise, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));
    signal = mat2cell(reshape(randomSeq, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));

    Y = cellfun(@(x, n) h * x + n, signal, noiseCells, "UniformOutput", false);
    X = cellfun(@(y) hInv * y, Y, "UniformOutput", false);

    sigMat = reshape(cell2mat(signal), numOfBits, 1);
    xMat = reshape(cell2mat(X), numOfBits, 1);

    meanVal = (max(xMat) + min(xMat)) / 2;

    xMat = xMat >= meanVal;

    ber{j}(i) = sum(xor((sigMat + 1) / 2, xMat));
  end
end

avgBer = zeros(1, length(SNR));

for i = 1: numOfRuns
  avgBer = avgBer + ber{i};
end

avgBer = ((avgBer / numOfRuns) / numOfBits) * 100;

figure;
plot(SNR, avgBer);
xlabel("SNR");
ylabel("BER");
xlim([-15,0]);