clc;
clear;

numOfBits = 1000;

randomSeq = randi([0, 1], numOfBits, 1) * 2 - 1;

L = 10;

hOrg = transpose(MultipathChannel(L));

h = tril(toeplitz(hOrg));
hInv = inv(h);

ebOverNo = [-15:0.1:0];
energyPerBit = 1;

for i = 1: length(ebOverNo)
  noise = awgn(randomSeq, ebOverNo(i), energyPerBit);

  noiseCells = mat2cell(reshape(noise, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));

  ##noise = mat2cell(reshape(noise, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));
  signal = mat2cell(reshape(randomSeq, 10, numOfBits / 10), [10], repelem([1], numOfBits / 10));

  Y = cellfun(@(x, n) h * x + n, signal, noiseCells, "UniformOutput", false);
  X = cellfun(@(y) hInv * y, Y, "UniformOutput", false);

  sigMat = reshape(cell2mat(signal), numOfBits, 1);
  xMat = reshape(cell2mat(X), numOfBits, 1);

  meanVal = (max(xMat) + min(xMat)) / 2;

  xMat = xMat >= meanVal;

  ber(i) = (sum(xor((sigMat + 1) / 2, xMat)) / numOfBits) * 100;
end

figure;
plot(ebOverNo, ber);
ylabel("BER");
xlim([-15,0]);