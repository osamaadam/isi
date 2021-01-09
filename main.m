clc;
clear;

numOfBits = 1000;

randomSeq = randi([0, 1], numOfBits, 1);

L = 10;

hOrg = MultipathChannel(L, numOfBits);

h = tril(toeplitz(hOrg));
hInv = inv(h);

ebOverNo = [-9:0];
energyPerBit = 1;
No = ones(1, length(ebOverNo))./((10.^(ebOverNo/10))*energyPerBit);

noise = (No / 2) .* randn(size(randomSeq));

noise = transpose(mat2cell(noise, [1000], repelem([1], L)));

Y = cellfun(@(x) h * randomSeq + double(x), noise, "UniformOutput", false);

X = cellfun(@(x) hInv * x, Y, "UniformOutput", false);