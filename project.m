clc
clear
close all

% Parameters
snr = 5;
snrShot = -50;
pShot = 0.00;
shotLength = 5;

% Load and Display the Image
load mig25.mat;
migData = X;
migImage = abs(fftshift(fft2(migData)));
migImage = floor(migImage/max(max(migImage))*255);

migImageBits = dec2bin(migImage);
migImageBits = reshape(migImageBits.', 1, []);
migImageBits = split(migImageBits, "");
migImageBits = cell2mat(migImageBits);
migImageBits = str2num(migImageBits).';

showImage(migImageBits, "Original Image");

% Simulate transmission with different settings
recv = simulateModel(migImageBits, snr, snrShot, pShot, shotLength, false, false);
recvHamming = simulateModel(migImageBits, snr, snrShot, pShot, shotLength, true, false);
recvInterlaced = simulateModel(migImageBits, snr, snrShot, pShot, shotLength, false, true);
recvHammingInterlaced = simulateModel(migImageBits, snr, snrShot, pShot, shotLength, true, true);

% Display results
showImage(recv, "Transmitted Image. BER = " + num2str(calcBer(migImageBits, recv)));
showImage(recvHamming, "Transmitted Image with Hamming Encoding. BER = " + num2str(calcBer(migImageBits, recvHamming)));
showImage(recvInterlaced, "Transmitted Image with Interlacing. BER = " + num2str(calcBer(migImageBits, recvInterlaced)));
showImage(recvHammingInterlaced, "Transmitted Image with Hamming Encoding and Interlacing. BER = " + num2str(calcBer(migImageBits, recvHammingInterlaced)));

snrs = [-15, -15, -5, 0, 5, 10, 15];
berNone = zeros(1, 7);
berHamInt = zeros(1, 7);
for i=1:7
    berNone(i) = calcBer(migImageBits, simulateModel(migImageBits, snrs(i), snrShot, pShot, shotLength, false, false));
    berHamInt(i) = calcBer(migImageBits, simulateModel(migImageBits, snrs(i), snrShot, pShot, shotLength, true, true));
end

semilogy(snrs, berNone, snrs, berHamInt);
xlabel("SNR");
ylabel("BER");
title("BER vs SNR");
legend("No encoding or interlacing", "Hamming encoded and interlaced");

% Displays image
function showImage(imgBits, ttl)
    imgData = num2str(imgBits);
    imgData = imgData(~isspace(imgData));
    imgData = reshape(imgData, 8, []).';
    imgData = bin2dec(imgData);
    imgData = reshape(imgData, 64, []);

    figure; 
    imshow(imgData);
    colormap(jet); 
    imagesc(imgData);
    shading flat;
    colorbar;
    title(ttl);
end

% Main simulation routine
function recbits = simulateModel(bits, snr, snrShot, pShot, shotLength, useHamming, useInterlacing)
    sndbits = bits;

    if useHamming
        sndbits = encodeHamming74(sndbits);
    end

    if useInterlacing
        sndbits = interlaceBits(sndbits, 8, 8);
    end

    recbits = simulateBitTransmissionWithShotNoise(sndbits, snr, snrShot, pShot, shotLength);

    if useInterlacing
        recbits = deinterlaceBits(recbits, 8, 8);
    end

    if useHamming
        recbits = decodeHamming74(recbits);
    end
end

% Calculates ber given vector of errors
function ber = calcBer(orig, err)
    ber = sum(xor(orig, err))/length(orig);
end

% Encodes vector of bits using Hamming 7,4 code
function encoded = encodeHamming74(bits)
    encode74Matrix = ...
    [1 0 0 1 0 1 1;
     0 1 0 1 0 1 0;
     0 0 1 1 0 0 1;
     0 0 0 0 1 1 1].';
    
    numBlocks = length(bits)/4;
    blocks = mod(encode74Matrix*reshape(bits, 4, numBlocks), 2);
    encoded = reshape(blocks, 1, numBlocks*7);
end

% Decodes Hamming 7,4 encoded bits
function decoded = decodeHamming74(encodedBits)
    syndromeMatrix = ...
    [1 1 1 1 0 0 0
     1 1 0 0 1 1 0
     1 0 1 0 1 0 1];

    numBlocks = length(encodedBits)/7;
    blocks = reshape(encodedBits, 7, numBlocks);
    errorPositions = [4 2 1]*mod(syndromeMatrix*blocks, 2);

    for i=find(errorPositions)
        pos = 8 - errorPositions(i);
        blocks(pos, i) = 1 - blocks(pos, i);
    end

    decoded = reshape(blocks([1 2 3 5], :), 1, numBlocks*4);
end

% Interlaces bits
function interlaced = interlaceBits(bits, blockLength, numBlocks)
    bitsPerChunk = blockLength*numBlocks;
    interlaced = zeros(1, length(bits));
    for i=1:bitsPerChunk:length(bits)
        range = i:(i + bitsPerChunk - 1);
        interlaced(range) = reshape(reshape(bits(range), blockLength, numBlocks).', 1, []);
    end
end

% Deinterlaces bits
function deinterlaced = deinterlaceBits(bits, blockLength, numBlocks)
    bitsPerChunk = blockLength*numBlocks;
    deinterlaced = zeros(1, length(bits));
    for i=1:bitsPerChunk:length(bits)
        range = i:(i + bitsPerChunk - 1);
        deinterlaced(range) = reshape(reshape(bits(range), numBlocks, blockLength).', 1, []);
    end
end

% Simulates QPSK transmission of bits with noise and shot noise
function recbits = simulateBitTransmissionWithShotNoise(bits, snr, snrShot, pShot, shotLength)
    % Initialize values
    samplesPerBit = 4;
    snrOffset = 10*log10(2/samplesPerBit);
    numBits = length(bits);
    numSymbols = ceil(numBits/2);
    numSamples = numSymbols * samplesPerBit;
    pairs = reshape(bits, 2, numSymbols);
    recbits = zeros(2, numSymbols);

    aoft = (pairs(1,:)*2-1) + 1i*(pairs(2,:)*2-1); % Change to gray coded complex symbol
    aoft = repmat(aoft,samplesPerBit,1); % Make into 4 samples per symbol
    aoft = reshape(aoft,1,numSamples); % Reshape into a row array
    aoft = complex(aoft);

    % Add noise with geometrically distribued shots of equal length
    i = 1;
    while i <= numSamples
        % Add regular noise
        samplesUntilShot = randGeo(pShot);
        sampleRange = i:min(i + samplesUntilShot, numSamples);
        aoft(sampleRange) = awgn(aoft(sampleRange), snr + snrOffset, 'measured');
        i = i + samplesUntilShot + 1;

        if i > numSamples
            break
        end

        % Add shot noise
        %shotLength = randGeo(pShotEnd);
        sampleRange = i:min(i + shotLength * 4, numSamples);
        aoft(sampleRange) = awgn(aoft(sampleRange), snrShot + snrOffset, 'measured');
        i = i + shotLength * 4 + 1;
    end
    
    % Reciever
    for i=1:numSymbols
        sumsym = 0;
        for i1 = 1:samplesPerBit
            sumsym = sumsym+aoft((i-1)*samplesPerBit+i1);
        end
    
        if real(sumsym) > 0
            recbits(1,i) = 1;
        else
            recbits(1,i) = 0;
        end

        if imag(sumsym) > 0
            recbits(2,i) = 1;
        else
            recbits(2,i) = 0;
        end
    end

    recbits = reshape(recbits, 1, numBits);
end

% Returns a random number from a geometric distribution given parameter p
function n = randGeo(p)
    if p == 0
        n = intmax("int32");
        return
    end
    n = floor(log(1 - rand())/log(1 - p));
end
