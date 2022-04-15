clc; clearvars; close all;

targetSig = [1, 2, -1, 3, 5, -6];
expandedSig = [zeros(1, 4), targetSig, zeros(1, 3)];

figure;
subplot(2, 1, 1);
plot(targetSig); grid on;
xlim([1, 20]);
subplot(2, 1, 2);
plot(expandedSig); grid on;
xlim([1, 20]);

disp(expandedSig);

[corr, delays] = xcorr(targetSig, expandedSig);

figure; 
subplot(2, 1, 1);
plot(abs(corr)); grid on;
subplot(2, 1, 2);
plot(abs(delays)); grid on;

[~, Xcoordinate] = max(corr);

cutoff = expandedSig(abs(delays(Xcoordinate)) + 1: ...
    abs(delays(Xcoordinate)) + length(targetSig));

disp(cutoff);