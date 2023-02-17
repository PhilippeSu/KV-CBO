%%

clear; clc; close all;

ref = imread('Eigenface_SVD_0outliers.png');
A = imread('Eigenface_SVD_15outliers_1.png');
B = imread('Eigenface_KVCBO_15outliers_1.png');
C = imread('Eigenface_FMS_15outliers_1.png');
D = imread('Eigenface_GD_15outliers_1.png');


fprintf('SVD\n*********\n')
[peaksnr, snr] = psnr(ref,A)
ssimA = ssim(ref,A)

fprintf('KVCBO\n*********\n')
[peaksnr, snr] = psnr(ref,B)
ssimA = ssim(ref,B)

fprintf('FMS\n*********\n')
[peaksnr, snr] = psnr(ref,C)
ssimA = ssim(ref,C)

fprintf('GD\n*********\n')
[peaksnr, snr] = psnr(ref,D)
ssimA = ssim(ref,D)
