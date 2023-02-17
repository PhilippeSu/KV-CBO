

ref = imread('pout.tif');
imshow(ref)
A = imnoise(ref,'gaussian', 0.2);
figure
imshow(A)

[peaksnr, snr] = psnr(A, ref);
  
fprintf('\n The Peak-SNR value is %0.4f', peaksnr);

fprintf('\n The SNR value is %0.4f \n', snr);
