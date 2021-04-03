clear; clc; close all;

%%
% Numerical test for the computation of eigenfaces. The faces need to be
% provided in the file faces.txt (stored columnswise, and seperated by ',')

% define the problem
[problem, param] = setUpClass.Eigenfaces();

% call the method
[Va, info] = KVCBO(problem, param);

% postprocessing
postProcessingEigenface(Va, info, problem, 45, 64)


%%
function postProcessingEigenface(Va, info, problem, pixelx, pixely)

t = datetime;
folder = sprintf('EF_%i_%i_%i_%i_%i_%i', day(t), month(t), year(t), hour(t), minute(t), round(second(t),0));
mkdir(folder);

%% show Eigenfaces
figure
P1 = reshape(Va, pixely, pixelx);
minP1 = min(P1(:));
P1 = double(P1 - minP1) ./ double( max(P1(:)) - minP1 );
P1 = im2double(P1);
subplot(1,3,1)
imshow(P1);
imwrite(P1, 'Eigenface1.png');
movefile('Eigenface1.png', folder);

% figure
P2 = reshape(-Va, pixely, pixelx);
minP2 = min(P2(:));
P2 = double(P2 - minP2) ./ double( max(P2(:)) - minP2 );
P2 = im2double(P2);
subplot(1,3,2)
imshow(P2);
imwrite(P2, 'Eigenface2.png');
movefile('Eigenface2.png', folder);

%% SVD Eigenface
Vmin = problem.Vmin;
PSVD = reshape(Vmin, pixely, pixelx);
minPSVD = min(PSVD(:));
PSVD = double(PSVD - minPSVD) ./ double(max(PSVD(:)) - minPSVD);
PSVD = im2double(PSVD);
subplot(1,3,3)
imshow(PSVD);
imwrite(PSVD, 'SVD.png');
movefile('SVD.png', folder);

%% PS2N
cd(folder)
a0 = imread('SVD.png');
b01 = imread('Eigenface1.png');
b02 = imread('Eigenface2.png');

fprintf('\n *****************************************\n\n')
fprintf(' Peak Signal-to-Noise: %.2f; %.2f', psnr(a0,b01), psnr(a0,b02))
cd('..')

%% Figures
figure;
subplot(1,3,1)
semilogy([info.iter], [info.cost], '.-')
xlabel('Iterations')
ylabel('Cost')

subplot(1,3,2)
semilogy([info.iter], [info.error], '.-');
xlabel('Iterations')
ylabel('Error')

subplot(1,3,3)
semilogy([info.iter], [info.gradnorm], '.-');
xlabel('Iterations')
ylabel('Gradnorm')

print('Figure','-depsc')
movefile('Figure.eps', folder);

%% 

fprintf('\n *****************************************\n\n')
fprintf(' Created folder \n\n %s/%s\n\n', pwd, folder)


end
