clear; clc; close all;

%  This script will create a new folder RESIZED with all the
%  resized black and white pictures and export them in a txt-file

%%

folder = 'DATA_Test';
extension = '.jpg';

%% counting & renaming

n = countFiles(folder, extension); 
% renameAllFiles(folder, extension);                    % rename the images as image_%i.extension

%% resizing & grayscale

height = 64; 
length = 45;

cd(folder)
mkdir('RESIZED')

files = dir(strcat('*', extension));
j = 1;

for file = files'
    
    x = imread(file.name);
    p = imresize(x, [height length]); 
    p = rgb2gray(p);
    p = im2double(p);
    
    % average pixels to reduce size
    % p = AveragePixels(p);
    
    tmp = sprintf('resized_%i.jpg',j); 
    j = j+1;
    
    imwrite(p, tmp);
    movefile(tmp,'RESIZED');
end

%% Export images as .txt

cd('RESIZED')
folder = pwd;

X = loadFile(height, length, n, folder);
writematrix(X, 'Filename.txt', 'Delimiter', ',');
cd ..

%% Check

P1 = reshape(X(:,1), height, length);
K = mat2gray(P1);
figure
imshow(K)

return

%% Attach two .txt files

cd(folder)
A = loadFromText('FileA.txt');
B = loadFromText('FileB.txt');

C = [A B];

writematrix(C, 'Filename.txt', 'Delimiter', ',');


function X = loadFile(height, length, n, folder)

    cd(folder)
    Pictures = zeros(height*length,n);

    for i = 1:n
        container = imread(sprintf('resized_%i.jpg',i));
        Pictures(:,i) = reshape(container, height*length, 1);
    end
    cd ..
  
    X = Pictures;
        
end

function n = countFiles(folder, extension)
  foldername = strcat(pwd,'/',folder);
  ext = strcat('/*', extension);
  a = dir([foldername ext]);
  n = size(a,1);
end


function [h,l] = getSize(filename, folder)
  cd(folder)
  tmp = imread(filename);
  h = size(tmp,1);
  l = size(tmp,2);
  if folder ~= pwd 
      cd ..
  end
end

function X = loadFromText(filename)
  headerLines = 0;
  X = importdata(filename, ',', headerLines);
end

function pr = AveragePixels(p)

  h = size(p,1); 
  l = size(p,2);
  pr = double(zeros(ceil(h/2), ceil(l/2)));
  
  for i=1:(h/2)
     for j=1:(l/2)   
         pr(i,j) = 1/4*(p(2*(i-1)+1,2*(j-1)+1) + p(2*(i-1)+2,2*(j-1)+1) + p(2*(i-1)+1,2*(j-1)+2) + p(2*(i-1)+2,2*(j-1)+2));
     end
  end
end

function renameAllFiles(folder, extension)
    cd(folder)
    files = dir(strcat('*', extension));
    j=1;
    for file = files'
        movefile(file.name, strcat(sprintf('image_%i',j),'.jpg')); j = j+1;
    end
    cd ..
end

