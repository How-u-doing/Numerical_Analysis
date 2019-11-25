% Image compression
% According to Data4Bio, YouTube lectures by University of Washington

% Load full image
disp('Loading full image...')
A=imread('Keira_Knightly.jpg');

figure(3) 
imshow(A)

% make black and white image (gray)
Abw=rgb2gray(A);
[nx, ny]=size(Abw);
figure(1), subplot(2,2,1), imshow(Abw)
title('Original Image','FontSize',16)

% compute FFT of our image using fft2
disp('Doing FFT analysis for sparsity check...')
Ahat=fft2(Abw);
F=log(abs(fftshift(Ahat)+1));
F=mat2gray(F);
figure(4)
imshow(F,[]); % display result

% zero out small coefficients and inverse transform
count_pic=2;
for threshold=.1*[0.001 0.005 0.01] * max(abs(Ahat(:)))
    idx=abs(Ahat)>threshold;
    AhatFilt=Ahat.*idx;
    percent=sum(idx(:))/(nx*ny)*100;
    Afilt=ifft2(AhatFilt);
    figure(1), subplot(2,2,count_pic)
    imshow(uint8(Afilt))
    count_pic=count_pic+1;
    drawnow
    title([num2str(percent) '% of FFT basis'], 'Fontsize',14)
end

% plot image as mountain surface
figure
Anew=imresize(Abw,0.1);
% overlook or look upward, you're capaple of seeing gorgeous Keira Knightly
surf(double(Anew))



