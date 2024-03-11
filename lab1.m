%% 2.1 Contrast Stretching
close all hidden;

%% a.
Pc = imread('mrt-train.jpg');
whos Pc
P = rgb2gray(Pc);

%% b.
figure;
imshow(P);

%% c.
disp('image P');
minP = double(min(P(:))); % min intensity
disp(['min intensity: ', num2str(minP)]);
maxP = double(max(P(:))); % max intensity
disp(['max intensity: ', num2str(maxP)]);

%% d.
P2 = imsubtract(P, minP);
P2 = P2 * (255 / (maxP - minP));
disp('image P2');
minP2 = min(P2(:)); % min intensity
disp(['min intensity: ', num2str(minP2)]);
maxP2 = max(P2(:)); % max intensity
disp(['max intensity: ', num2str(maxP2)]);

%% e.
figure;
imshow(P2);

%% 2.2 Histogram Equalization

%% a.
figure;
imhist(P, 10);
figure;
imhist(P, 256);

%% b.
P3 = histeq(P, 255);
figure;
imhist(P3, 10);
figure;
imhist(P3, 256);

%% c.
P4 = histeq(P3, 255);
figure;
imhist(P4, 10);
figure;
imhist(P4, 256);

%% 2.3 Linear Spatial Filtering

%% a.
h1 = fspecial('gaussian', 5, 1);
sum(h1, 'all')
h2 = fspecial('gaussian', 5, 2);
sum(h2, 'all')
figure;
mesh(h1);
figure;
mesh(h2);

%% b.
ntugn = imread('ntu-gn.jpg');
figure;
imshow(ntugn);

%% c.
ntugn1 = uint8(conv2(ntugn, h1));
figure;
imshow(ntugn1);
ntugn2 = uint8(conv2(ntugn, h2));
figure;
imshow(ntugn2);

%% d.
ntusp = imread('ntu-sp.jpg');
figure;
imshow(ntusp);

%% e.
ntusp1 = uint8(conv2(ntusp, h1));
figure;
imshow(ntusp1);
ntusp2 = uint8(conv2(ntusp, h2));
figure;
imshow(ntusp2);

%% 2.4 Median Filtering
ntugn3 = medfilt2(ntugn, [3 3]);
figure;
imshow(ntugn3);
ntugn4 = medfilt2(ntugn, [5 5]);
figure;
imshow(ntugn4);
ntusp3 = medfilt2(ntusp, [3 3]);
figure;
imshow(ntusp3);
ntusp4 = medfilt2(ntusp, [5 5]);
figure;
imshow(ntusp4);

%% 2.5 Suppressing Noise Interference Patterns

%% a.
pckint = imread('pck-int.jpg');
figure;
imshow(pckint);

%% b.
F = fft2(pckint);
S = abs(F).^2;
figure;
imagesc(fftshift(S.^0.1));
colormap('default');

%% c.
figure;
imagesc(S.^0.1);
[x, y] = ginput(2);

%% d.
F(y(1)-2:y(1)+2, x(1)-2:x(1)+2) = 0;
F(y(2)-2:y(2)+2, x(2)-2:x(2)+2) = 0;
S = abs(F).^2;
figure;
imagesc(fftshift(S.^0.1));
colormap('default');

%% e.
FI = uint8(ifft2(F));
figure;
imshow(FI);

%%
P = imread('pck-int.jpg');
F = fft2(P);
F(y(1)-6:y(1)+6, x(1)-6:x(1)+6) = 0;
F(y(2)-6:y(2)+6, x(2)-6:x(2)+6) = 0;
S = abs(F).^2;
figure;
imagesc(fftshift(S.^0.1));
colormap('default');
FI = uint8(ifft2(F));
figure;
imshow(FI);

%% f.
%close all hidden;
primatecaged = rgb2gray(imread('primate-caged.jpg'));
figure;
imshow(primatecaged);
F2 = fft2(primatecaged);
S2 = abs(F2).^2;
figure;
imagesc(fftshift(S2.^0.1));

%%
F2 = fft2(primatecaged);
S2 = abs(F2).^2;
figure;
imagesc(S2.^0.1);
[X, Y] = ginput(8);

%%
F2(Y(1)-2:Y(1)+2, X(1)-2:X(1)+2) = 0;
F2(Y(2)-2:Y(2)+2, X(2)-2:X(2)+2) = 0;
F2(Y(3)-2:Y(3)+2, X(3)-2:X(3)+2) = 0;
F2(Y(4)-2:Y(4)+2, X(4)-2:X(4)+2) = 0;
F2(Y(5)-2:Y(5)+2, X(5)-2:X(5)+2) = 0;
F2(Y(6)-2:Y(6)+2, X(6)-2:X(6)+2) = 0;
F2(Y(7)-2:Y(7)+2, X(7)-2:X(7)+2) = 0;
F2(Y(8)-2:Y(8)+2, X(8)-2:X(8)+2) = 0;
S2 = abs(F2).^2;
figure;
imagesc(fftshift(S2.^0.1));
colormap('default');
FI2 = uint8(ifft2(F2));
figure;
imshow(FI2);

%% 2.6 Undoing Perspective Distortion of Planar Surface 

%% a.
P = imread('book.jpg');
figure;
imshow(P);

%% b.
[X,Y] = ginput(4);
imageX = [0 210 210 0];
imageY = [0 0 297 297];

%% c.
A = [
  [X(1),Y(1),1,0,0,0, -imageX(1)*X(1),-imageX(1)*Y(1)];
  [0,0,0,X(1),Y(1),1, -imageY(1)*X(1),-imageY(1)*Y(1)];
  [X(2),Y(2),1,0,0,0, -imageX(2)*X(2),-imageX(2)*Y(2)];
  [0,0,0,X(2),Y(2),1, -imageY(2)*X(2),-imageY(2)*Y(2)];
  [X(3),Y(3),1,0,0,0, -imageX(3)*X(3),-imageX(3)*Y(3)];
  [0,0,0,X(3),Y(3),1, -imageY(3)*X(3),-imageY(3)*Y(3)];
  [X(4),Y(4),1,0,0,0, -imageX(4)*X(4),-imageX(4)*Y(4)];
  [0,0,0,X(4),Y(4),1, -imageY(4)*X(4),-imageY(4)*Y(4)];
];
v = [imageX(1); imageY(1); imageX(2); imageY(2); imageX(3); imageY(3); imageX(4); imageY(4)];
u = A \ v;
U = reshape([u;1], 3, 3)';
w = U*[X'; Y'; ones(1,4)];
w = w ./ (ones(3,1) * w(3,:));

%% d.
T = maketform('projective', U');
P2 = imtransform(P, T, 'XData', [0 210], 'YData', [0 297]); 

%% e.
figure;
imshow(P2);