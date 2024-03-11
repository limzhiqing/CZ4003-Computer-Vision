%% 3.1. Edge Detection
close all hidden;

%% a.
mcr_rgb = imread("macritchie.jpg");
figure;
imshow(mcr_rgb);
mcr_g = rgb2gray(mcr_rgb);
figure;
imshow(mcr_g);

%% b. 
sobel_h = [-1 -2 -1; 
            0  0  0; 
            1  2  1];
sobel_v = [-1  0  1; 
           -2  0  2; 
           -1  0  1;];

%%
mcr_sh = conv2(mcr_g, sobel_h); 
figure;
imshow(uint8(mcr_sh));

%%
mcr_sv = conv2(mcr_g, sobel_v);
figure;
imshow(uint8(mcr_sv));

%% c.
mcr_shvsc = mcr_sh.^2 + mcr_sv.^2;
figure;
imshow(uint8(mcr_shvsc));

%% d.
E = mcr_shvsc;
Et1 = E>10000;
figure;
imshow(Et1);
Et2 = E>20000;
figure;
imshow(Et2);
Et3 = E>40000;
figure;
imshow(Et3);
Et4 = E>80000;
figure;
imshow(Et4);

%% e.i.
I = mcr_g;
tl = 0.04;
th = 0.1;
sigma = 1.0;
E1 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(E1);
sigma = 2.3;
E2 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(E2);
sigma = 3.7;
E3 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(E3);
sigma = 5.0;
E4 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(E4);

%% e.ii.
tl = 0.01;
th = 0.1;
sigma = 1.0;
Etl1 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(Etl1);
tl = 0.04;
Etl2 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(Etl2);
tl = 0.07;
Etl3 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(Etl3);
tl = 0.09;
Etl4 = edge(I, 'canny', [tl th], sigma);
figure;
imshow(Etl4);

%% 3.2. Line Finding using Hough Transform
close all hidden;

%% a.
I = mcr_g;
tl = 0.04;
th = 0.1;
sigma = 1.0;
E = edge(I, 'canny', [tl th], sigma);

%% b.
[H, xp] = radon(E);
figure;
imshow(uint8(H));

%% c.
[radius, theta] = find(H ==max(H(:)));
disp(radius);
disp(theta);

%% d.
radius = xp(radius);

[A, B] = pol2cart(theta*pi/180, radius);
B = -B;

whos mcr_g;

%%
C = A*(A+179) + B*(B+145);

disp(A);
disp(B);
disp(C);

%% e.
xl = 0;
xr = 357;
yl = (C - A * xl) / B;
yr = (C - A * xr) / B;
disp(yl);
disp(yr);

%% f.
figure;
imshow(I);
line([xl xr], [yl yr]);

%%
close all hidden;
I = mcr_g;
tl = 0.04;
th = 0.1;
sigma = 9.0;
E = edge(I, 'canny', [tl th], sigma);
figure;
imshow(E);
[H, xp] = radon(E);
[radius, theta] = find(H ==max(H(:)));
radius = xp(radius);
theta = theta + 0.85;
[A, B] = pol2cart(theta*pi/180, radius);
B = -B;
figure;
imshow(H);
C = A*(A+179) + B*(B+145);
xl = 0;
xr = 357;
yl = (C - A * xl) / B;
yr = (C - A * xr) / B;
figure;
imshow(I);
line([xl xr], [yl yr]);

%% 3.3. 3D Stereo
close all hidden;

%% a.
% see below

%% b.
cl_rgb = imread("corridorl.jpg");
figure;
imshow(cl_rgb);
cl_g = im2double(rgb2gray(cl_rgb));
figure;
imshow(cl_g);

cr_rgb = imread("corridorr.jpg");
figure;
imshow(cr_rgb);
cr_g = im2double(rgb2gray(cr_rgb));
figure;
imshow(cr_g);

%% c.
close all hidden;
D = disparity_map(cl_g, cr_g, 11, 11);
figure;
imshow(-D,[-15 15]);
cd = imread('corridor_disp.jpg');
figure;
imshow(cd);

%% d.
close all hidden;
tcl= imread("triclopsi2l.jpg");
tclg = im2double(rgb2gray(tcl));
figure;
imshow(tclg);
tcr = imread("triclopsi2r.jpg");
tcrg = im2double(rgb2gray(tcr));
figure;
imshow(tcrg);
D2 = disparity_map(tclg, tcrg, 11, 11);
figure;
imshow(-D2,[-15 15]);
tcd =  imread("triclopsid.jpg");
figure;
imshow(tcd);

%% disparity map function
function disparity_map = disparity_map(left, right, dimx, dimy)
    [x y] = size(left);
    tempCorr = zeros(x, y);
    dispLimit = 15;
    templateRegion = ones(dimx, dimy);
    pxCorr = zeros(x, y, dispLimit);
    for d = 0 : dispLimit - 1
       rmCol = x * d + 1 : x * y;
       tempCorr(:) = Inf;
       corr_dlr = left(rmCol) - right(rmCol - x * d);
       tempCorr(rmCol) = corr_dlr.^2;
       pxCorr(:, :, d + 1) = conv2(tempCorr, templateRegion, 'same');
    end
    [v, dm] = min(pxCorr, [], 3);
    disparity_map = dm * -1;
end
