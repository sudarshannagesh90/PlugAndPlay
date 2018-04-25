clc
close all
clear all
%%
im  = im2double(imread('saturn.png'));
im1 = zeros(256,256,3);

for ind = 1:3
    im1(:,:,ind) = imresize(im(:,:,ind),[256 256]);
end
im1 = uint8(im1*255);
imshow(im1,[])
%%
h5create('imagenet_64.h5','/data',[3 256 256 1])
data(1,:,:,1) = im1(:,:,1)';
data(2,:,:,1) = im1(:,:,2)';
data(3,:,:,1) = im1(:,:,3)';
h5write('imagenet_64.h5', '/data', data)
%%