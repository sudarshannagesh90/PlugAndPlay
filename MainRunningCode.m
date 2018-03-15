clc
close all
clear all
%% Input parameters 
r        = im2double(imread('cameraman.tif'));
sigma_w  = 0.1;
A        = 'fft';
seedNum  = 1;
%% Generate w,g, and y form input-parameters
[M,N]    = size(r);                                                        % Size of the input-image 
g        = (sqrt(r)/2).*randn(M,N)+1j*(sqrt(r)/2).*randn(M,N);             % g ~ CN(0,D(r))




