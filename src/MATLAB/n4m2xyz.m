% Using a standard UIUC .dat file finds the modified coefficients described
% at http://www.pdas.com/naca456thick4m.html, and plots an xyz file using
% 'real' co-ordinates at a desired resolution of points.
%% initialization
clc;
clear all;
close all;
%% vars
% get data from .dat file
data_struct     =   importdata(uigetfile('*.dat'),' ',1);
data            =   data_struct.data;
% define desired airfoil resolution
n_pts           =   1e3;
% find maximum thickness array index
max_t_ind       =   find(data(:,1) <= 0.4,1);
% get 4 points before max thickness
x_b             =   data(max_t_ind:max_t_ind+3,1);
y_b             =   data(max_t_ind:max_t_ind+3,2);
% get 4 points after max thickness
x_a             =   data([1 max_t_ind-2:max_t_ind],1);
y_a             =   data([1 max_t_ind-2:max_t_ind],2);
%% analyze
% solve for the a0..a3 coefs for the before equation
A_b     =   [sqrt(x_b(1)) x_b(1) x_b(1)^2 x_b(1)^3;
            sqrt(x_b(2)) x_b(2) x_b(2)^2 x_b(2)^3;
            sqrt(x_b(3)) x_b(3) x_b(3)^2 x_b(3)^3;
            sqrt(x_b(4)) x_b(4) x_b(4)^2 x_b(4)^3;];
as      =   A_b\y_b;
% solve for the d0..d3 coefs for the after equation
A_a     =   [1 1-x_a(1) (1-x_a(1))^2 (1-x_a(1))^3;
            1 1-x_a(2) (1-x_a(2))^2 (1-x_a(2))^3;
            1 1-x_a(3) (1-x_a(3))^2 (1-x_a(3))^3;
            1 1-x_a(4) (1-x_a(4))^2 (1-x_a(4))^3;];
ds      =   A_a\y_a;
% govern num of points for each part, using 1:2
n_bs    = floor(n_pts/3);
n_as    = n_pts-n_bs;
% using trigonometric distribution, get x_bs and x_as
tht_c   = acos(0.6);
x_bs    = 1-cos(tht_c/(n_bs-1)*(0:(n_bs-1)));
y_bs    = as(1)*sqrt(x_bs)+as(2)*x_bs+as(3)*x_bs.^2+as(4)*x_bs.^3;
tht_s   = asin(0.4);
x_as    = sin(tht_s+(pi/2-tht_s)*(0:(n_as-1))/(n_as-1));
y_as    = ds(1)+ds(2)*(1-x_as)+ds(3)*(1-x_as).^2+ds(4)*(1-x_as).^3;
%% save to file
% convert to 'real' coordinates
z_s     = -[y_bs y_as 0]';
x_s     = -[x_bs x_as 1]';
y_s     = zeros(size(x_s));
output  = [x_s, y_s, z_s];
output2 = [x_s-2.86, y_s+4, z_s];
save([data_struct.textdata{1},data_struct.textdata{2},'R.txt'],...
    'output','-ASCII');
save([data_struct.textdata{1},data_struct.textdata{2},'T.txt'],...
    'output2','-ASCII');