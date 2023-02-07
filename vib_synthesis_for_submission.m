%%
import cls_machine.*
mac = cls_machine();
% load mode data
load("mode_struct.mat");
% load force data
load("force_struct.mat");
obj_f = f_no_load;  % no-load force
% obj_f = f_on_load;  % on-load force
%%
freq_order = 2; % Define the freq_order. Single order, e.g. 1, 2, .... Or, multiple orders input as array, e.g. 1:100.
speed = 1500; % Define the rotating speed in r/min.
i_max = 200; % Define the considered maximum spatial orders.
order_r = zeros(2*i_max, 1);
acc_rst = zeros(length(freq_order),1);
for i = 1:i_max
    order_r(2*i-1) = 12 * (i - 1) - 2;
    order_r(2*i) = 12 * (i - 1) + 2;
end
order_r = cat(1, [1, 3]', order_r(2:end));
sig_v = 0:(pi / 60):(2 * pi); % Define the phase shift.
% choose the mesh
mesh = 1; % If mesh1 is chosen, mesh=1
mev1 = mev1_mesh1; % If mesh1 is chosen, mev1 = mev_mesh1
mev2 = mev2_mesh1; % If mesh1 is chosen, mev2 = mev_mesh1
% Response values and modal frequencies under five meshes
phi_l_mode1 = [5.51, -5.51, 6.09, -6.19, -6.19];
freq_mode1 = [1299, 1306, 1310, 1298, 1297];
phi_l_mode2 = [-14.62, -14.57, -14.36, -14.31, 14.32];
freq_mode2 = [1328, 1339, 1342, 1336, 1336];
for i = 1:length(freq_order)
    v = freq_order(i);
    % calculate the response of the first 2-order mode
    a_v_u_r_a1_n1 = ...
        extn_acc2(freq_mode1(mesh), mev1.mev_rd, obj_f.fd_rd, mac, order_r, order_r, v, speed, sig_v, phi_l_mode1(mesh)) + ...
        extn_acc2(freq_mode1(mesh), mev1.mev_rd, obj_f.fd_rd, mac, order_r, -order_r, v, speed, sig_v, phi_l_mode1(mesh));
    a_v_u_r_a2_n1 = ...
        extn_acc2(freq_mode1(mesh), mev1.mev_tg, obj_f.fd_tg, mac, order_r, order_r, v, speed, sig_v, phi_l_mode1(mesh)) + ...
        extn_acc2(freq_mode1(mesh), mev1.mev_tg, obj_f.fd_tg, mac, order_r, -order_r, v, speed, sig_v, phi_l_mode1(mesh));
    % calculate the response of the fsecond 2-order mode
    a_v_u_r_a1_n2 = ...
        extn_acc2(freq_mode2(mesh), mev2.mev_rd, obj_f.fd_rd, mac, order_r, order_r, v, speed, sig_v, phi_l_mode2(mesh)) + ...
        extn_acc2(freq_mode2(mesh), mev2.mev_rd, obj_f.fd_rd, mac, order_r, -order_r, v, speed, sig_v, phi_l_mode2(mesh));
    a_v_u_r_a2_n2 = ...
        extn_acc2(freq_mode2(mesh), mev2.mev_tg, obj_f.fd_tg, mac, order_r, order_r, v, speed, sig_v, phi_l_mode2(mesh)) + ...
        extn_acc2(freq_mode2(mesh), mev2.mev_tg, obj_f.fd_tg, mac, order_r, -order_r, v, speed, sig_v, phi_l_mode2(mesh));
    acc_sum = sum(a_v_u_r_a1_n1, 1) + sum(a_v_u_r_a2_n1, 1) + sum(a_v_u_r_a1_n2, 1) + sum(a_v_u_r_a2_n2, 1);
    acc_rst(i) = max(real(acc_sum));
end
[acc_max_amp, pos] = max(real(acc_sum));
acc_max_ang = angle(acc_sum(pos));
acc_max_vec = a_v_u_r_a1_n1(:, pos);
acc_max_vec_real = real(acc_max_vec);
plot(freq_order,acc_rst,Marker="*");

%%
function [a_v_u_r_a1_n] = extn_acc2(fr, mev_rd, fd_rd, mac, r, u, v, speed, sig_v, phi_l)
% *********************************************************
% fr: modal frequncies in Hz;
% mev_rd: sampling sequence of mode field, with 36000 spatial points in this case study;
% fd_rd: sampling sequence of force field, with 36000 spatial points and 120 tempral points in this case study;
% *********************************************************
% fft of mode field
[fft_mev1.mev_rd_amp, fft_mev1.mev_rd_ang] = fft_plot(36000, 36000, mev_rd); % conduct 1-d fft
% convert basis modal data
wr1 = fr * 2 * pi;
gr1 = 0.03; % define modal damp
kr1 = wr1^2; % modal stiffness
cr1 = gr1 / wr1 * kr1;
mr1 = 1;
% define basic parameters of machines
w_e = speed / 60 * (mac.num_poles * 0.5) * 2 * pi; % fundamental frequencies
r_s = mac.r_outgap; % inner radius of stator
l_stk = mac.l_stk; % stack length
% fft of force field
l_lay = 40; % axial layers of elements
X1 = fft2(fd_rd*l_stk/l_lay); % conduct 2-d fft
X1 = fftshift(X1) / 36000 / 120;
px = 36000 / 2 + 1;
pt = mac.num_fea_step / 2 + 1;
% get Fourier coefficients of forces
f_v_1u_a1 = abs(X1(px+u, pt+v));
f_v_0u_a1 = abs(X1(px-u, pt+v));
gam_v_1u_a1 = angle(X1(px+u, pt+v));
gam_v_0u_a1 = angle(X1(px-u, pt+v));
% get Fourier coefficients of modes
phi_r_a1 = fft_mev1.mev_rd_amp(r+1) * 0.5;
chi_r_a1 = fft_mev1.mev_rd_ang(r+1);
% calculate acceleration
w = v * w_e;
h_n = 1 / (kr1 - w^2 * mr1 + 1j * w * cr1);
a_v_u_r_a1_n = -w^2 * 2 * pi * r_s * h_n * phi_l * phi_r_a1 ...
    .* (sinc(r+u) .* exp(1j*chi_r_a1) + sinc(r-u) .* exp(-1j*chi_r_a1)) ...
    .* (f_v_1u_a1 .* exp(1j*gam_v_1u_a1+1j*sig_v) + f_v_0u_a1 .* exp(-1j*gam_v_0u_a1-1j*sig_v));
a_v_u_r_a1_n = a_v_u_r_a1_n * 1e-3; % unit transform, from mm to m.
end

function [amp, ang] = fft_plot(varargin1, varargin2, varargin3, varargin4) %function for 1-d fft
if nargin == 4
    N = varargin1;
    L = varargin2;
    sg = varargin3;
    flag = varargin4;
else
    N = varargin1;
    L = varargin2;
    sg = varargin3;
    flag = 0;
end
sp = fft(sg, L);
sp_amp = abs(sp);
sp_ang = angle(sp);
f_x = 0:N / L:N - N / L;
f_x = f_x(1:L/2+1);
a_y = sp_amp(1:L/2+1);
a_y = a_y / L;
a_y(2:end) = 2 * a_y(2:end);
b_y = sp_ang(1:L/2+1);
if flag == 1
    figure;
    bar(f_x, a_y);
elseif flag == 2
    figure;
    bar(f_x, a_y);
    figure;
    bar(f_x, b_y);
end
amp = a_y;
ang = b_y;
end