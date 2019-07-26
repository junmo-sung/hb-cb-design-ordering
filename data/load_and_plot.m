close all
clear
clc

addpath iter_500

figureSize = [8,5]; %inches
color_a = [     0  0.4470  0.7410];
color_b = [0.8500  0.3250  0.0980];
color_c = [0.9290, 0.6940, 0.1250];
color_d = [0.4940, 0.1840, 0.5560];
color_e = [0.4660, 0.6740, 0.1880];
color_f = [0.3010, 0.7450, 0.9330];
color_g = [0.6350, 0.0780, 0.1840];

PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M32_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M32_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M96_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M96_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M160_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M160_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M192_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M192_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M224_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M224_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.mat', 'SP');

PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np2_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np2_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np4_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np4_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np8_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np8_Q16.mat', 'SP');
PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np10_Q16 = load('PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np10_Q16.mat', 'SP');

adaptive_Nt64_Nr16_Lt8_Lr4_G96_Np6 = load('adaptive_Nt64_Nr16_Lt8_Lr4_G96_Np6.mat');
% figure
% hold on
% plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_md)), 'x--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)
% plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_sung)), 'o--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)
% plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_lee)), 's--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)
% hold off


%% NMSE vs SNR
figure1 = figure;
set(gcf, 'defaulttextinterpreter', 'latex')
hold on
set(gca, 'FontSize', 13)
set(gcf, 'Units', 'inches')
pos = get(gcf, 'position');
pos(3:4) = figureSize;
set(gcf, 'position', pos)

plot(inf, inf, 'ko', 'MarkerSize', 10)
plot(inf, inf, 'ks', 'MarkerSize', 10)
plot(inf, inf, 'kx', 'MarkerSize', 10)
plot(inf, inf, '--', 'color', color_a, 'LineWidth', 1.1) 
plot(inf, inf, '-', 'color', color_b, 'LineWidth', 1.1) 
plot(inf, inf, ':', 'color', color_c, 'LineWidth', 1.1) 
plot(inf, inf, 'd-.', 'color', color_d, 'MarkerSize', 10, 'LineWidth', 1.1) 

plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_md)), 's--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_sung)), 'o--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_lee)), 'x--', 'color', color_a, 'LineWidth', 1.1, 'MarkerSize', 10)

plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_md)), 's-', 'color', color_b, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_sung)), 'o-', 'color', color_b, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_lee)), 'x-', 'color', color_b, 'LineWidth', 1.1, 'MarkerSize', 10)

plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_md)), 's:', 'color', color_c, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_sung)), 'o:', 'color', color_c, 'LineWidth', 1.1, 'MarkerSize', 10)
plot(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.SNR_db_array, 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_lee)), 'x:', 'color', color_c, 'LineWidth', 1.1, 'MarkerSize', 10)

plot(adaptive_Nt64_Nr16_Lt8_Lr4_G96_Np6.SNR_dBa, 10*log10(mean(adaptive_Nt64_Nr16_Lt8_Lr4_G96_Np6.NMSE)), 'd-.', 'color', color_d, 'LineWidth', 1.1, 'MarkerSize', 10)

annotation(figure1,'ellipse',...
    [0.420831223628691 0.628034482758621 0.0369746835443038 0.206896551724138],...
    'LineWidth',1);
annotation(figure1,'ellipse',...
    [0.572729957805904 0.474965517241379 0.0433037974683576 0.0658620689655169],...
    'LineWidth',1);
annotation(figure1,'ellipse',...
    [0.726738396624468 0.280758620689655 0.0433037974683576 0.103448275862068],...
    'LineWidth',1);
annotation(figure1,'textarrow',[0.508438818565401 0.462025316455696],...
    [0.841379310344828 0.813793103448276],'String',{'M_x=2 (M=64)'}, 'FontSize', 13);
annotation(figure1,'textarrow',[0.514767932489451 0.571729957805907],...
    [0.441379310344828 0.496551724137931],'String',{'M_x=4 (M=128)'}, 'FontSize', 13);
annotation(figure1,'textarrow',[0.672995780590717 0.729957805907172],...
    [0.231034482758621 0.286206896551724],'String',{'M_x=8 (M=256)'}, 'FontSize', 13);
annotation(figure1,'textarrow',[0.793402777777778 0.755208333333333],...
    [0.601777777777778 0.558333333333333],'String',{'M=864'}, 'FontSize', 13);

hold off
grid on
xlabel('SNR [dB]')
ylabel('NMSE [dB]')
% title('Nt64 / Nr16 / Lt8 / Lr8 / G64 / Np6 / Q16')
hLeg = legend('Proposed', ...
       'Random', ...
       'MTC', ...
       '$M_x=2$', ...
       '$M_x=4$', ...
       '$M_x=8$', ...
       'Adaptive CS', ...
       'Location', 'southwest');
set(hLeg, 'Interpreter', 'latex')
set(hLeg, 'FontSize', 14)

%% NMSE vs M_x (
NMSE_OMP_md = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M32_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M96_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M160_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M192_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M224_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_md))];
NMSE_OMP_lee = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M32_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M96_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M160_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M192_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M224_Np6_Q16.SP.NMSE_OMP_lee)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_lee))];
NMSE_OMP_sung = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M32_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M64_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M96_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M160_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M192_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M224_Np6_Q16.SP.NMSE_OMP_sung)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M256_Np6_Q16.SP.NMSE_OMP_sung))];

SNR_idx = 6;
figureSize = [8,3.5]; %inches

figure
set(gcf, 'defaulttextinterpreter', 'latex')
hold on
set(gca, 'FontSize', 13)
set(gcf, 'Units', 'inches')
pos = get(gcf, 'position');
pos(3:4) = figureSize;
set(gcf, 'position', pos)

plot(1:8, NMSE_OMP_sung(:, SNR_idx), 'ko-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(1:8, NMSE_OMP_md(:, SNR_idx), 'kx-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(1:8, NMSE_OMP_lee(:, SNR_idx), 'ks-', 'LineWidth', 1.1, 'MarkerSize', 10)

hold off
grid on
xlabel('$M_x$')
ylabel('NMSE [dB]')
ylim([-21, -5])
% title('Nt64 / Nr16 / Lt8 / Lr8 / G64 / Np6 / Q16')
hLeg = legend('Proposed', ...
       'Random', ...
       'MTC', ...
       'Location', 'northeast');
set(hLeg, 'Interpreter', 'latex')
set(hLeg, 'FontSize', 14)

SNR_idx = 3;
figure
set(gcf, 'defaulttextinterpreter', 'latex')
hold on
set(gca, 'FontSize', 13)
set(gcf, 'Units', 'inches')
pos = get(gcf, 'position');
pos(3:4) = figureSize;
set(gcf, 'position', pos)

plot(1:8, NMSE_OMP_sung(:, SNR_idx), 'ko-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(1:8, NMSE_OMP_md(:, SNR_idx), 'kx-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(1:8, NMSE_OMP_lee(:, SNR_idx), 'ks-', 'LineWidth', 1.1, 'MarkerSize', 10)

hold off
grid on
xlabel('$M_x$')
ylabel('NMSE [dB]')
% ylim([-21, -5])
% title('Nt64 / Nr16 / Lt8 / Lr8 / G64 / Np6 / Q16')
hLeg = legend('Proposed', ...
       'Random', ...
       'MTC', ...
       'Location', 'northeast');
set(hLeg, 'Interpreter', 'latex')
set(hLeg, 'FontSize', 14)

%% sparsity
NMSE_OMP_md = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np2_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np4_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np8_Q16.SP.NMSE_OMP_md)); ...
                10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np10_Q16.SP.NMSE_OMP_md))];
NMSE_OMP_lee = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np2_Q16.SP.NMSE_OMP_lee)); ...
                 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np4_Q16.SP.NMSE_OMP_lee)); ...
                 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_lee)); ...
                 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np8_Q16.SP.NMSE_OMP_lee)); ...
                 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np10_Q16.SP.NMSE_OMP_lee))];
NMSE_OMP_sung = [ 10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np2_Q16.SP.NMSE_OMP_sung)); ...
                  10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np4_Q16.SP.NMSE_OMP_sung)); ...
                  10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np6_Q16.SP.NMSE_OMP_sung)); ...
                  10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np8_Q16.SP.NMSE_OMP_sung)); ...
                  10*log10(mean(PS_Nt64_Nr16_Lt8_Lr4_Mult1p5_M128_Np10_Q16.SP.NMSE_OMP_sung))];

SNR_idx = 6;

figure
set(gcf, 'defaulttextinterpreter', 'latex')
hold on
set(gca, 'FontSize', 13)
set(gcf, 'Units', 'inches')
pos = get(gcf, 'position');
pos(3:4) = figureSize;
set(gcf, 'position', pos)

plot(2:2:10, NMSE_OMP_sung(:, SNR_idx), 'ko-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(2:2:10, NMSE_OMP_md(:, SNR_idx), 'kx-', 'LineWidth', 1.1, 'MarkerSize', 10)
plot(2:2:10, NMSE_OMP_lee(:, SNR_idx), 'ks-', 'LineWidth', 1.1, 'MarkerSize', 10)

hold off
grid on
xlabel('$N_p$')
ylabel('NMSE [dB]')
% title('Nt64 / Nr16 / Lt8 / Lr8 / G64 / Np6 / Q16')
hLeg = legend('Proposed', ...
       'Random', ...
       'MTC', ...
       'Location', 'northwest');
set(hLeg, 'Interpreter', 'latex')
set(hLeg, 'FontSize', 14)