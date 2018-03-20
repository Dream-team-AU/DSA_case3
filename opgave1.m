%% Midlingsfiltre
%  HTR 20/03-2018

%% Generelt setup:
clear; close all; clc; format compact
Nfft = 2048;
k = 0:Nfft-1;
N = 2500;          %    <-- prøv N = 1e5, hvis teorien skal passe...
n = 0:N-1;


%% Indlæsning af data, samt skabelse af 2 dele af disse

load('vejecelle_data.mat');
x = vejecelle_data;
x1 = vejecelle_data(1:1000);
n1 = (1:1000);
N1 = 1000;
x2 = vejecelle_data(1050:2500);
n2 = (1050:2500);
N2 = 1450;

figure
subplot(1,2,1:2)
plot(n,x), grid
xlabel('n'), ylabel('x(n)'), title('hele signalet DC-signal')

%% MA-filter (ikke-rekursivt) for første del
M = 16;    % filterkoefficienter
hMA = 1/M*ones(1,M); % MA-filter, filterkoefficienter

hMA_imp_resp  = hMA;                        % impulsrespons
hMA_step_resp = filter(hMA,1,ones(1,2*M));  % steprespons
L_MA_trans_resp = M-1;                      % længde af transientrespons
HMA1 = fft(hMA,Nfft);                        % frekvensrespons

yMA1 = filter(hMA,1,x1);            % filtrerer inputsignal for første del
yMA2 = filter(hMA,1,x2);            % filtrerer inputsignal for anden del

var_x1 = var(x1(M:N1));    % varians i 1 signal i del efter transientrespons
var_yMA1 = var(yMA1(M:N1));% varians i 1 signal i del efter transientrespons

var_x2 = var(x2(M:N2));    % varians i 2 signal i del efter transientrespons
var_yMA2 = var(yMA2(M:N2));% varians i 2 signal i del efter transientrespons

%% --- plotting for første signal ---
figure('name', 'første MA-filter')
subplot(2,6,1:3)
plot(n1,x1), grid
xlabel('n'), ylabel('x(n)'), title('første del af DC-signal')

subplot(2,6,4:6)
histogram(x1)
title('histogram')

subplot(2,6,7:10)
plot(n1,x1), grid, hold on
plot(n1,yMA1,'linewidth',2)
xlabel('n'), ylabel('x(n)'), title(['MA-filter, M = ' num2str(M)])
legend('input','output')

subplot(2,6,11:12)
text(0,0.5,...
    {['første MA-filter, M = ' num2str(M)],...
     ['Transientrespons: ' num2str(L_MA_trans_resp) ' samples'],...
     ['Støjeffekt i inputsignal (varians):  ' num2str(var_x1)],...
     ['Støjeffekt i outputsignal (varians): ' num2str(var_yMA1)],...
     ['Reduktion i støjeffekt: ' num2str((var_x1/var_yMA1)) ' gange'],...
     ['Reduktion i støjeffekt: ' num2str(10*log10(var_x1/var_yMA1)) ' dB']})
 axis off

%% --- plotting for andet signal ---
figure('name', 'andet MA-filter')
subplot(2,6,1:3)
plot(n2,x2), grid
xlabel('n'), ylabel('x(n)'), title('anden del af DC-signal')

subplot(2,6,4:6)
histogram(x2)
title('histogram')

subplot(2,6,7:10)
plot(n2,x2), grid, hold on
plot(n2,yMA2,'linewidth',2)
xlabel('n'), ylabel('x(n)'), title(['MA-filter, M = ' num2str(M)])
legend('input','output')

subplot(2,6,11:12)
text(0,0.5,...
    {['første MA-filter, M = ' num2str(M)],...
     ['Transientrespons: ' num2str(L_MA_trans_resp) ' samples'],...
     ['Støjeffekt i inputsignal (varians):  ' num2str(var_x2)],...
     ['Støjeffekt i outputsignal (varians): ' num2str(var_yMA2)],...
     ['Reduktion i støjeffekt: ' num2str((var_x2/var_yMA2)) ' gange'],...
     ['Reduktion i støjeffekt: ' num2str(10*log10(var_x2/var_yMA2)) ' dB']})
 axis off
