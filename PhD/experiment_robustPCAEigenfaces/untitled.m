clear; clc; close all;

Tpulse = 20*exp(-3);
Fs = 10*exp(3);

t = -1:1/Fs:1;

x = rectpuls(t,Tpulse);

plot(t,x)

rng default

y = 0.1*randn(size(x));

s = x + y;

hold on
plot(t,s)

snr(s,y)
