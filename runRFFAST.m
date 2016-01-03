clear all; close all;
create_test_signal;
subplot(2,1,1);
plot(abs(fft(x)));
title('fft(x)');
subplot(2,1,2);
plot(abs(RFFAST(x)));
title('RFFAST(x)');

