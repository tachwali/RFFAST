create_test_signal;
subplot(2,1,1);
title('fft(x)')
plot(abs(fft(x)));
subplot(2,1,2);
title('RFFAST(x)')
plot(abs(RFFAST(x)))
