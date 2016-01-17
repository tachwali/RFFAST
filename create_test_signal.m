N=49*50*51;
t=0:N-1;
x= exp(2*pi*j*60000*t/N);%exp(2*pi*j*800*t/N) + exp(2*pi*j*300*t/N) + exp(2*pi*j*900*t/N);
%x=awgn(x,20);