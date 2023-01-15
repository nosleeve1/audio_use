[x,fs] = audioread('test.wav');
x = x(:,1);
% x = x + 0.5*sin(2*pi*1500*[1:len]'/fs);%+1*sin(2*pi*2500*[1:len]'/fs).*y;
% x = sin(2*pi*500*[1:128*4]'/fs);

f_shift = 0;
pos = 0;
n = 128;
s = [];
x = single(x);
len = 2*n;
win = sin(pi/(len-1)*[0:len-1]');
weight = (win(1:len/2).^2+win(len/2+1:end).^2);
save_piece = single(zeros(n,1));
last_frame = single(zeros(n,1));
f_shift = single(f_shift);
fs = single(fs);
pos = single(pos);
n = single(n);
s = single(s);
win = single(win);
weigth = single(weight);
tic
phsw = single(0);
p = 10;
ar = zeros(fix(length(x)/n),p+1);
zf = zeros(16,1);
for i = n:n:length(x)
    frame = x(i-n+1:i);
    [use_frame,save_piece,last_frame,pos,zf] = fre_shift(frame,last_frame,save_piece,win,weight,pos,f_shift,fs,phsw,zf);
    use_frame = use_frame(:);    
    s = [s;use_frame];
end
toc
subplot(211);
plot(x,'b');
ylim([-1,1]);
subplot(212);
plot(s,'r');
ylim([-1,1]);
audiowrite('x.wav',x,fs);
audiowrite('s.wav',s,fs);


% dlmwrite('input.txt', x, 'delimiter' , '\t' , 'precision', '%0.7f');
% load('output.txt');
% sound(output,fs);
% audiowrite('output.wav',output,fs);
