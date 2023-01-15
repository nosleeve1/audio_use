function [use_frame,save_piece,last_frame,pos,zf] = fre_shift(frame,last_frame,save_piece,win,weight,pos,f_shift,fs,phsw,zf)
%%% frame:一列数据，pos:位置
%%% f_shift为移动频率,fs为采样频率 
buffer = [last_frame;frame];
len = single(length(buffer));
buffer = buffer.*win;
buffer_c = hand_hilbert(buffer);
n = [0:len-1]'+ pos * ones(len,1);
t = n*f_shift/fs;
if phsw
    delta = pos*f_shift/fs-0.5;
else
    delta = single(0);
end
shift_buffer = real(buffer_c.*exp(1j*2*pi*(t+sin(delta*pi)))).*win;
use_frame = (shift_buffer(1:end/2) + save_piece)./weight;
save_piece = shift_buffer(end/2+1:end);
last_frame = frame;

if pos+len-1 > abs(fs/f_shift)
    pos = fix(pos - abs(fs/f_shift));
end
pos = pos + len/2;

end

function ha = hand_hilbert(a)
len = single(length(a));
A = fft(a);
A(1) = 0;
A(2:ceil(len/2)) = A(2:ceil(len/2)) * -1i;
A(ceil(len/2)+1:len) = A(ceil(len/2)+1:len) * 1i;
if mod(len,2) == 0
    A(len/2+1) = 0;
end
ha = a + ifft2(A) * 1i;
end
