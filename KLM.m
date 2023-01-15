function output = KLM(audio,fs)

%% 初始化参数
wlen = wlen_select(0.02,0.03,fs);
N = length(audio);
fn = fix(N/wlen);
S = zeros(wlen,fn);
for i = 1:fn
    S(:,i) = audio((i-1)*wlen+1:i*wlen);
end
tic
arOrder = 6; % AR模型阶数
numIter = 5; % 求解语音信号的AR模型参数时的迭代次数
H = [zeros(1, arOrder - 1), 1]; % 观测增益矩阵
R = var(audio(2*wlen:5*wlen)); % 噪声方差

errCov = R * eye(arOrder); % 后验估计误差协方差矩阵
output = zeros(length(audio),1); % 为输出信号分配内存
output(1:arOrder) = audio(1 : arOrder, 1)'; % 初始化输出信号，初始化为带噪语音样本前阶数个样本
estOutput = audio(1 : arOrder,1); % 初始化后验估计
Q = zeros(1,fn);
arCoeff = zeros(fn,arOrder+1);
for k = 1:fn % 对每一帧进行卡尔曼滤波
	oldOutput = estOutput; % 保存该帧迭代前的后验估计
    % 初始化开始进行卡尔曼滤波的位置
    if k == 1 
        iiStart = arOrder + 1; % 如果是第一帧，则从第arOrder+1个点开始处理
    end
%     [arCoeff(k,:), Q(k)] = (lpc(S(:,k), arOrder)); % 每帧带噪语音的AR模型系数arCoeff和方差Q
    [arCoeff(k,:), Q(k)] = lpc_copy(S(:,k),arOrder);
    % 迭代求解语音信号的AR模型参数
    for iter = 1 : numIter
		% 状态变换矩阵
        A = [zeros(arOrder - 1, 1), eye(arOrder - 1); fliplr(-arCoeff(k,2:end))]; % fliplr:左右翻转
        
        for ii = iiStart : wlen
			% 计算先验估计
			aheadEstOutput = A * estOutput;	
			% 计算先验估计误差的协方差矩阵p-
			aheadErrCov  = A * errCov * A' + H' * Q(k) * H;
			% 计算卡尔曼增益
			K = (aheadErrCov * H') / (H * aheadErrCov * H' + R);
			% 计算后验估计
			estOutput = aheadEstOutput + K * (S(ii,k) - H * aheadEstOutput);
			% 更新输出结果
			index = ii - iiStart + arOrder + 1 + (k - 1) * wlen;
			output(index - arOrder + 1 : index) = estOutput;
			% 计算后验估计误差的协方差矩阵p
			errCov  = (eye(arOrder) - K * H) * aheadErrCov;
        end
        iiStart = 1;
		if iter < numIter % 如果AR模型参数迭代还未完成，则恢复第一次迭代前的后验估计，用来进行下一次迭代
			estOutput = oldOutput; % 不然会有吱吱声
		end
		% 使用LPC算法计算上述滤波结果的AR模型系数以及方差
        [arCoeff(k , :), Q(k)] = lpc_copy(output((k - 1) * wlen + 1 : k * wlen), arOrder);
    end
end
toc
end
