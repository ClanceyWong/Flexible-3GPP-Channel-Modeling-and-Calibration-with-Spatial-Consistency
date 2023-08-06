clear
show_progress = 0;
set(0,'defaultAxesFontName','Times') % Default Font Type
set(0,'DefaultAxesFontSize',12);     % Default Font Size
set(0,'defaultTextFontName','Times')
set(0,'DefaultTextFontSize',12);
d_cor  = 10;                  % 自相关距离
ACFexp = @(d)exp(-(d)/d_cor); % 自相关函数ACF

% train
prop = (1-exp(-50/d_cor))/200;
k = 0:200;                         % 训练集不同距离的数量
dist_train = -d_cor*log(1-prop*k); % 训练集的距离
directions_train = 112;             % 训练集的每个距离的方向数量
[theta, phi] = fibonacci_sphere(directions_train);
dx = cos( phi ) .* cos( theta )*dist_train;  % 转换为直角坐标系的位置
dy = sin( phi ) .* cos( theta )*dist_train;
dz = sin( theta )*dist_train;
xyz(1,:,:) = dx';
xyz(2,:,:) = dy';
xyz(3,:,:) = dz';
xyz = reshape(xyz,3,[]);
rind = randperm(size(xyz,2));
acf = repmat(ACFexp(dist_train),1,directions_train);
xyz = xyz(:,rind);
ddd = sqrt(sum(xyz.^2));
dt = ddd(1:2:end);
dr = ddd(2:2:end);
xyz = reshape(xyz,6,[])';
acf = acf(rind);
acf = reshape(acf,2,[])';
acf = acf(:,1).*acf(:,2);

% valid
dist_valid = 1:50;
direction_valid = 80;
[theta_v, phi_v] = fibonacci_sphere(direction_valid);
dx_v = cos( phi_v ) .* cos( theta_v )*dist_valid;
dy_v = sin( phi_v ) .* cos( theta_v )*dist_valid;
dz_v = sin( theta_v )*dist_valid;
xyz_v(1,:,:) = dx_v';
xyz_v(2,:,:) = dy_v';
xyz_v(3,:,:) = dz_v';
xyz_v = reshape(xyz_v,3,[]);
rind = randperm(size(xyz_v,2));
acf_v  = repmat(ACFexp(dist_valid),1,direction_valid);
xyz_v          = xyz_v(:,rind);
xyz_v          = reshape(xyz_v,6,[])';
acf_v          = acf_v(rind);
acf_v          = reshape(acf_v,2,[])';
acf_v          = acf_v(:,1).*acf_v(:,2);

times = 200; % 重复次数
% num_sinusoids = (100:100:1000); % 正弦波数量
num_sinusoids = 500;
error_log = zeros(times,length(num_sinusoids));
error4plot = nan(200,1000);
rng(769);

for time = 1:times
    fprintf(['Time ',num2str(time),', 10log10(error_test): '])
    for N = num_sinusoids
        SOS_model = sos_create([6,N,1],'LearningRate',0.05,'OptimizationMethod','Adam','WeightDecay',0.00001); % 'normal' or 'Adam'
        %training
        option.batch_size = 201;  % 批大小
        option.iteration = 1;
        rep = 1000; % epoches
        min_error_valid = inf; % 初始化
        totalCost = zeros(1,rep);
        totalError = zeros(1,rep);
        iteration = 0;
        breakflag = 0;
        while iteration <= rep
            iteration = iteration +1; % 第 iteration 次迭代
            SOS_model = sos_train(SOS_model,option,xyz,acf); % 训练模型
            
            totalCost(iteration) = sum(SOS_model.cost)/length(SOS_model.cost); % 计算代价
            cost = totalCost(iteration);
            [error_train, SOS_model] = sos_test(SOS_model, xyz, acf); % 在训练集上的误差
            if N == 500
                error4plot(time, iteration) = error_train;
            end
            
            model_forValid = SOS_model;
            [error_valid, model_forValid] = sos_test(model_forValid, xyz_v, acf_v); % 在验证集上的误差
            totalError(iteration) = error_train;
            breakflag = breakflag + 1;
            if error_valid < min_error_valid
                min_error_valid = error_valid; % 更新最小验证误差
                min_error_train = error_train; % 最小验证误差对应的训练误差
                storedSOS_model = SOS_model; % 保存SOS模型
                it = iteration; % 当前迭代次数
                breakflag = 0;
            end
            if breakflag >= 20 % 连续20验证集误差未减小
                break;
            end
        end
        % fprintf(['Time: ',num2str(time),', the number of sinusoids is ',num2str(N),':\n','    valid error: ',num2str(10*log10(min_error_valid)),' dB;    cost: ',num2str(cost),'.\n']);
        
        %% testing
        clear xyzt
        direction_test = 200;
        [theta, phi] = fibonacci_sphere(direction_test);
        dist_test = 1:1:50;
        dxt = cos( phi ) .* cos( theta )*dist_test;
        dyt = sin( phi ) .* cos( theta )*dist_test;
        dzt = sin( theta )*dist_test;
        xyzt(1,:,:) = dxt';
        xyzt(2,:,:) = dyt';
        xyzt(3,:,:) = dzt';
        xyzt1           = reshape(xyzt,3,[])';
        rind1          = randperm(size(xyzt1,1));
        xyzt1          = xyzt1(rind1,:);
        dist1          = sqrt(sum(xyzt1.*xyzt1,2));
        acft1          = ACFexp(dist1);
        xyzt2           = reshape(xyzt,3,[])';
        rind2          = randperm(size(xyzt2,1));
        xyzt2          = xyzt2(rind2,:);
        dist2          = sqrt(sum(xyzt2.*xyzt2,2));
        acft2          = ACFexp(dist2);
        
        [error_test,storedSOS_model] = sos_test(storedSOS_model,[xyzt1,xyzt2],(acft1.*acft2));
        error_log(time,N/100) = 10*log10(error_test);
        fprintf([num2str(error_log(time,N/100)),' dB, ']); % '10log10(error_test) = ',
    end
    fprintf('\n');
end

if N == 500
    plot(1:numel(nanmean(error4plot)), nanmean(error4plot)); hold on;
end

error_log_6d_exp = error_log;
save error_log_6d_exp error_log_6d_exp;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%