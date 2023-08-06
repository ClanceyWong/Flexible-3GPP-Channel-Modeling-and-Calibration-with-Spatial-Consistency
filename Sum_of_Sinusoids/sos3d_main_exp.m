clear
show_progress = 0;
set(0,'defaultAxesFontName','Times') % Default Font Type
set(0,'DefaultAxesFontSize',12);     % Default Font Size
set(0,'defaultTextFontName','Times')
set(0,'DefaultTextFontSize',12);
d_cor  = 10;                  % 自相关距�?
ACFexp = @(d)exp(-(d)/d_cor); % 自相关函数ACF

% train
prop = (1-exp(-50/d_cor))/200;
k = 0:200;                         % 训练集不同距离的数量
dist_train = -d_cor*log(1-prop*k); % 训练集的距离
directions_train = 56;             % 训练集的每个距离的方向数�?
[theta, phi] = fibonacci_sphere(directions_train);
dx = cos( phi ) .* cos( theta )*dist_train;  % 转换为直角坐标系的位�?
dy = sin( phi ) .* cos( theta )*dist_train;
dz = sin( theta )*dist_train;
xyz(1,:,:) = dx';
xyz(2,:,:) = dy';
xyz(3,:,:) = dz';
xyz = reshape(xyz,3,[])';
acf = repmat(ACFexp(dist_train),1,directions_train)';

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
xyz_v = reshape(xyz_v,3,[])';
acf_v = repmat(ACFexp(dist_valid),1,direction_valid)';

times = 10; % 重复次数
% num_sinusoids = (100:100:1000); % 正弦波数�?
num_sinusoids = 500; % 正弦波数�?
error_log = zeros(times,length(num_sinusoids));
error4plot = nan(200,1000);
rng(769);

for time = 1:times
    fprintf(['Time ',num2str(time),', 10log10(error_test): '])
    for N = num_sinusoids
        SOS_model = sos_create([3,N,1],'LearningRate',0.02,'OptimizationMethod','Adam','WeightDecay',0.00001); % 'normal' or 'Adam'
        %training
        option.batch_size = 201;  % 批大�?
        option.iteration = 1;
        rep = 1000; % epoches
        min_error_valid = inf; % 初始�?
        totalCost = zeros(1,rep);
        totalError = zeros(1,rep);
        iteration = 0;
        breakflag = 0;
        while iteration <= rep
            iteration = iteration +1; % �? iteration 次迭�?
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
                min_error_valid = error_valid; % 更新�?小验证误�?
                min_error_train = error_train; % �?小验证误差对应的训练误差
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
        xyzt = reshape(xyzt,3,[])';
        acft = repmat(ACFexp(dist_test),1,direction_test)';
        [error_test,storedSOS_model] = sos_test(storedSOS_model,xyzt,acft);
        error_log(time,N/100) = 10*log10(error_test);
        fprintf([num2str(error_log(time,N/100)),' dB, ']); % '10log10(error_test) = ',
    end
    fprintf('\n');
end

if N == 500
    plot(1:numel(nanmean(error4plot)), nanmean(error4plot)); hold on;
end

error_log_3d_exp_new = error_log;
save error_log_3d_exp_new error_log_3d_exp_new;
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % figure
% % scatter3(storedNN.W{1}(:,1),storedNN.W{1}(:,2),storedNN.W{1}(:,3))
% figure(8);hold off;
% plot(dist_train(1:4:end),ACFexp(dist_train(1:4:end)),'Linewidth',2 );hold on;grid on;
% rrr = reshape(storedSOS_model.a{3}',numel(dist_train),[])';
% meanR = mean(rrr); maxR = max(rrr)-mean(rrr); minR = min(rrr)-mean(rrr);
% plot(dist_train,meanR,'r-','Linewidth',2);
% errorbar(dist_train(1:8:end),meanR(1:8:end),maxR(1:8:end),minR(1:8:end),'.r-');
% legend('ideal ACF','ave. estimated ACF','range of estimated ACF','location','northeast');
% xlabel('distance [m]'); ylabel('ACF');
% drawnow
% 
% %% Figure 1. Desired 2-D ACF on the x ? y plane (top) and approximated ACF using 300 sinusoids (bottom)
% clear qxyz
% xy      = -50:1:50;
% zz      = 0;
% [x,y,z] = meshgrid(xy,xy,zz);
% qxyz(1,:,:) = x;
% qxyz(2,:,:) = y;
% qxyz(3,:,:) = z;
% qxyz = reshape(qxyz,3,[])';
% dist_train = sqrt(sum(qxyz.^2,2));
% qacf          = ACFexp(dist_train)';
% acf_d  = reshape(qacf,size(x,1),size(x,2));
% [error_train,storedSOS_model] = sos_test(storedSOS_model,qxyz,qacf);
% acf_e = reshape(storedSOS_model.a{3},size(x,1),size(x,2));
% figure(9); 
% set(gcf,'Position',[200 150 1000 500]);
% subplot(121);hold off;
% mesh(x,y,acf_d);
% xlim([min(xy),max(xy)]); ylim([min(xy),max(xy)]); zlim([min(zz),1]); 
% xlabel('\Delta{\itx} [m]');ylabel('\Delta{\ity} [m]');zlabel('{\it\rho}(\Delta{\itx},\Delta{\ity},0)'); view([60, 20]);
% hold on; grid on;drawnow
% subplot(122);hold off;
% mesh(x,y,acf_e);
% xlim([min(xy),max(xy)]); ylim([min(xy),max(xy)]); zlim([min(acf_e(:)),1]); 
% xlabel('\Delta{\itx} [m]');ylabel('\Delta{\ity} [m]');zlabel("{\itf}(\Delta{\itx},\Delta{\ity},0)"); view([60, 20]);
% hold on; grid on;drawnow
% % mesh(x,y,z);
% 
% %% Figure 2. Example of a spatially correlated 2-D process
% clear pxyz
% xy      = 0:100;
% zz      = 0:100;
% [x,y,z] = meshgrid(xy,xy,zz);
% pxyz(1,:,:,:) = x;
% pxyz(2,:,:,:) = y;
% pxyz(3,:,:,:) = z;
% pxyz = reshape(pxyz,3,[]);
% N = size(storedSOS_model.W{1},1);
% phase = (2*pi*rand(N,1)-pi);
% % phase3 = [sqrt(0.5)*randn(10*N/2,1)-pi/2 ; sqrt(0.5)*randn(10*N/2,1)+pi/2];
% % % phase3 = phase3(phase3<=pi&phase3>=-pi);
% % phase3 = sort(phase3);
% % nu     = numel(phase3);
% % stp    = floor(nu/N);
% % phase3 = phase3(1:stp:end);
% % phase3 = phase3(randperm(N));
% % figure(9);ecdf(phase3);hold on; ecdf(phase)
% 
% shawing = sum(sqrt(2*storedSOS_model.W{2}').*cos(storedSOS_model.W{1}*pxyz + phase));
% 
% % plot
% shawing      = reshape(shawing,size(x,1),size(x,2),size(x,3));
% shawingplot  = shawing;
% shawingplot(shawingplot<-3.5) = -3.5;
% shawingplot(shawingplot> 3.5) =  3.5;
% figure(10);hold off;
% slice(x(:,:,100:101),y(:,:,100:101),z(:,:,100:101),shawingplot(:,:,100:101),xy,xy,[100]); shading interp;hold on;
% slice(x(100:101,:,:),y(100:101,:,:),z(100:101,:,:),shawingplot(100:101,:,:),[100],xy,zz); shading interp;hold on;
% slice(x(:,1:2,:),y(:,1:2,:),z(:,1:2,:),shawingplot(:,1:2,:),xy,[0],zz); shading interp;hold on;
% xlim([min(xy),max(xy)]); ylim([min(xy),max(xy)]); zlim([min(zz),max(zz)]); 
% xlabel('x-position (m)');ylabel('y-position (m)');zlabel('z-position (m)'); view([0, 90]);
% hold on; box on; colorbar; colormap('jet');drawnow
% 
% 
% %% Figure 4. Estimated ACF from spatially correlated random values
% clear sxyz
% xy      = 0:0.25:50;
% [x,y,z] = meshgrid(xy,xy,xy);
% % xy      = 1:1000;
% % zz      = 1:50;
% % [x,y,z] = meshgrid(xy,xy,zz);
% sxyz(1,:,:,:) = x;
% sxyz(2,:,:,:) = y;
% sxyz(3,:,:,:) = z;
% sxyz = reshape(sxyz,3,[]);
% dist_group = 1:2:70;
% group_no   = numel(dist_group);
% num_pos    = 10000;
% 
% timetotal = 100;
% pearson        = zeros(timetotal,group_no);
% for times = 1:timetotal
%     phase   = (2*pi*rand(N,1)-pi);
%     for g = 1:group_no
%         sxyz1 = sxyz(:,randi(size(sxyz,2),1,num_pos));
%         shawing1 = sum(sqrt(2*storedSOS_model.W{2}').*cos(storedSOS_model.W{1}*sxyz1 + phase));
%         
%         theta1 = rand(1,num_pos)*pi;
%         phi1   = rand(1,num_pos)*2*pi;
%         direc1 = [sin(theta1).*cos(phi1);
%                   sin(theta1).*sin(phi1);
%                   cos(theta1)]; 
%         sxyz2 = sxyz1 + dist_group(g)*direc1;
%         shawing2 = sum(sqrt(2*storedSOS_model.W{2}').*cos(storedSOS_model.W{1}*sxyz2 + phase));
% %         pearson(times,g) = 1 - sum((shawing1-shawing2).^2)/2/num_pos;
%         pearson(times,g) = sum((shawing1-mean(shawing1)).*(shawing2-mean(shawing2)))/sqrt(var(shawing1)*var(shawing2))/num_pos;
% 
%     end
% end
% pearmean = mean(pearson);
% pearmax  = max(pearson-mean(pearson));
% pearmin  = min(pearson-mean(pearson));
% figure(11);plot(0:70,ACFexp(0:70),'Linewidth',1.5);hold on ;
% plot(dist_group,pearmean,'s-','Linewidth',1.5);grid on;
% errorbar(dist_group,pearmean,pearmax,pearmax,'r.','Linewidth',1);drawnow
% pearvar  = sqrt(var(pearson));
% errorbar(dist_group,pearmean,pearvar,'.','Linewidth',1);
% legend('ideal ACF','estimated ACF','[min., max.]','[{\it-\sigma}, {\it\sigma}]','location','northeast');
% xlabel('distance [m]'); ylabel('ACF');
% 
% 
% %% Figure 5. CDF of the output values vs. the Gaussian CDF
% clear sxyz
% xy      = 0:1000;
% zz      = 0:50;
% [x,y,z] = meshgrid(xy,xy,zz);
% sxyz(1,:,:,:) = x;
% sxyz(2,:,:,:) = y;
% sxyz(3,:,:,:) = z;
% sxyz = reshape(sxyz,3,[]);
% sxyz = sxyz(:,randperm(size(sxyz,2)));
% sxyz = sxyz(:,1:10000);
% phase   = (2*pi*rand(N,1)-pi);
% shawing0 = sum(sqrt(2*storedSOS_model.W{2}').*cos(storedSOS_model.W{1}*sxyz + phase));
% rn      = randn(1,numel(shawing0));
% prob    = (1:numel(rn))/numel(rn);
% rn      = sort(rn); shawing0 = sort(shawing0(:));
% ind = [];
% for p = -3:0.2:3
%     ind = [ind,find(abs(shawing0-p)==min(abs((shawing0-p))))];
% end
% figure(12);hold off; plot(rn(ind),prob(ind),'-','Linewidth',2);hold on;
% plot(shawing0(ind),prob(ind),'s-','Linewidth',1);hold on;grid on;
% xlim([-3,3]); xlabel('random number'); ylabel('empirical CDF');
% legend('Gaussian CDF','3-D Apprximation','location','southeast');drawnow
% 
% 
% %%
% clear sxyz
% N = size(storedSOS_model.W{1},1);
% xy      = 0:1000;
% zz      = 0:50;
% [x,y,z] = meshgrid(xy,xy,zz);
% sxyz(1,:,:,:) = x;
% sxyz(2,:,:,:) = y;
% sxyz(3,:,:,:) = z;
% sxyz = reshape(sxyz,3,[]);
% sxyz = sxyz(:,randperm(size(sxyz,2)));
% 
% total_t = 10;
% shawing0 = zeros(total_t,10000);
% for time_rn = 1:total_t
%     clear sxyz_
%     sxyz_ = sxyz(:,randi(size(sxyz,2),1,10000));
% 
%     phase   = (2*pi*rand(N,1)-pi);
%     shawing0(time_rn,:) = sum(sqrt(2*storedSOS_model.W{2}').*cos(storedSOS_model.W{1}*sxyz_ + phase));
%     shawing0(time_rn,:) = sort(shawing0(time_rn,:));
% end
% rn      = randn(1,size(shawing0,2));
% prob    = (1:numel(rn))/numel(rn);
% rn      = sort(rn);
% 
% figure(12);hold off; plot(rn(1:10:end),prob(1:10:end),'k-','Linewidth',2);hold on;
% for iii = 1:10
% plot(shawing0(iii,:),prob,'r--','Linewidth',1);hold on;grid on;
% end
% plot(rn(1:10:end),prob(1:10:end),'k-','Linewidth',2);hold on;
% % plot(max(shawing0(ind)),prob(ind(1,:)),'y--','Linewidth',1);
% % plot(min(shawing0(ind)),prob(ind(1,:)),'y--','Linewidth',1);
% xlim([-3,3]); xlabel('random number'); ylabel('empirical CDF');
% legend('Gaussian CDF','3-D Apprximation','location','southeast');drawnow