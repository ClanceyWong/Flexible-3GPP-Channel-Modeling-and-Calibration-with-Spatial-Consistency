function SOS_model = sos_forward(SOS_model,batch_xyz,batch_acf)    
    s = size(SOS_model.cost) + 1;  % s为cost矩阵的行和列 2维向量并各加1,这步配合第51行的 SOS_model.cost(s)
                                   % 实际效果其实就是每次在cost行向量后挤入一个新值
    batch_xyz = batch_xyz';
    batch_acf = batch_acf';
    m = size(batch_xyz,2);    % size（，2）得到矩阵列数
    SOS_model.a{1} = batch_xyz;
    cost2 = 0;              % cost2指cost的第二个和式即添加的正则项
    for k = 2 : SOS_model.depth
        y = SOS_model.W{k-1} * SOS_model.a{k-1};
        if k == SOS_model.depth % 输出层激活函数选择
            switch SOS_model.output_function
                case 'sigmoid'
                    SOS_model.a{k} = sigmoid(y);
                case 'tanh'
                    SOS_model.a{k} = tanh(y);
                case 'relu'
                    SOS_model.a{k} = max(y,0);
                case 'softmax'
                    SOS_model.a{k} = softmax(y);
                case 'sum'
                    SOS_model.a{k} = y;
            end
        else 
            switch SOS_model.active_function  % 隐层激活函数选择
                case 'sigmoid'
                    SOS_model.a{k} = sigmoid(y);
                case 'tanh'
                    SOS_model.a{k} = tanh(y);
                case 'relu'
                    SOS_model.a{k} = max(y,0);
                case 'cos'
                    SOS_model.a{k}  = cos(y);
                    SOS_model.ad{k} = -sin(y);
            end
        end
        cost2 = cost2 + sum(sum(SOS_model.W{k-1}.^2)); % 正则项计算
    end
    if strcmp(SOS_model.objective_function,'MSE')
        SOS_model.cost(s) = 0.5 / m * sum(sum((SOS_model.a{k} - batch_acf).^2)) + 0.5 * SOS_model.weight_decay * cost2;
    elseif strcmp(SOS_model.objective_function,'Cross Entropy')
        SOS_model.cost(s) = -0.5*sum(sum(batch_acf.*log(SOS_model.a{k})))/m + 0.5 * SOS_model.weight_decay * cost2;
    end
end