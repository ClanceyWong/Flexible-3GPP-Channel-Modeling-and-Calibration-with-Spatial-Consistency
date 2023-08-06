function SOS_model = sos_backpropagation(SOS_model,batch_acf)
    batch_acf = batch_acf';
    m = size(SOS_model.a{1},2);  % m 为 batch 的数目
    SOS_model.theta{1} = 0;      % 用来记录每一层的导数大小，这也是后向传播的意义
    switch SOS_model.output_function
        case 'sigmoid'
            SOS_model.theta{SOS_model.depth} = -(batch_acf-SOS_model.a{SOS_model.depth}) .* SOS_model.a{SOS_model.depth} .* (1 - SOS_model.a{SOS_model.depth});
        case 'tanh'
            SOS_model.theta{SOS_model.depth} = -(batch_acf-SOS_model.a{SOS_model.depth}) .* (1 - SOS_model.a{SOS_model.depth}.^2);
        case 'softmax'
            SOS_model.theta{SOS_model.depth} = SOS_model.a{SOS_model.depth} - batch_acf;
        case 'sum'
            SOS_model.theta{SOS_model.depth} = SOS_model.a{SOS_model.depth} - batch_acf;
    end

    % 第二项是引入正则项后求导产生的，即正则项系数 SOS_model.weight_decay 乘以 w
    SOS_model.W_grad{SOS_model.depth-1} = SOS_model.theta{SOS_model.depth}*SOS_model.a{SOS_model.depth-1}'/m + SOS_model.weight_decay*SOS_model.W{SOS_model.depth-1};
    % 上面式子都对 m 个梯度求了平均
    switch SOS_model.active_function
        case 'cos'
            for ll = 2 : SOS_model.depth - 1
                k = SOS_model.depth - ll + 1; % 这里就用到了后向传播的典型意义，从后向前递推
                SOS_model.theta{k} = ((SOS_model.W{k}'*SOS_model.theta{k+1})) .* SOS_model.ad{k};
                SOS_model.W_grad{k-1} = SOS_model.theta{k}*SOS_model.a{k-1}'/m + SOS_model.weight_decay*SOS_model.W{k-1};
            end
    end
end