function SOS_model = sos_train(SOS_model,option,train_xyz,train_acf)
    iteration = option.iteration;     % 迭代次数，即操作几轮
    batch_size = option.batch_size;   % 累积BP算法中的批量大小
    m = size(train_xyz,1);            % size(，1)返回矩阵行数
    num_batches = m / batch_size;
    for k = 1 : iteration
        index = randperm(m);             % 打乱1到m的顺序 
        if  strcmp(SOS_model.optimization_method,'Adam')
            SOS_model.AdamTime = SOS_model.AdamTime+1;
        end
        for batche = 1 : num_batches       % 即对60000个训练数据打乱后打成小包
            batch_xyz = train_xyz(index((batche-1)*batch_size + (1:batch_size)),:);
            batch_acf = train_acf(index((batche-1)*batch_size + (1:batch_size)),:);
            SOS_model = sos_forward(SOS_model,batch_xyz,batch_acf); % 前向计算
            SOS_model = sos_backpropagation(SOS_model,batch_acf);   % 后向计算
            SOS_model = sos_applygradient(SOS_model);               % 梯度下降
        end
    end
end